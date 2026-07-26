/// <reference lib="webworker" />

import { GravityStore, TerrainStore } from "./data-store";
import type {
  ErrorMessage,
  LaunchInput,
  ProgressMessage,
  ResultMessage,
  RunMessage,
  SolverConfig,
  SolverResult,
  TrajectoryPoint,
} from "./types";

const ctx: DedicatedWorkerGlobalScope = self as unknown as DedicatedWorkerGlobalScope;

const MU_MOON = 4.902800118e12;
const MOON_RADIUS_M = 1_737_400;
const MOON_ROTATION_PERIOD_S = 2_360_591.5104;
const OMEGA = (2 * Math.PI) / MOON_ROTATION_PERIOD_S;
const BALL_RADIUS_M = 0.036888;
const MU_EARTH = 3.986004418e14;
const MU_SUN = 1.32712440018e20;
const EARTH_DISTANCE_M = 384_400_000;
const SUN_DISTANCE_M = 149_597_870_700;

type Vec3 = [number, number, number];
type State = [number, number, number, number, number, number];

const add3 = (a: Vec3, b: Vec3): Vec3 => [
  a[0] + b[0],
  a[1] + b[1],
  a[2] + b[2],
];
const subtract3 = (a: Vec3, b: Vec3): Vec3 => [
  a[0] - b[0],
  a[1] - b[1],
  a[2] - b[2],
];
const scale3 = (a: Vec3, factor: number): Vec3 => [
  a[0] * factor,
  a[1] * factor,
  a[2] * factor,
];
const norm3 = (a: Vec3): number => Math.hypot(a[0], a[1], a[2]);
const distance3 = (a: Vec3, b: Vec3): number => norm3(subtract3(a, b));
const cross3 = (a: Vec3, b: Vec3): Vec3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
const degrees = (value: number): number => (value * 180) / Math.PI;
const radians = (value: number): number => (value * Math.PI) / 180;

const geodetic = (position: Vec3): [number, number, number] => {
  const radius = norm3(position);
  const latitude = degrees(Math.asin(position[2] / radius));
  const longitude = ((degrees(Math.atan2(position[1], position[0])) % 360) + 360) % 360;
  return [latitude, longitude, radius];
};

const localBasis = (latitudeDeg: number, longitudeDeg: number): [Vec3, Vec3, Vec3] => {
  const latitude = radians(latitudeDeg);
  const longitude = radians(longitudeDeg);
  const cl = Math.cos(latitude);
  const sl = Math.sin(latitude);
  const co = Math.cos(longitude);
  const so = Math.sin(longitude);
  const radial: Vec3 = [cl * co, cl * so, sl];
  const east: Vec3 = [-so, co, 0];
  const north: Vec3 = [-sl * co, -sl * so, cl];
  return [radial, east, north];
};

const rotateZ = (vector: Vec3, angle: number): Vec3 => {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  return [c * vector[0] - s * vector[1], s * vector[0] + c * vector[1], vector[2]];
};

const thirdBody = (position: Vec3, bodyPosition: Vec3, mu: number): Vec3 => {
  const relative = subtract3(bodyPosition, position);
  const relativeNorm = norm3(relative);
  const bodyNorm = norm3(bodyPosition);
  return subtract3(
    scale3(relative, mu / relativeNorm ** 3),
    scale3(bodyPosition, mu / bodyNorm ** 3),
  );
};

const initialState = (
  launch: LaunchInput,
  launchRadiusM: number,
): { state: State; site: Vec3 } => {
  const [radial, east, north] = localBasis(
    launch.latitudeDeg,
    launch.longitudeDegEast,
  );
  const azimuth = radians(launch.azimuthDeg);
  const elevation = radians(launch.elevationDeg);
  const horizontal = add3(scale3(north, Math.cos(azimuth)), scale3(east, Math.sin(azimuth)));
  const velocity = scale3(
    add3(scale3(horizontal, Math.cos(elevation)), scale3(radial, Math.sin(elevation))),
    launch.speedMS,
  );
  const site = scale3(radial, launchRadiusM);
  return {
    state: [site[0], site[1], site[2], velocity[0], velocity[1], velocity[2]],
    site,
  };
};

const derivative = async (
  timeS: number,
  state: State,
  gravity: GravityStore,
  config: SolverConfig,
): Promise<State> => {
  const position: Vec3 = [state[0], state[1], state[2]];
  const velocity: Vec3 = [state[3], state[4], state[5]];
  const [latitude, longitude, radius] = geodetic(position);
  const [gr, gtheta, gphi] = await gravity.sphericalAcceleration(
    radius,
    latitude,
    longitude,
  );
  const [radial, east, north] = localBasis(latitude, longitude);
  let acceleration = add3(
    add3(scale3(radial, gr), scale3(north, -gtheta)),
    scale3(east, gphi),
  );

  const omegaVector: Vec3 = [0, 0, OMEGA];
  acceleration = add3(acceleration, scale3(cross3(omegaVector, velocity), -2));
  acceleration = add3(
    acceleration,
    scale3(cross3(omegaVector, cross3(omegaVector, position)), -1),
  );
  if (config.includeEarthTide) {
    const earth = rotateZ([EARTH_DISTANCE_M, 0, 0], -OMEGA * timeS);
    acceleration = add3(acceleration, thirdBody(position, earth, MU_EARTH));
  }
  if (config.includeSolarTide) {
    const sun = rotateZ([0, SUN_DISTANCE_M, 0], -OMEGA * timeS);
    acceleration = add3(acceleration, thirdBody(position, sun, MU_SUN));
  }
  return [
    velocity[0],
    velocity[1],
    velocity[2],
    acceleration[0],
    acceleration[1],
    acceleration[2],
  ];
};

const combineState = (state: State, derivatives: State[], weights: number[], step: number): State =>
  state.map(
    (value, index) =>
      value +
      step *
        derivatives.reduce(
          (sum, derivativeValue, derivativeIndex) =>
            sum + weights[derivativeIndex] * derivativeValue[index],
          0,
        ),
  ) as State;

const rk4 = async (
  timeS: number,
  state: State,
  stepS: number,
  gravity: GravityStore,
  config: SolverConfig,
): Promise<State> => {
  const k1 = await derivative(timeS, state, gravity, config);
  const k2 = await derivative(
    timeS + stepS / 2,
    combineState(state, [k1], [0.5], stepS),
    gravity,
    config,
  );
  const k3 = await derivative(
    timeS + stepS / 2,
    combineState(state, [k2], [0.5], stepS),
    gravity,
    config,
  );
  const k4 = await derivative(
    timeS + stepS,
    combineState(state, [k3], [1], stepS),
    gravity,
    config,
  );
  return combineState(state, [k1, k2, k3, k4], [1 / 6, 1 / 3, 1 / 3, 1 / 6], stepS);
};

const stepWithError = async (
  timeS: number,
  state: State,
  stepS: number,
  gravity: GravityStore,
  config: SolverConfig,
): Promise<{ state: State; localPositionErrorM: number }> => {
  const full = await rk4(timeS, state, stepS, gravity, config);
  const half = await rk4(timeS, state, stepS / 2, gravity, config);
  const twoHalf = await rk4(timeS + stepS / 2, half, stepS / 2, gravity, config);
  const error = Math.hypot(
    twoHalf[0] - full[0],
    twoHalf[1] - full[1],
    twoHalf[2] - full[2],
  ) / 15;
  return { state: twoHalf, localPositionErrorM: error };
};

const pointFromState = async (
  timeS: number,
  state: State,
  terrain: TerrainStore,
): Promise<TrajectoryPoint> => {
  const position: Vec3 = [state[0], state[1], state[2]];
  const [latitude, longitude, radius] = geodetic(position);
  const terrainElevation = await terrain.elevationM(latitude, longitude);
  const altitude = radius - MOON_RADIUS_M;
  return {
    timeS,
    xBodyM: state[0],
    yBodyM: state[1],
    zBodyM: state[2],
    vxBodyMS: state[3],
    vyBodyMS: state[4],
    vzBodyMS: state[5],
    latitudeDeg: latitude,
    longitudeDegEast: longitude,
    altitudeM: altitude,
    terrainElevationM: terrainElevation,
    clearanceM: altitude - terrainElevation - BALL_RADIUS_M,
  };
};

const hermite = (
  left: TrajectoryPoint,
  right: TrajectoryPoint,
  fraction: number,
): { position: Vec3; velocity: Vec3; timeS: number } => {
  const dt = right.timeS - left.timeS;
  const u = fraction;
  const u2 = u * u;
  const u3 = u2 * u;
  const h00 = 2 * u3 - 3 * u2 + 1;
  const h10 = u3 - 2 * u2 + u;
  const h01 = -2 * u3 + 3 * u2;
  const h11 = u3 - u2;
  const p0: Vec3 = [left.xBodyM, left.yBodyM, left.zBodyM];
  const p1: Vec3 = [right.xBodyM, right.yBodyM, right.zBodyM];
  const v0: Vec3 = [left.vxBodyMS, left.vyBodyMS, left.vzBodyMS];
  const v1: Vec3 = [right.vxBodyMS, right.vyBodyMS, right.vzBodyMS];
  const position = add3(
    add3(scale3(p0, h00), scale3(v0, h10 * dt)),
    add3(scale3(p1, h01), scale3(v1, h11 * dt)),
  );
  const dh00 = (6 * u2 - 6 * u) / dt;
  const dh10 = 3 * u2 - 4 * u + 1;
  const dh01 = (-6 * u2 + 6 * u) / dt;
  const dh11 = 3 * u2 - 2 * u;
  const velocity = add3(
    add3(scale3(p0, dh00), scale3(v0, dh10)),
    add3(scale3(p1, dh01), scale3(v1, dh11)),
  );
  return { position, velocity, timeS: left.timeS + u * dt };
};

const goldenMinimum = (
  left: TrajectoryPoint,
  right: TrajectoryPoint,
  site: Vec3,
): { distanceM: number; timeS: number; relativeSpeedMS: number } => {
  const ratio = (Math.sqrt(5) - 1) / 2;
  let a = 0;
  let b = 1;
  let c = b - ratio * (b - a);
  let d = a + ratio * (b - a);
  const objective = (fraction: number): number =>
    distance3(hermite(left, right, fraction).position, site);
  let fc = objective(c);
  let fd = objective(d);
  for (let iteration = 0; iteration < 42; iteration += 1) {
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - ratio * (b - a);
      fc = objective(c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + ratio * (b - a);
      fd = objective(d);
    }
  }
  const fraction = (a + b) / 2;
  const value = hermite(left, right, fraction);
  return {
    distanceM: distance3(value.position, site),
    timeS: value.timeS,
    relativeSpeedMS: norm3(value.velocity),
  };
};

const closestReturn = (
  points: TrajectoryPoint[],
  site: Vec3,
  minimumTimeS: number,
): { distanceM: number; timeS: number; relativeSpeedMS: number } | null => {
  let best: { distanceM: number; timeS: number; relativeSpeedMS: number } | null = null;
  for (let index = 1; index < points.length; index += 1) {
    if (points[index].timeS < minimumTimeS) continue;
    const leftDistance = distance3(
      [points[index - 1].xBodyM, points[index - 1].yBodyM, points[index - 1].zBodyM],
      site,
    );
    const rightDistance = distance3(
      [points[index].xBodyM, points[index].yBodyM, points[index].zBodyM],
      site,
    );
    if (
      index + 1 < points.length &&
      rightDistance >
        distance3(
          [points[index + 1].xBodyM, points[index + 1].yBodyM, points[index + 1].zBodyM],
          site,
        ) &&
      leftDistance > rightDistance
    ) {
      continue;
    }
    const candidate = goldenMinimum(points[index - 1], points[index], site);
    if (!best || candidate.distanceM < best.distanceM) best = candidate;
  }
  return best;
};

const refineImpact = async (
  left: TrajectoryPoint,
  right: TrajectoryPoint,
  terrain: TerrainStore,
): Promise<TrajectoryPoint> => {
  let low = 0;
  let high = 1;
  let impact = right;
  for (let iteration = 0; iteration < 24; iteration += 1) {
    const fraction = (low + high) / 2;
    const value = hermite(left, right, fraction);
    const state: State = [
      ...value.position,
      ...value.velocity,
    ] as State;
    const point = await pointFromState(value.timeS, state, terrain);
    if (point.clearanceM > 0) low = fraction;
    else {
      high = fraction;
      impact = point;
    }
  }
  return impact;
};

const solve = async (
  launch: LaunchInput,
  config: SolverConfig,
  baseUrl: string,
): Promise<SolverResult> => {
  const started = performance.now();
  const terrain = new TerrainStore(baseUrl);
  const gravity = new GravityStore(baseUrl);
  const [, gravityMeta] = await Promise.all([terrain.initialize(), gravity.initialize()]);
  const launchTerrain = await terrain.elevationM(
    launch.latitudeDeg,
    launch.longitudeDegEast,
  );
  const launchRadius =
    MOON_RADIUS_M + launchTerrain + launch.heightAboveTerrainM;
  if (launch.heightAboveTerrainM <= BALL_RADIUS_M) {
    throw new Error("Ball-centre height must exceed the baseball radius.");
  }
  const initial = initialState(launch, launchRadius);
  let state = initial.state;
  let timeS = 0;
  let numericalError = 0;
  const points: TrajectoryPoint[] = [await pointFromState(0, state, terrain)];
  let impact: TrajectoryPoint | null = null;
  const steps = Math.ceil(config.durationS / config.stepS);

  for (let index = 0; index < steps; index += 1) {
    const step = Math.min(config.stepS, config.durationS - timeS);
    if (step <= 0) break;
    const next = await stepWithError(timeS, state, step, gravity, config);
    numericalError += next.localPositionErrorM;
    const nextTime = timeS + step;
    const point = await pointFromState(nextTime, next.state, terrain);
    if (point.clearanceM <= 0) {
      impact = await refineImpact(points[points.length - 1], point, terrain);
      points.push(impact);
      state = next.state;
      timeS = impact.timeS;
      break;
    }
    points.push(point);
    state = next.state;
    timeS = nextTime;
    if (index % Math.max(1, Math.floor(steps / 100)) === 0) {
      const progress: ProgressMessage = {
        type: "progress",
        fraction: index / steps,
        label: `Propagating ${Math.round((index / steps) * 100)}%`,
      };
      ctx.postMessage(progress);
    }
    const radius = Math.hypot(state[0], state[1], state[2]);
    const speed2 = state[3] ** 2 + state[4] ** 2 + state[5] ** 2;
    if (radius > 10 * MOON_RADIUS_M && 0.5 * speed2 - MU_MOON / radius > 0) break;
  }

  const closest = closestReturn(points, initial.site, config.minimumReturnTimeS);
  const minimumClearance = Math.min(...points.map((point) => point.clearanceM));
  const maximumAltitude = Math.max(...points.map((point) => point.altitudeM));
  const minimumRadius = Math.min(
    ...points.map((point) =>
      Math.hypot(point.xBodyM, point.yBodyM, point.zBodyM),
    ),
  );
  let status: SolverResult["status"];
  let message: string;
  const effectivePositionErrorBound = Math.max(
    config.assumedTotalPositionErrorBoundM,
    numericalError,
  );
  if (
    closest &&
    closest.distanceM + effectivePositionErrorBound <=
      config.returnToleranceM
  ) {
    status = "RETURN";
    message =
      "Conditional return: the miss plus the larger of the supplied assumed total position-error bound and the RK4 numerical diagnostic lies inside the requested return sphere. The solver does not derive a physical uncertainty bound.";
  } else if (
    closest &&
    closest.distanceM - effectivePositionErrorBound <=
      config.returnToleranceM
  ) {
    status = "RETURN_UNCERTAIN";
    message =
      "The conditional interval crosses the return boundary. The supplied position-error bound is an assumption, and the RK4 quantity is a numerical diagnostic rather than a certified bound.";
  } else if (impact) {
    status = "IMPACT";
    message = "The ball intersects the LDEM64 surface before a qualified return.";
  } else {
    const final = points[points.length - 1];
    const radius = Math.hypot(final.xBodyM, final.yBodyM, final.zBodyM);
    const speed2 = final.vxBodyMS ** 2 + final.vyBodyMS ** 2 + final.vzBodyMS ** 2;
    if (radius > 2 * MOON_RADIUS_M && 0.5 * speed2 - MU_MOON / radius > 0) {
      status = "ESCAPE";
      message = "The trajectory is unbound and leaves the lunar vicinity.";
    } else {
      status = "NO_RETURN_WITHIN_WINDOW";
      message = "No qualified return occurs inside the selected time window.";
    }
  }
  const domainOk = minimumRadius >= gravityMeta.brillouin_radius_m;
  if (!domainOk) {
    message +=
      " Part of the arc is below the conservative Brillouin sphere, where the external harmonic field is an explicitly flagged extrapolation.";
  }

  return {
    status,
    message,
    points,
    launchTerrainElevationM: launchTerrain,
    launchRadiusM: launchRadius,
    gravityDegree: gravityMeta.degree,
    closestReturnTimeS: closest?.timeS ?? null,
    closestReturnDistanceM: closest?.distanceM ?? null,
    closestReturnRelativeSpeedMS: closest?.relativeSpeedMS ?? null,
    impactTimeS: impact?.timeS ?? null,
    impactLatitudeDeg: impact?.latitudeDeg ?? null,
    impactLongitudeDegEast: impact?.longitudeDegEast ?? null,
    impactSpeedMS: impact
      ? Math.hypot(impact.vxBodyMS, impact.vyBodyMS, impact.vzBodyMS)
      : null,
    minimumClearanceM: minimumClearance,
    maximumAltitudeM: maximumAltitude,
    numericalPositionUncertaintyM: numericalError,
    entireArcOutsideBrillouinSphere: domainOk,
    elapsedMS: performance.now() - started,
  };
};

ctx.onmessage = async (event: MessageEvent<RunMessage>) => {
  if (event.data.type !== "run") return;
  try {
    const result = await solve(event.data.launch, event.data.config, event.data.baseUrl);
    const message: ResultMessage = { type: "result", result };
    ctx.postMessage(message);
  } catch (error) {
    const message: ErrorMessage = {
      type: "error",
      message: error instanceof Error ? error.message : String(error),
    };
    ctx.postMessage(message);
  }
};

export {};
