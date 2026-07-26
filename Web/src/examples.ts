import type { LaunchInput, SolverConfig, SolverStatus } from "./types";

export interface TrajectoryExample {
  status: Exclude<SolverStatus, "NUMERICAL_FAILURE">;
  launch: LaunchInput;
  config: SolverConfig;
}

const config = (
  durationS: number,
  returnToleranceM = 10,
  assumedTotalPositionErrorBoundM = 10,
): SolverConfig => ({
  durationS,
  stepS: 4,
  returnToleranceM,
  assumedTotalPositionErrorBoundM,
  minimumReturnTimeS: 1000,
  includeEarthTide: true,
  includeSolarTide: true,
});

export const trajectoryExamples: TrajectoryExample[] = [
  {
    status: "RETURN",
    launch: {
      latitudeDeg: 5.4296875,
      longitudeDegEast: 201.3671875,
      heightAboveTerrainM: 19_243,
      speedMS: 1661.4559,
      elevationDeg: 0.531603,
      azimuthDeg: 88.383067,
    },
    config: config(6686.34, 12, 10),
  },
  {
    status: "RETURN_UNCERTAIN",
    launch: {
      latitudeDeg: 5.4296875,
      longitudeDegEast: 201.3671875,
      heightAboveTerrainM: 19_243,
      speedMS: 1661.4559,
      elevationDeg: 0.531603,
      azimuthDeg: 88.383067,
    },
    config: config(6686.34, 10, 10),
  },
  {
    status: "IMPACT",
    launch: {
      latitudeDeg: 0,
      longitudeDegEast: 0,
      heightAboveTerrainM: 1,
      speedMS: 100,
      elevationDeg: 45,
      azimuthDeg: 90,
    },
    config: config(1200),
  },
  {
    status: "ESCAPE",
    launch: {
      latitudeDeg: 0,
      longitudeDegEast: 0,
      heightAboveTerrainM: 1000,
      speedMS: 2600,
      elevationDeg: 90,
      azimuthDeg: 90,
    },
    config: config(12_000),
  },
  {
    status: "NO_RETURN_WITHIN_WINDOW",
    launch: {
      latitudeDeg: 0,
      longitudeDegEast: 0,
      heightAboveTerrainM: 100_000,
      speedMS: 1635,
      elevationDeg: 0,
      azimuthDeg: 90,
    },
    config: config(2000),
  },
];
