export type SolverStatus =
  | "RETURN"
  | "RETURN_UNCERTAIN"
  | "IMPACT"
  | "ESCAPE"
  | "NO_RETURN_WITHIN_WINDOW"
  | "NUMERICAL_FAILURE";

export interface LaunchInput {
  latitudeDeg: number;
  longitudeDegEast: number;
  heightAboveTerrainM: number;
  speedMS: number;
  elevationDeg: number;
  azimuthDeg: number;
}

export interface SolverConfig {
  durationS: number;
  stepS: number;
  returnToleranceM: number;
  assumedTotalPositionErrorBoundM: number;
  minimumReturnTimeS: number;
  includeEarthTide: boolean;
  includeSolarTide: boolean;
}

export interface TrajectoryPoint {
  timeS: number;
  xBodyM: number;
  yBodyM: number;
  zBodyM: number;
  vxBodyMS: number;
  vyBodyMS: number;
  vzBodyMS: number;
  latitudeDeg: number;
  longitudeDegEast: number;
  altitudeM: number;
  terrainElevationM: number;
  clearanceM: number;
}

export interface SolverResult {
  status: SolverStatus;
  message: string;
  points: TrajectoryPoint[];
  launchTerrainElevationM: number;
  launchRadiusM: number;
  gravityDegree: number;
  closestReturnTimeS: number | null;
  closestReturnDistanceM: number | null;
  closestReturnRelativeSpeedMS: number | null;
  impactTimeS: number | null;
  impactLatitudeDeg: number | null;
  impactLongitudeDegEast: number | null;
  impactSpeedMS: number | null;
  minimumClearanceM: number;
  maximumAltitudeM: number;
  numericalPositionUncertaintyM: number;
  entireArcOutsideBrillouinSphere: boolean;
  elapsedMS: number;
}

export interface TerrainManifest {
  dataset: string;
  reference_radius_m: number;
  pixels_per_degree: number;
  tile_degrees: number;
  tile_rows: number;
  tile_columns: number;
  scale_m_per_dn: number;
  predictor: string;
}

export interface GravityManifest {
  model: string;
  degree: number;
  reference_radius_m: number;
  gm_m3_s2: number;
  stores_noncentral_correction: boolean;
  brillouin_radius_m: number;
  pixels_per_degree: number;
  tile_degrees: number;
  tile_rows: number;
  tile_columns: number;
  altitude_shells_m: number[];
  component_order: ["radial", "theta", "phi"];
}

export interface RunMessage {
  type: "run";
  launch: LaunchInput;
  config: SolverConfig;
  baseUrl: string;
}

export interface ProgressMessage {
  type: "progress";
  fraction: number;
  label: string;
}

export interface ResultMessage {
  type: "result";
  result: SolverResult;
}

export interface ErrorMessage {
  type: "error";
  message: string;
}
