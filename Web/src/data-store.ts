import type { GravityManifest, TerrainManifest } from "./types";

const decompressGzip = async (response: Response): Promise<ArrayBuffer> => {
  if (!response.ok) {
    throw new Error(`Data request failed: ${response.status} ${response.url}`);
  }
  if (!response.body) throw new Error("Streaming response is unavailable");
  const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
  return new Response(stream).arrayBuffer();
};

const wrapLongitude = (longitudeDeg: number): number =>
  ((longitudeDeg % 360) + 360) % 360;

export class TerrainStore {
  private manifest?: TerrainManifest;
  private cache = new Map<string, Int16Array>();

  constructor(private readonly baseUrl: string) {}

  async initialize(): Promise<TerrainManifest> {
    const response = await fetch(`${this.baseUrl}data/terrain/manifest.json`);
    if (!response.ok) throw new Error("LDEM64 terrain manifest is unavailable");
    this.manifest = (await response.json()) as TerrainManifest;
    return this.manifest;
  }

  private async tile(latIndex: number, lonIndex: number): Promise<Int16Array> {
    const key = `${latIndex}:${lonIndex}`;
    const cached = this.cache.get(key);
    if (cached) return cached;
    const filename = `t_${latIndex.toString().padStart(2, "0")}_${lonIndex
      .toString()
      .padStart(2, "0")}.i16.gz`;
    const buffer = await decompressGzip(
      await fetch(`${this.baseUrl}data/terrain/${filename}`),
    );
    const encoded = new Int16Array(buffer);
    const tile = new Int16Array(encoded.length);
    const columns = this.manifest!.tile_columns;
    for (let row = 0; row < this.manifest!.tile_rows; row += 1) {
      let value = 0;
      const offset = row * columns;
      for (let column = 0; column < columns; column += 1) {
        value = (value + encoded[offset + column]) << 16 >> 16;
        tile[offset + column] = value;
      }
    }
    this.cache.set(key, tile);
    return tile;
  }

  async elevationM(latitudeDeg: number, longitudeDeg: number): Promise<number> {
    if (!this.manifest) await this.initialize();
    const meta = this.manifest!;
    const lat = Math.max(
      -90 + 0.5 / meta.pixels_per_degree,
      Math.min(90 - 0.5 / meta.pixels_per_degree, latitudeDeg),
    );
    const lon = wrapLongitude(longitudeDeg);
    const latIndex = Math.min(
      180 / meta.tile_degrees - 1,
      Math.floor((lat + 90) / meta.tile_degrees),
    );
    const lonIndex = Math.floor(lon / meta.tile_degrees);
    const south = -90 + latIndex * meta.tile_degrees;
    const north = south + meta.tile_degrees;
    const west = lonIndex * meta.tile_degrees;
    const row = (north - lat) * meta.pixels_per_degree - 0.5;
    const col = (lon - west) * meta.pixels_per_degree - 0.5;
    // Keep the un-clamped neighbouring indices. `sample` maps them into the
    // adjacent latitude/longitude tile, so interpolation remains continuous
    // across every 10-degree seam (including the 0/360-degree meridian).
    const r0 = Math.floor(row);
    const c0 = Math.floor(col);
    const fr = row - Math.floor(row);
    const fc = col - Math.floor(col);

    const tile00 = await this.tile(latIndex, lonIndex);
    const sample = async (rr: number, cc: number): Promise<number> => {
      let li = latIndex;
      let lo = lonIndex;
      let r = rr;
      let c = cc;
      if (r < 0) {
        li = Math.min(180 / meta.tile_degrees - 1, li + 1);
        r += meta.tile_rows;
      } else if (r >= meta.tile_rows) {
        li = Math.max(0, li - 1);
        r -= meta.tile_rows;
      }
      if (c < 0) {
        lo = (lo - 1 + 360 / meta.tile_degrees) % (360 / meta.tile_degrees);
        c += meta.tile_columns;
      } else if (c >= meta.tile_columns) {
        lo = (lo + 1) % (360 / meta.tile_degrees);
        c -= meta.tile_columns;
      }
      const data = li === latIndex && lo === lonIndex ? tile00 : await this.tile(li, lo);
      return data[r * meta.tile_columns + c];
    };
    const [z00, z01, z10, z11] = await Promise.all([
      sample(r0, c0),
      sample(r0, c0 + 1),
      sample(r0 + 1, c0),
      sample(r0 + 1, c0 + 1),
    ]);
    return (
      ((1 - fr) * (1 - fc) * z00 +
        (1 - fr) * fc * z01 +
        fr * (1 - fc) * z10 +
        fr * fc * z11) *
      meta.scale_m_per_dn
    );
  }
}

export class GravityStore {
  private manifest?: GravityManifest;
  private cache = new Map<string, Float32Array>();

  constructor(private readonly baseUrl: string) {}

  async initialize(): Promise<GravityManifest> {
    const response = await fetch(`${this.baseUrl}data/gravity/manifest.json`);
    if (!response.ok) throw new Error("GRAIL gravity manifest is unavailable");
    this.manifest = (await response.json()) as GravityManifest;
    return this.manifest;
  }

  private async tile(
    shellIndex: number,
    latIndex: number,
    lonIndex: number,
  ): Promise<Float32Array> {
    const key = `${shellIndex}:${latIndex}:${lonIndex}`;
    const cached = this.cache.get(key);
    if (cached) return cached;
    const filename = `g_${shellIndex.toString().padStart(2, "0")}_${latIndex
      .toString()
      .padStart(2, "0")}_${lonIndex.toString().padStart(2, "0")}.f32.gz`;
    const buffer = await decompressGzip(
      await fetch(`${this.baseUrl}data/gravity/${filename}`),
    );
    const tile = new Float32Array(buffer);
    this.cache.set(key, tile);
    return tile;
  }

  private async sampleShell(
    shellIndex: number,
    latitudeDeg: number,
    longitudeDeg: number,
  ): Promise<[number, number, number]> {
    const meta = this.manifest!;
    const lat = Math.max(
      -90 + 0.5 / meta.pixels_per_degree,
      Math.min(90 - 0.5 / meta.pixels_per_degree, latitudeDeg),
    );
    const lon = wrapLongitude(longitudeDeg);
    const latIndex = Math.min(
      180 / meta.tile_degrees - 1,
      Math.floor((lat + 90) / meta.tile_degrees),
    );
    const lonIndex = Math.floor(lon / meta.tile_degrees);
    const south = -90 + latIndex * meta.tile_degrees;
    const north = south + meta.tile_degrees;
    const west = lonIndex * meta.tile_degrees;
    const row = (north - lat) * meta.pixels_per_degree - 0.5;
    const col = (lon - west) * meta.pixels_per_degree - 0.5;
    const r0 = Math.floor(row);
    const c0 = Math.floor(col);
    const fr = row - r0;
    const fc = col - c0;

    const sample = async (rr: number, cc: number): Promise<[number, number, number]> => {
      let li = latIndex;
      let lo = lonIndex;
      let r = rr;
      let c = cc;
      if (r < 0) {
        li = Math.min(180 / meta.tile_degrees - 1, li + 1);
        r += meta.tile_rows;
      } else if (r >= meta.tile_rows) {
        li = Math.max(0, li - 1);
        r -= meta.tile_rows;
      }
      if (c < 0) {
        lo = (lo - 1 + 360 / meta.tile_degrees) % (360 / meta.tile_degrees);
        c += meta.tile_columns;
      } else if (c >= meta.tile_columns) {
        lo = (lo + 1) % (360 / meta.tile_degrees);
        c -= meta.tile_columns;
      }
      const data = await this.tile(shellIndex, li, lo);
      const offset = (r * meta.tile_columns + c) * 3;
      return [data[offset], data[offset + 1], data[offset + 2]];
    };
    const corners = await Promise.all([
      sample(r0, c0),
      sample(r0, c0 + 1),
      sample(r0 + 1, c0),
      sample(r0 + 1, c0 + 1),
    ]);
    return [0, 1, 2].map((component) => {
      const [q00, q01, q10, q11] = corners.map((value) => value[component]);
      return (
        (1 - fr) * (1 - fc) * q00 +
        (1 - fr) * fc * q01 +
        fr * (1 - fc) * q10 +
        fr * fc * q11
      );
    }) as [number, number, number];
  }

  async sphericalAcceleration(
    radiusM: number,
    latitudeDeg: number,
    longitudeDeg: number,
  ): Promise<[number, number, number]> {
    if (!this.manifest) await this.initialize();
    const meta = this.manifest!;
    const altitude = radiusM - meta.reference_radius_m;
    const shells = meta.altitude_shells_m;
    let low = 0;
    while (low + 1 < shells.length && altitude > shells[low + 1]) low += 1;
    const high = Math.min(shells.length - 1, low + 1);
    const denominator = Math.max(1, shells[high] - shells[low]);
    const fraction = Math.max(0, Math.min(1, (altitude - shells[low]) / denominator));
    const a = await this.sampleShell(low, latitudeDeg, longitudeDeg);
    let interpolated: [number, number, number];
    if (high === low) {
      interpolated = [...a];
    } else {
      const b = await this.sampleShell(high, latitudeDeg, longitudeDeg);
      interpolated = [
        a[0] + fraction * (b[0] - a[0]),
        a[1] + fraction * (b[1] - a[1]),
        a[2] + fraction * (b[2] - a[2]),
      ];
    }
    if (meta.stores_noncentral_correction) {
      interpolated[0] -= meta.gm_m3_s2 / (radiusM * radiusM);
    }
    return interpolated;
  }

  metadata(): GravityManifest {
    if (!this.manifest) throw new Error("Gravity store is not initialized");
    return this.manifest;
  }
}
