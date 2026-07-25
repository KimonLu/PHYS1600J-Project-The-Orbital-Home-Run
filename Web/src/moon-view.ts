import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import type { SolverResult, TrajectoryPoint } from "./types";

const MOON_RADIUS_M = 1_737_400;
const MOON_ROTATION_PERIOD_S = 27.321661 * 86_400;
const MOON_ANGULAR_SPEED_RAD_S = (2 * Math.PI) / MOON_ROTATION_PERIOD_S;
const SPIN_AXIS = new THREE.Vector3(0, 1, 0);

const bodyToScene = (x: number, y: number, z: number): THREE.Vector3 =>
  new THREE.Vector3(x, z, -y).divideScalar(MOON_RADIUS_M);

const inertialPosition = (point: TrajectoryPoint): THREE.Vector3 =>
  bodyToScene(point.xBodyM, point.yBodyM, point.zBodyM).applyAxisAngle(
    SPIN_AXIS,
    MOON_ANGULAR_SPEED_RAD_S * point.timeS,
  );

export interface PlaybackState {
  playing: boolean;
  timeS: number;
  durationS: number;
  speed: number;
}

export class MoonView {
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(34, 1, 0.01, 100);
  private readonly renderer = new THREE.WebGLRenderer({
    antialias: true,
    alpha: true,
    powerPreference: "high-performance",
  });
  private readonly controls: OrbitControls;
  private readonly moonGroup = new THREE.Group();
  private readonly pressedKeys = new Set<string>();
  private readonly resizeObserver: ResizeObserver;
  private orbit?: THREE.Line;
  private launchMarker?: THREE.Mesh;
  private eventMarker?: THREE.Mesh;
  private ballMarker?: THREE.Mesh;
  private trajectory: TrajectoryPoint[] = [];
  private inertialPoints: THREE.Vector3[] = [];
  private simulationTimeS = 0;
  private durationS = 0;
  private playbackSpeed = 100;
  private playing = false;
  private lastFrameMS = performance.now();

  constructor(
    private readonly container: HTMLElement,
    baseUrl: string,
    private readonly onPlayback?: (state: PlaybackState) => void,
  ) {
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 0.95;
    this.renderer.domElement.tabIndex = 0;
    container.append(this.renderer.domElement);

    this.camera.position.set(2.7, 1.45, 2.8);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.set(0, 0, 0);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.05;
    this.controls.minDistance = 1.35;
    this.controls.maxDistance = 12;
    this.controls.update();

    const loader = new THREE.TextureLoader();
    const color = loader.load(`${baseUrl}assets/lroc_color_2k.jpg`);
    color.colorSpace = THREE.SRGBColorSpace;
    color.wrapS = THREE.RepeatWrapping;
    const displacement = loader.load(`${baseUrl}assets/ldem_3_8bit.jpg`);
    displacement.wrapS = THREE.RepeatWrapping;
    const geometry = new THREE.SphereGeometry(1, 256, 128);
    const material = new THREE.MeshStandardMaterial({
      map: color,
      displacementMap: displacement,
      displacementScale: 0.012,
      displacementBias: -0.006,
      roughness: 0.96,
      metalness: 0,
    });
    const moon = new THREE.Mesh(geometry, material);
    moon.rotation.y = -Math.PI / 2;
    this.moonGroup.add(moon);
    this.scene.add(this.moonGroup);

    // Uniform illumination deliberately avoids an artificial day/night terminator.
    this.scene.add(new THREE.AmbientLight(0xffffff, 2.8));
    this.addSpinAxis();
    this.bindKeyboardPan();

    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(container);
    this.resize();
    this.animate();
  }

  setPlaybackSpeed(speed: number): void {
    if (!Number.isFinite(speed) || speed <= 0) return;
    this.playbackSpeed = speed;
    this.emitPlayback();
  }

  togglePlayback(): void {
    if (this.trajectory.length === 0) return;
    if (this.simulationTimeS >= this.durationS) this.simulationTimeS = 0;
    this.playing = !this.playing;
    this.lastFrameMS = performance.now();
    this.updateAnimationObjects();
    this.emitPlayback();
  }

  restartPlayback(): void {
    if (this.trajectory.length === 0) return;
    this.simulationTimeS = 0;
    this.playing = true;
    this.lastFrameMS = performance.now();
    this.updateAnimationObjects();
    this.emitPlayback();
  }

  private addSpinAxis(): void {
    const geometry = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, -1.42, 0),
      new THREE.Vector3(0, 1.42, 0),
    ]);
    const material = new THREE.LineDashedMaterial({
      color: 0xffd36a,
      dashSize: 0.055,
      gapSize: 0.03,
      transparent: true,
      opacity: 0.88,
      depthTest: false,
    });
    const line = new THREE.Line(geometry, material);
    line.computeLineDistances();
    line.renderOrder = 10;
    this.scene.add(line);

    const coneGeometry = new THREE.ConeGeometry(0.035, 0.11, 16);
    const coneMaterial = new THREE.MeshBasicMaterial({
      color: 0xffd36a,
      depthTest: false,
      transparent: true,
      opacity: 0.9,
    });
    const north = new THREE.Mesh(coneGeometry, coneMaterial);
    north.position.y = 1.46;
    north.renderOrder = 10;
    this.scene.add(north);
    const south = new THREE.Mesh(coneGeometry, coneMaterial.clone());
    south.rotation.z = Math.PI;
    south.position.y = -1.46;
    south.renderOrder = 10;
    this.scene.add(south);
  }

  private bindKeyboardPan(): void {
    const canvas = this.renderer.domElement;
    canvas.addEventListener("pointerdown", () => canvas.focus());
    canvas.addEventListener("keydown", (event) => {
      const key = event.key.toLowerCase();
      if (!["w", "a", "s", "d"].includes(key)) return;
      event.preventDefault();
      this.pressedKeys.add(key);
    });
    canvas.addEventListener("keyup", (event) => {
      this.pressedKeys.delete(event.key.toLowerCase());
    });
    canvas.addEventListener("blur", () => this.pressedKeys.clear());
  }

  private panFromKeyboard(deltaS: number): void {
    if (this.pressedKeys.size === 0) return;
    const distance = this.camera.position.distanceTo(this.controls.target);
    const amount = distance * 0.55 * deltaS;
    const right = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 0);
    const up = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 1);
    const movement = new THREE.Vector3();
    if (this.pressedKeys.has("a")) movement.addScaledVector(right, -amount);
    if (this.pressedKeys.has("d")) movement.addScaledVector(right, amount);
    if (this.pressedKeys.has("w")) movement.addScaledVector(up, amount);
    if (this.pressedKeys.has("s")) movement.addScaledVector(up, -amount);
    this.camera.position.add(movement);
    this.controls.target.add(movement);
  }

  private resize(): void {
    const width = Math.max(1, this.container.clientWidth);
    const height = Math.max(1, this.container.clientHeight);
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(width, height, false);
  }

  private animate = (frameMS = performance.now()): void => {
    const realDeltaS = Math.min(0.1, Math.max(0, (frameMS - this.lastFrameMS) / 1000));
    this.lastFrameMS = frameMS;
    this.panFromKeyboard(realDeltaS);
    if (this.playing && this.durationS > 0) {
      this.simulationTimeS += realDeltaS * this.playbackSpeed;
      if (this.simulationTimeS >= this.durationS) {
        this.simulationTimeS = this.durationS;
        this.playing = false;
      }
      this.updateAnimationObjects();
      this.emitPlayback();
    }
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    requestAnimationFrame(this.animate);
  };

  private updateAnimationObjects(): void {
    this.moonGroup.rotation.y = MOON_ANGULAR_SPEED_RAD_S * this.simulationTimeS;
    if (!this.ballMarker || this.trajectory.length === 0) return;
    let high = this.trajectory.findIndex((point) => point.timeS >= this.simulationTimeS);
    if (high < 0) high = this.trajectory.length - 1;
    const low = Math.max(0, high - 1);
    const left = this.trajectory[low];
    const right = this.trajectory[high];
    const span = right.timeS - left.timeS;
    const fraction = span > 0 ? (this.simulationTimeS - left.timeS) / span : 0;
    this.ballMarker.position
      .copy(this.inertialPoints[low])
      .lerp(this.inertialPoints[high], Math.min(1, Math.max(0, fraction)));
  }

  private emitPlayback(): void {
    this.onPlayback?.({
      playing: this.playing,
      timeS: this.simulationTimeS,
      durationS: this.durationS,
      speed: this.playbackSpeed,
    });
  }

  update(result: SolverResult): void {
    if (this.orbit) {
      this.scene.remove(this.orbit);
      this.orbit.geometry.dispose();
      (this.orbit.material as THREE.Material).dispose();
    }
    for (const marker of [this.launchMarker, this.eventMarker, this.ballMarker]) {
      if (marker) {
        marker.removeFromParent();
        marker.geometry.dispose();
        (marker.material as THREE.Material).dispose();
      }
    }

    this.trajectory = result.points;
    this.inertialPoints = result.points.map(inertialPosition);
    this.durationS = result.points.at(-1)?.timeS ?? 0;
    this.simulationTimeS = 0;
    this.playing = true;
    this.lastFrameMS = performance.now();
    this.moonGroup.rotation.y = 0;

    const geometry = new THREE.BufferGeometry().setFromPoints(this.inertialPoints);
    const material = new THREE.LineBasicMaterial({
      color: result.status === "IMPACT" ? 0xffa15c : 0x7de3f4,
      transparent: true,
      opacity: 0.96,
    });
    this.orbit = new THREE.Line(geometry, material);
    this.scene.add(this.orbit);

    this.launchMarker = new THREE.Mesh(
      new THREE.SphereGeometry(0.018, 24, 12),
      new THREE.MeshBasicMaterial({ color: 0x7ff0b5 }),
    );
    this.launchMarker.position.copy(
      bodyToScene(
        result.points[0].xBodyM,
        result.points[0].yBodyM,
        result.points[0].zBodyM,
      ),
    );
    this.moonGroup.add(this.launchMarker);

    this.ballMarker = new THREE.Mesh(
      new THREE.SphereGeometry(0.024, 24, 12),
      new THREE.MeshBasicMaterial({ color: 0xffffff }),
    );
    this.ballMarker.position.copy(this.inertialPoints[0]);
    this.scene.add(this.ballMarker);

    if (result.impactTimeS !== null || result.closestReturnTimeS !== null) {
      const eventTime =
        result.impactTimeS ?? result.closestReturnTimeS ?? result.points.at(-1)!.timeS;
      const eventIndex = result.points.reduce(
        (closest, point, index) =>
          Math.abs(point.timeS - eventTime) <
          Math.abs(result.points[closest].timeS - eventTime)
            ? index
            : closest,
        0,
      );
      this.eventMarker = new THREE.Mesh(
        new THREE.SphereGeometry(0.022, 24, 12),
        new THREE.MeshBasicMaterial({
          color: result.impactTimeS !== null ? 0xff765f : 0xffd36a,
        }),
      );
      this.eventMarker.position.copy(this.inertialPoints[eventIndex]);
      this.scene.add(this.eventMarker);
    }
    this.updateAnimationObjects();
    this.emitPlayback();
  }
}
