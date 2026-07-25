import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import type { SolverResult } from "./types";

const toScene = (x: number, y: number, z: number): THREE.Vector3 =>
  new THREE.Vector3(x, z, -y).divideScalar(1_737_400);

export class MoonView {
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(34, 1, 0.01, 100);
  private readonly renderer = new THREE.WebGLRenderer({
    antialias: true,
    alpha: true,
    powerPreference: "high-performance",
  });
  private readonly controls: OrbitControls;
  private orbit?: THREE.Line;
  private launchMarker?: THREE.Mesh;
  private eventMarker?: THREE.Mesh;
  private readonly resizeObserver: ResizeObserver;

  constructor(
    private readonly container: HTMLElement,
    baseUrl: string,
  ) {
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 0.9;
    container.append(this.renderer.domElement);

    this.camera.position.set(2.7, 1.45, 2.8);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.05;
    this.controls.minDistance = 1.35;
    this.controls.maxDistance = 8;

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
    this.scene.add(moon);

    this.scene.add(new THREE.HemisphereLight(0x8298a8, 0x05080a, 0.34));
    const sun = new THREE.DirectionalLight(0xfff4df, 4.2);
    sun.position.set(3.5, 2.2, 4.5);
    this.scene.add(sun);
    const rim = new THREE.DirectionalLight(0x5ba9c9, 0.8);
    rim.position.set(-4, 0.5, -3);
    this.scene.add(rim);

    const atmosphere = new THREE.Mesh(
      new THREE.SphereGeometry(1.012, 96, 48),
      new THREE.MeshBasicMaterial({
        color: 0x88c7de,
        transparent: true,
        opacity: 0.025,
        side: THREE.BackSide,
      }),
    );
    this.scene.add(atmosphere);

    this.resizeObserver = new ResizeObserver(() => this.resize());
    this.resizeObserver.observe(container);
    this.resize();
    this.animate();
  }

  private resize(): void {
    const width = Math.max(1, this.container.clientWidth);
    const height = Math.max(1, this.container.clientHeight);
    this.camera.aspect = width / height;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(width, height, false);
  }

  private animate = (): void => {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    requestAnimationFrame(this.animate);
  };

  update(result: SolverResult): void {
    if (this.orbit) {
      this.scene.remove(this.orbit);
      this.orbit.geometry.dispose();
      (this.orbit.material as THREE.Material).dispose();
    }
    for (const marker of [this.launchMarker, this.eventMarker]) {
      if (marker) {
        this.scene.remove(marker);
        marker.geometry.dispose();
        (marker.material as THREE.Material).dispose();
      }
    }
    const points = result.points.map((point) =>
      toScene(point.xBodyM, point.yBodyM, point.zBodyM),
    );
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({
      color: result.status === "IMPACT" ? 0xffa15c : 0x7de3f4,
      transparent: true,
      opacity: 0.96,
    });
    this.orbit = new THREE.Line(geometry, material);
    this.scene.add(this.orbit);

    const markerGeometry = new THREE.SphereGeometry(0.018, 24, 12);
    this.launchMarker = new THREE.Mesh(
      markerGeometry,
      new THREE.MeshBasicMaterial({ color: 0x7ff0b5 }),
    );
    this.launchMarker.position.copy(points[0]);
    this.scene.add(this.launchMarker);

    if (result.impactTimeS !== null || result.closestReturnTimeS !== null) {
      const eventTime =
        result.impactTimeS ?? result.closestReturnTimeS ?? result.points.at(-1)!.timeS;
      const eventPoint = result.points.reduce((closest, point) =>
        Math.abs(point.timeS - eventTime) < Math.abs(closest.timeS - eventTime)
          ? point
          : closest,
      );
      this.eventMarker = new THREE.Mesh(
        new THREE.SphereGeometry(0.022, 24, 12),
        new THREE.MeshBasicMaterial({
          color: result.impactTimeS !== null ? 0xff765f : 0xffd36a,
        }),
      );
      this.eventMarker.position.copy(
        toScene(eventPoint.xBodyM, eventPoint.yBodyM, eventPoint.zBodyM),
      );
      this.scene.add(this.eventMarker);
    }
  }
}
