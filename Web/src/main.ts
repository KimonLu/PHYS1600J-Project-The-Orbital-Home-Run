import "./style.css";
import { MoonView } from "./moon-view";
import type {
  ErrorMessage,
  LaunchInput,
  ProgressMessage,
  ResultMessage,
  RunMessage,
  SolverConfig,
  SolverResult,
} from "./types";

const baseUrl = new URL(import.meta.env.BASE_URL, window.location.href).href;

document.querySelector<HTMLDivElement>("#app")!.innerHTML = `
  <main class="shell">
    <header class="masthead">
      <div>
        <p class="eyebrow">PHYS1600J · LUNAR TRAJECTORY LABORATORY</p>
        <h1>Orbital Home Run</h1>
        <p class="lede">A terrain-aware test of whether a launched baseball can orbit the Moon, clear real topography, and return to its rotating launch site.</p>
      </div>
      <div class="model-badges" aria-label="Scientific model">
        <span>LOLA <strong>LDEM64</strong></span>
        <span>GRAIL <strong id="degree-badge">loading</strong></span>
        <span>Integrator <strong>RK4 step-doubling</strong></span>
      </div>
    </header>

    <section class="workspace">
      <div class="visual-column">
        <article class="panel moon-panel">
          <div class="panel-heading overlay">
            <div>
              <span class="section-index">01</span>
              <h2>Body-fixed 3D trajectory</h2>
            </div>
            <p>Drag to orbit · scroll to zoom</p>
          </div>
          <div id="moon-view" class="moon-view" aria-label="Interactive 3D Moon and trajectory"></div>
          <div class="legend overlay-legend">
            <span><i class="dot launch"></i>launch</span>
            <span><i class="line orbit"></i>trajectory</span>
            <span><i class="dot event"></i>return / impact</span>
          </div>
        </article>

        <div class="analysis-grid">
          <article class="panel chart-panel">
            <div class="panel-heading">
              <div><span class="section-index">02</span><h2>Equirectangular ground track</h2></div>
              <span class="unit">longitude / latitude</span>
            </div>
            <canvas id="ground-track" aria-label="Ground track on lunar map"></canvas>
          </article>
          <article class="panel chart-panel">
            <div class="panel-heading">
              <div><span class="section-index">03</span><h2>Terrain clearance</h2></div>
              <span class="unit">km above LDEM64</span>
            </div>
            <canvas id="clearance-chart" aria-label="Terrain clearance versus time"></canvas>
          </article>
        </div>
      </div>

      <aside class="control-column">
        <form id="solver-form" class="panel controls">
          <div class="panel-heading">
            <div><span class="section-index">INPUT</span><h2>Launch state</h2></div>
            <button type="button" id="reset-button" class="text-button">Reset</button>
          </div>

          <fieldset>
            <legend>Launch site · body-fixed</legend>
            <div class="field-grid">
              <label>Latitude <span>deg</span><input name="latitude" type="number" min="-90" max="90" step="0.0000001" value="5.4296875" required></label>
              <label>Longitude E <span>deg</span><input name="longitude" type="number" step="0.0000001" value="201.3671875" required></label>
              <label class="wide">Ball-centre height above terrain <span>m</span><input name="height" type="number" min="0.04" step="1" value="19243" required></label>
            </div>
          </fieldset>

          <fieldset>
            <legend>Surface-relative velocity</legend>
            <div class="field-grid">
              <label class="wide">Speed <span>m s⁻¹</span><input name="speed" type="number" min="0" step="0.0001" value="1661.4559" required></label>
              <label>Elevation <span>deg</span><input name="elevation" type="number" min="-90" max="90" step="0.000001" value="0.531603" required></label>
              <label>Azimuth <span>deg</span><input name="azimuth" type="number" step="0.000001" value="88.383067" required></label>
            </div>
            <p class="field-note">Azimuth is clockwise from local north; 90° is east.</p>
          </fieldset>

          <fieldset>
            <legend>Decision and integration</legend>
            <div class="field-grid">
              <label>Time window <span>s</span><input name="duration" type="number" min="10" step="0.01" value="6686.34" required></label>
              <label>Integrator step <span>s</span><input name="step" type="number" min="0.25" max="20" step="0.25" value="4" required></label>
              <label>Return radius <span>m</span><input name="tolerance" type="number" min="0" step="0.1" value="10" required></label>
              <label>Model band <span>m</span><input name="uncertainty" type="number" min="0" step="1" value="10" required></label>
            </div>
            <div class="switch-row">
              <label class="switch"><input name="earth" type="checkbox" checked><span></span>Earth tide</label>
              <label class="switch"><input name="sun" type="checkbox" checked><span></span>Solar tide</label>
            </div>
          </fieldset>

          <button class="run-button" type="submit">
            <span>Run trajectory</span>
            <small>stream tiles · propagate · classify</small>
          </button>
          <div id="progress" class="progress" hidden>
            <div><span id="progress-label">Preparing model</span><strong id="progress-value">0%</strong></div>
            <div class="progress-track"><i id="progress-bar"></i></div>
          </div>
        </form>

        <section class="panel result-panel" aria-live="polite">
          <div class="result-heading">
            <span id="status-dot" class="status-dot idle"></span>
            <div><p>TRAJECTORY STATUS</p><h2 id="status">Ready to solve</h2></div>
          </div>
          <p id="status-message" class="status-message">Enter a launch state and run the high-fidelity browser model.</p>
          <dl class="metrics">
            <div><dt>Closest return</dt><dd id="metric-miss">—</dd></div>
            <div><dt>Return time</dt><dd id="metric-time">—</dd></div>
            <div><dt>Minimum clearance</dt><dd id="metric-clearance">—</dd></div>
            <div><dt>Maximum altitude</dt><dd id="metric-altitude">—</dd></div>
            <div><dt>Impact location</dt><dd id="metric-impact">—</dd></div>
            <div><dt>Numerical bound</dt><dd id="metric-numerical">—</dd></div>
          </dl>
          <div id="domain-warning" class="domain-warning" hidden>
            <strong>Harmonic-domain warning</strong>
            Part of this arc lies below the conservative Brillouin sphere. Collision remains LDEM64-based, but high-degree gravity is a flagged downward continuation.
          </div>
        </section>
      </aside>
    </section>

    <section class="method-strip">
      <div><span>MODEL CHAIN</span><strong>GRGM atlas → rotating-frame dynamics → LDEM64 collision → continuous return search</strong></div>
      <div><span>INTERPRETATION</span><strong>A nominal return is not a catch: terminal speed and uncertainty are reported separately.</strong></div>
      <a href="https://svs.gsfc.nasa.gov/4720" target="_blank" rel="noreferrer">NASA texture credit ↗</a>
    </section>

    <footer>
      <p>Scientific visualization for “The Orbital Home Run.” Source models: NASA LOLA LDEM64 and GRAIL GRGM1200B.</p>
      <p>Not for navigation, landing, or safety-critical use.</p>
    </footer>
  </main>
`;

const form = document.querySelector<HTMLFormElement>("#solver-form")!;
const progress = document.querySelector<HTMLDivElement>("#progress")!;
const progressLabel = document.querySelector<HTMLSpanElement>("#progress-label")!;
const progressValue = document.querySelector<HTMLElement>("#progress-value")!;
const progressBar = document.querySelector<HTMLElement>("#progress-bar")!;
const runButton = form.querySelector<HTMLButtonElement>(".run-button")!;
const moonView = new MoonView(
  document.querySelector<HTMLElement>("#moon-view")!,
  baseUrl,
);
const worker = new Worker(new URL("./solver.worker.ts", import.meta.url), {
  type: "module",
});

const lunarImage = new Image();
lunarImage.src = `${baseUrl}assets/lroc_color_2k.jpg`;
lunarImage.onload = () => drawEmptyGroundTrack();

const formatDistance = (value: number | null): string => {
  if (value === null || !Number.isFinite(value)) return "—";
  if (Math.abs(value) >= 1000) return `${(value / 1000).toFixed(3)} km`;
  return `${value.toFixed(value < 10 ? 3 : 1)} m`;
};

const formatTime = (value: number | null): string => {
  if (value === null || !Number.isFinite(value)) return "—";
  const minutes = Math.floor(value / 60);
  const seconds = value - minutes * 60;
  return `${minutes} min ${seconds.toFixed(2)} s`;
};

const resizeCanvas = (canvas: HTMLCanvasElement): CanvasRenderingContext2D => {
  const ratio = Math.min(window.devicePixelRatio, 2);
  const width = Math.max(1, canvas.clientWidth);
  const height = Math.max(1, canvas.clientHeight);
  canvas.width = Math.round(width * ratio);
  canvas.height = Math.round(height * ratio);
  const context = canvas.getContext("2d")!;
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  return context;
};

const drawEmptyGroundTrack = (): void => {
  const canvas = document.querySelector<HTMLCanvasElement>("#ground-track")!;
  const context = resizeCanvas(canvas);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  context.fillStyle = "#111a20";
  context.fillRect(0, 0, width, height);
  if (lunarImage.complete) {
    context.globalAlpha = 0.56;
    context.drawImage(lunarImage, 0, 0, width, height);
    context.globalAlpha = 1;
  }
};

const drawGroundTrack = (result: SolverResult): void => {
  drawEmptyGroundTrack();
  const canvas = document.querySelector<HTMLCanvasElement>("#ground-track")!;
  const context = canvas.getContext("2d")!;
  const ratio = Math.min(window.devicePixelRatio, 2);
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  context.strokeStyle = "rgba(255,255,255,0.13)";
  context.lineWidth = 1;
  for (let lon = 0; lon <= 360; lon += 60) {
    const x = (lon / 360) * width;
    context.beginPath();
    context.moveTo(x, 0);
    context.lineTo(x, height);
    context.stroke();
  }
  for (let lat = -60; lat <= 60; lat += 30) {
    const y = ((90 - lat) / 180) * height;
    context.beginPath();
    context.moveTo(0, y);
    context.lineTo(width, y);
    context.stroke();
  }
  context.strokeStyle = result.status === "IMPACT" ? "#ff9c64" : "#79e4f5";
  context.lineWidth = 2;
  context.shadowBlur = 8;
  context.shadowColor = context.strokeStyle;
  let drawing = false;
  let previousX = 0;
  for (const point of result.points) {
    const x = (point.longitudeDegEast / 360) * width;
    const y = ((90 - point.latitudeDeg) / 180) * height;
    if (!drawing || Math.abs(x - previousX) > width / 2) {
      context.beginPath();
      context.moveTo(x, y);
      drawing = true;
    } else {
      context.lineTo(x, y);
      context.stroke();
    }
    previousX = x;
  }
  context.shadowBlur = 0;
  const first = result.points[0];
  context.fillStyle = "#7ff0b5";
  context.beginPath();
  context.arc(
    (first.longitudeDegEast / 360) * width,
    ((90 - first.latitudeDeg) / 180) * height,
    4,
    0,
    Math.PI * 2,
  );
  context.fill();
};

const drawClearance = (result: SolverResult): void => {
  const canvas = document.querySelector<HTMLCanvasElement>("#clearance-chart")!;
  const context = resizeCanvas(canvas);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  const margin = { left: 48, right: 12, top: 16, bottom: 27 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  context.fillStyle = "#0c151b";
  context.fillRect(0, 0, width, height);
  const values = result.points.map((point) => point.clearanceM / 1000);
  const maxTime = result.points.at(-1)!.timeS;
  const minimum = Math.min(0, ...values);
  const maximum = Math.max(0.001, ...values);
  const padding = Math.max(0.02, (maximum - minimum) * 0.08);
  const yMin = minimum - padding;
  const yMax = maximum + padding;
  const xOf = (time: number): number => margin.left + (time / maxTime) * plotWidth;
  const yOf = (value: number): number =>
    margin.top + ((yMax - value) / (yMax - yMin)) * plotHeight;
  context.strokeStyle = "rgba(255,255,255,0.12)";
  context.fillStyle = "#8496a0";
  context.font = "11px IBM Plex Mono, monospace";
  context.textAlign = "right";
  for (let index = 0; index <= 4; index += 1) {
    const value = yMin + ((yMax - yMin) * index) / 4;
    const y = yOf(value);
    context.beginPath();
    context.moveTo(margin.left, y);
    context.lineTo(width - margin.right, y);
    context.stroke();
    context.fillText(value.toFixed(1), margin.left - 7, y + 4);
  }
  if (yMin < 0 && yMax > 0) {
    context.strokeStyle = "#ff765f";
    context.setLineDash([4, 4]);
    context.beginPath();
    context.moveTo(margin.left, yOf(0));
    context.lineTo(width - margin.right, yOf(0));
    context.stroke();
    context.setLineDash([]);
  }
  context.strokeStyle = result.status === "IMPACT" ? "#ff9c64" : "#79e4f5";
  context.lineWidth = 2;
  context.beginPath();
  result.points.forEach((point, index) => {
    const x = xOf(point.timeS);
    const y = yOf(point.clearanceM / 1000);
    if (index === 0) context.moveTo(x, y);
    else context.lineTo(x, y);
  });
  context.stroke();
  context.fillStyle = "#8496a0";
  context.textAlign = "center";
  context.fillText("time after launch (min)", margin.left + plotWidth / 2, height - 7);
};

const updateResult = (result: SolverResult): void => {
  const status = document.querySelector<HTMLElement>("#status")!;
  const dot = document.querySelector<HTMLElement>("#status-dot")!;
  status.textContent = result.status.replaceAll("_", " ");
  dot.className = `status-dot ${result.status.toLowerCase()}`;
  document.querySelector<HTMLElement>("#status-message")!.textContent = result.message;
  document.querySelector<HTMLElement>("#degree-badge")!.textContent =
    `GRGM${result.gravityDegree}`;
  document.querySelector<HTMLElement>("#metric-miss")!.textContent = formatDistance(
    result.closestReturnDistanceM,
  );
  document.querySelector<HTMLElement>("#metric-time")!.textContent = formatTime(
    result.closestReturnTimeS,
  );
  document.querySelector<HTMLElement>("#metric-clearance")!.textContent =
    formatDistance(result.minimumClearanceM);
  document.querySelector<HTMLElement>("#metric-altitude")!.textContent =
    formatDistance(result.maximumAltitudeM);
  document.querySelector<HTMLElement>("#metric-impact")!.textContent =
    result.impactLatitudeDeg === null
      ? "—"
      : `${result.impactLatitudeDeg.toFixed(3)}°, ${result.impactLongitudeDegEast!.toFixed(3)}°E`;
  document.querySelector<HTMLElement>("#metric-numerical")!.textContent =
    formatDistance(result.numericalPositionUncertaintyM);
  document.querySelector<HTMLElement>("#domain-warning")!.hidden =
    result.entireArcOutsideBrillouinSphere;
  moonView.update(result);
  drawGroundTrack(result);
  drawClearance(result);
};

const readNumber = (data: FormData, key: string): number => {
  const value = Number(data.get(key));
  if (!Number.isFinite(value)) throw new Error(`Invalid ${key}`);
  return value;
};

form.addEventListener("submit", (event) => {
  event.preventDefault();
  const data = new FormData(form);
  const launch: LaunchInput = {
    latitudeDeg: readNumber(data, "latitude"),
    longitudeDegEast: readNumber(data, "longitude"),
    heightAboveTerrainM: readNumber(data, "height"),
    speedMS: readNumber(data, "speed"),
    elevationDeg: readNumber(data, "elevation"),
    azimuthDeg: readNumber(data, "azimuth"),
  };
  const config: SolverConfig = {
    durationS: readNumber(data, "duration"),
    stepS: readNumber(data, "step"),
    returnToleranceM: readNumber(data, "tolerance"),
    modelUncertaintyM: readNumber(data, "uncertainty"),
    minimumReturnTimeS: 1000,
    includeEarthTide: data.get("earth") === "on",
    includeSolarTide: data.get("sun") === "on",
  };
  progress.hidden = false;
  progressLabel.textContent = "Loading scientific tiles";
  progressValue.textContent = "0%";
  progressBar.style.width = "0%";
  runButton.disabled = true;
  const message: RunMessage = { type: "run", launch, config, baseUrl };
  worker.postMessage(message);
});

worker.onmessage = (
  event: MessageEvent<ProgressMessage | ResultMessage | ErrorMessage>,
) => {
  if (event.data.type === "progress") {
    const value = Math.round(event.data.fraction * 100);
    progressLabel.textContent = event.data.label;
    progressValue.textContent = `${value}%`;
    progressBar.style.width = `${value}%`;
  } else if (event.data.type === "result") {
    progress.hidden = true;
    runButton.disabled = false;
    updateResult(event.data.result);
  } else {
    progress.hidden = true;
    runButton.disabled = false;
    document.querySelector<HTMLElement>("#status")!.textContent = "MODEL ERROR";
    document.querySelector<HTMLElement>("#status-message")!.textContent =
      event.data.message;
    document.querySelector<HTMLElement>("#status-dot")!.className =
      "status-dot numerical_failure";
  }
};

document.querySelector<HTMLButtonElement>("#reset-button")!.addEventListener("click", () => {
  form.reset();
});

window.addEventListener("resize", () => {
  drawEmptyGroundTrack();
});
drawEmptyGroundTrack();
