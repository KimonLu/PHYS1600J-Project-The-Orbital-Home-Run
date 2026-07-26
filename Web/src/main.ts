import "./style.css";
import { trajectoryExamples } from "./examples";
import { GroundTrackView } from "./ground-track-view";
import {
  statusMessage,
  statusText,
  text,
  type Language,
  type TranslationKey,
} from "./i18n";
import { MoonView, type PlaybackState } from "./moon-view";
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
const MU_MOON_M3_S2 = 4.902800118e12;
const MEAN_MOON_RADIUS_M = 1_737_400;
let language: Language =
  window.localStorage.getItem("orbital-home-run-language") === "zh" ? "zh" : "en";
let latestResult: SolverResult | null = null;
let latestPlayback: PlaybackState = {
  playing: false,
  timeS: 0,
  durationS: 0,
  speed: 100,
};

document.querySelector<HTMLDivElement>("#app")!.innerHTML = `
  <main class="shell">
    <header class="masthead">
      <div>
        <p class="eyebrow" data-i18n="eyebrow"></p>
        <h1 data-i18n="title"></h1>
        <p class="lede" data-i18n="lede"></p>
      </div>
      <div class="masthead-tools">
        <div class="language-switch" role="group" aria-label="Language">
          <button type="button" data-language="en">EN</button>
          <button type="button" data-language="zh">中文</button>
        </div>
        <div class="model-badges" data-i18n-aria="scientificModel">
          <a href="https://pds-geosciences.wustl.edu/missions/lro/lola.htm" target="_blank" rel="noreferrer">
            <span>LOLA</span><strong>LDEM64 ↗</strong>
          </a>
          <a href="https://pgda.gsfc.nasa.gov/products/75" target="_blank" rel="noreferrer">
            <span>GRAIL</span><strong id="degree-badge">GRGM600 ↗</strong>
          </a>
          <a href="https://www.gnu.org/software/gsl/doc/html/ode-initval.html" target="_blank" rel="noreferrer">
            <span data-i18n="integrator"></span><strong>RK4 step-doubling ↗</strong>
          </a>
          <a href="https://svs.gsfc.nasa.gov/4720" target="_blank" rel="noreferrer">
            <span data-i18n="surfaceTexture"></span><strong>NASA SVS Moon Kit ↗</strong>
          </a>
        </div>
      </div>
    </header>

    <section class="workspace">
      <div class="visual-column">
        <article class="panel moon-panel">
          <div class="panel-heading overlay">
            <div>
              <span class="section-index">01</span>
              <h2 data-i18n="trajectory3d"></h2>
            </div>
            <p data-i18n="trajectoryHint"></p>
          </div>
          <div id="moon-view" class="moon-view" aria-label="Interactive 3D Moon and trajectory"></div>
          <div class="playback-controls">
            <button type="button" id="playback-toggle" disabled data-i18n="pause"></button>
            <button type="button" id="playback-restart" disabled data-i18n="restart"></button>
            <label>
              <span data-i18n="playbackSpeed"></span>
              <select id="playback-speed">
                <option value="1">1×</option>
                <option value="10">10×</option>
                <option value="100" selected>100×</option>
                <option value="500">500×</option>
                <option value="1000">1000×</option>
              </select>
            </label>
            <output id="playback-time">00:00 / 00:00</output>
          </div>
          <div class="legend overlay-legend">
            <span><i class="dot launch"></i><span data-i18n="launch"></span></span>
            <span><i class="line orbit"></i><span data-i18n="trajectory"></span></span>
            <span><i class="dot ball"></i><span data-i18n="ball"></span></span>
            <span><i class="line axis"></i><span data-i18n="spinAxis"></span></span>
            <span><i class="dot event"></i><span data-i18n="event"></span></span>
          </div>
        </article>

        <div class="analysis-grid">
          <article class="panel chart-panel">
            <div class="panel-heading">
              <div><span class="section-index">02</span><h2 data-i18n="groundTrack"></h2></div>
              <span class="unit" data-i18n="groundTrackUnit"></span>
            </div>
            <canvas id="ground-track" aria-label="Zoomable ground track on lunar map"></canvas>
          </article>
          <article class="panel chart-panel">
            <div class="panel-heading">
              <div><span class="section-index">03</span><h2 data-i18n="clearance"></h2></div>
              <span class="unit" data-i18n="clearanceUnit"></span>
            </div>
            <canvas id="clearance-chart" aria-label="Terrain clearance versus time"></canvas>
          </article>
        </div>
      </div>

      <aside class="control-column">
        <form id="solver-form" class="panel controls">
          <div class="panel-heading">
            <div><span class="section-index" data-i18n="input"></span><h2 data-i18n="launchState"></h2></div>
            <button type="button" id="reset-button" class="text-button" data-i18n="reset"></button>
          </div>

          <fieldset class="example-fieldset">
            <legend data-i18n="example"></legend>
            <select id="example-select" aria-label="Example trajectory">
              <option value="" data-i18n="customExample"></option>
            </select>
          </fieldset>

          <fieldset>
            <legend data-i18n="launchSite"></legend>
            <div class="field-grid">
              <label><span class="label-text" data-i18n="latitude"></span><span>deg</span><input name="latitude" type="number" min="-90" max="90" step="0.0000001" value="5.4296875" required></label>
              <label><span class="label-text" data-i18n="longitude"></span><span>deg</span><input name="longitude" type="number" step="0.0000001" value="201.3671875" required></label>
              <label class="wide"><span class="label-text" data-i18n="height"></span><span>m</span><input name="height" type="number" min="0.037" step="any" value="19243" required></label>
            </div>
          </fieldset>

          <fieldset>
            <legend data-i18n="velocity"></legend>
            <div class="field-grid">
              <label class="wide"><span class="label-text" data-i18n="speed"></span><span>m s⁻¹</span><input name="speed" type="number" min="0" step="any" value="1661.4559" required></label>
              <label><span class="label-text" data-i18n="elevation"></span><span>deg</span><input name="elevation" type="number" min="-90" max="90" step="any" value="0.531603" required></label>
              <label><span class="label-text" data-i18n="azimuth"></span><span>deg</span><input name="azimuth" type="number" step="any" value="88.383067" required></label>
            </div>
            <p class="field-note" data-i18n="azimuthNote"></p>
            <p class="field-note speed-reference" id="speed-reference"></p>
          </fieldset>

          <fieldset>
            <legend data-i18n="decision"></legend>
            <div class="field-grid">
              <label><span class="label-text" data-i18n="duration"></span><span>s</span><input name="duration" type="number" min="1000" step="any" value="6686.34" required></label>
              <label><span class="label-text" data-i18n="step"></span><span>s</span><input name="step" type="number" min="0.25" max="20" step="any" value="4" required></label>
              <label><span class="label-text" data-i18n="tolerance"></span><span>m</span><input name="tolerance" type="number" min="0" step="any" value="10" required></label>
              <label><span class="label-text" data-i18n="uncertainty"></span><span>m</span><input name="uncertainty" type="number" min="0" step="any" value="10" required></label>
            </div>
            <div class="switch-row">
              <label class="switch"><input name="earth" type="checkbox" checked><span></span><span class="switch-text" data-i18n="earthTide"></span></label>
              <label class="switch"><input name="sun" type="checkbox" checked><span></span><span class="switch-text" data-i18n="solarTide"></span></label>
            </div>
            <p class="field-note" data-i18n="thirdBodyNote"></p>
          </fieldset>

          <button class="run-button" type="submit">
            <span data-i18n="run"></span>
            <small data-i18n="runDetail"></small>
          </button>
          <div id="progress" class="progress" hidden>
            <div><span id="progress-label" data-i18n="preparing"></span><strong id="progress-value">0%</strong></div>
            <div class="progress-track"><i id="progress-bar"></i></div>
          </div>
        </form>

        <section class="panel result-panel" aria-live="polite">
          <div class="result-heading">
            <span id="status-dot" class="status-dot idle"></span>
            <div><p data-i18n="trajectoryStatus"></p><h2 id="status" data-i18n="ready"></h2></div>
          </div>
          <p id="status-message" class="status-message" data-i18n="readyMessage"></p>
          <dl class="metrics">
            <div><dt data-i18n="closestReturn"></dt><dd id="metric-miss">—</dd></div>
            <div><dt data-i18n="returnTime"></dt><dd id="metric-time">—</dd></div>
            <div><dt data-i18n="minimumClearance"></dt><dd id="metric-clearance">—</dd></div>
            <div><dt data-i18n="maximumAltitude"></dt><dd id="metric-altitude">—</dd></div>
            <div><dt data-i18n="impactLocation"></dt><dd id="metric-impact">—</dd></div>
            <div><dt data-i18n="numericalBound"></dt><dd id="metric-numerical">—</dd></div>
          </dl>
          <div id="domain-warning" class="domain-warning" hidden>
            <strong data-i18n="domainTitle"></strong>
            <span data-i18n="domainMessage"></span>
          </div>
        </section>
      </aside>
    </section>

    <section class="method-strip">
      <div><span data-i18n="modelChain"></span><strong data-i18n="modelChainText"></strong></div>
      <a href="https://github.com/KimonLu/PHYS1600J-Project-The-Orbital-Home-Run" target="_blank" rel="noreferrer" data-i18n="githubRepository"></a>
    </section>
  </main>
`;

const form = document.querySelector<HTMLFormElement>("#solver-form")!;
const progress = document.querySelector<HTMLDivElement>("#progress")!;
const progressLabel = document.querySelector<HTMLSpanElement>("#progress-label")!;
const progressValue = document.querySelector<HTMLElement>("#progress-value")!;
const progressBar = document.querySelector<HTMLElement>("#progress-bar")!;
const runButton = form.querySelector<HTMLButtonElement>(".run-button")!;
const playbackToggle = document.querySelector<HTMLButtonElement>("#playback-toggle")!;
const playbackRestart = document.querySelector<HTMLButtonElement>("#playback-restart")!;
const playbackSpeed = document.querySelector<HTMLSelectElement>("#playback-speed")!;
const playbackTime = document.querySelector<HTMLOutputElement>("#playback-time")!;
const exampleSelect = document.querySelector<HTMLSelectElement>("#example-select")!;
const heightInput = form.elements.namedItem("height") as HTMLInputElement;

const formatClock = (seconds: number): string => {
  const safeSeconds = Math.max(0, seconds);
  const hours = Math.floor(safeSeconds / 3600);
  const minutes = Math.floor((safeSeconds - hours * 3600) / 60);
  const remainder = Math.floor(safeSeconds % 60);
  return `${hours > 0 ? `${hours.toString().padStart(2, "0")}:` : ""}${minutes
    .toString()
    .padStart(2, "0")}:${remainder.toString().padStart(2, "0")}`;
};

const updatePlayback = (state: PlaybackState): void => {
  latestPlayback = state;
  playbackToggle.disabled = state.durationS <= 0;
  playbackRestart.disabled = state.durationS <= 0;
  playbackToggle.textContent = text(language, state.playing ? "pause" : "play");
  playbackTime.textContent = `${formatClock(state.timeS)} / ${formatClock(state.durationS)}`;
};

const moonView = new MoonView(
  document.querySelector<HTMLElement>("#moon-view")!,
  baseUrl,
  updatePlayback,
);
const worker = new Worker(new URL("./solver.worker.ts", import.meta.url), {
  type: "module",
});

const lunarImage = new Image();
lunarImage.src = `${baseUrl}assets/lroc_color_2k.jpg`;
const groundTrackView = new GroundTrackView(
  document.querySelector<HTMLCanvasElement>("#ground-track")!,
  lunarImage,
);
lunarImage.onload = () => groundTrackView.draw();

const formatDistance = (value: number | null): string => {
  if (value === null || !Number.isFinite(value)) return "—";
  if (Math.abs(value) >= 1000) return `${(value / 1000).toFixed(3)} km`;
  return `${value.toFixed(Math.abs(value) < 10 ? 3 : 1)} m`;
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

const drawClearance = (result: SolverResult): void => {
  const canvas = document.querySelector<HTMLCanvasElement>("#clearance-chart")!;
  const context = resizeCanvas(canvas);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  const margin = { left: 52, right: 14, top: 18, bottom: 31 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  context.fillStyle = "#0c151b";
  context.fillRect(0, 0, width, height);
  const values = result.points.map((point) => point.clearanceM / 1000);
  const maxTime = Math.max(1, result.points.at(-1)!.timeS);
  const minimum = Math.min(0, ...values);
  const maximum = Math.max(0.001, ...values);
  const padding = Math.max(0.02, (maximum - minimum) * 0.08);
  const yMin = minimum - padding;
  const yMax = maximum + padding;
  const xOf = (time: number): number => margin.left + (time / maxTime) * plotWidth;
  const yOf = (value: number): number =>
    margin.top + ((yMax - value) / (yMax - yMin)) * plotHeight;
  context.strokeStyle = "rgba(255,255,255,0.12)";
  context.fillStyle = "#91a5ae";
  context.font = '12px "IBM Plex Mono", monospace';
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
  context.fillStyle = "#91a5ae";
  context.textAlign = "center";
  context.fillText(text(language, "timeAxis"), margin.left + plotWidth / 2, height - 9);
};

const renderResultSummary = (result: SolverResult): void => {
  document.querySelector<HTMLElement>("#status")!.textContent = statusText(
    language,
    result.status,
  );
  document.querySelector<HTMLElement>("#status-dot")!.className =
    `status-dot ${result.status.toLowerCase()}`;
  document.querySelector<HTMLElement>("#status-message")!.textContent = statusMessage(
    language,
    result.status,
    !result.entireArcOutsideBrillouinSphere,
  );
  document.querySelector<HTMLElement>("#degree-badge")!.textContent =
    `GRGM${result.gravityDegree} ↗`;
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
};

const updateResult = (result: SolverResult): void => {
  latestResult = result;
  renderResultSummary(result);
  moonView.update(result);
  groundTrackView.update(result);
  drawClearance(result);
};

const updateSpeedReference = (): void => {
  const heightM = Math.max(0, Number(heightInput.value) || 0);
  const radiusM = MEAN_MOON_RADIUS_M + heightM;
  const circular = Math.sqrt(MU_MOON_M3_S2 / radiusM) / 1000;
  const escape = Math.sqrt((2 * MU_MOON_M3_S2) / radiusM) / 1000;
  document.querySelector<HTMLElement>("#speed-reference")!.textContent =
    language === "en"
      ? `At mean datum + entered height: circular (first cosmic) ${circular.toFixed(3)} km/s · escape ${escape.toFixed(3)} km/s.`
      : `按平均月面半径加当前高度估算：第一宇宙速度 ${circular.toFixed(3)} km/s · 第二宇宙速度 ${escape.toFixed(3)} km/s。`;
};

const populateExamples = (): void => {
  const selected = exampleSelect.value;
  exampleSelect.querySelectorAll("option:not(:first-child)").forEach((option) => option.remove());
  for (const example of trajectoryExamples) {
    const option = document.createElement("option");
    option.value = example.status;
    option.textContent = statusText(language, example.status);
    exampleSelect.append(option);
  }
  exampleSelect.value = selected;
};

const applyLanguage = (nextLanguage: Language): void => {
  language = nextLanguage;
  window.localStorage.setItem("orbital-home-run-language", language);
  document.documentElement.lang = language === "zh" ? "zh-CN" : "en";
  document.title =
    language === "zh"
      ? "月球轨道本垒打 — 月球轨迹实验室"
      : "Orbital Home Run — Lunar Trajectory Laboratory";
  document.querySelectorAll<HTMLElement>("[data-i18n]").forEach((element) => {
    element.textContent = text(language, element.dataset.i18n as TranslationKey);
  });
  document.querySelectorAll<HTMLElement>("[data-i18n-aria]").forEach((element) => {
    element.setAttribute(
      "aria-label",
      text(language, element.dataset.i18nAria as TranslationKey),
    );
  });
  document.querySelectorAll<HTMLButtonElement>("[data-language]").forEach((button) => {
    button.classList.toggle("active", button.dataset.language === language);
    button.setAttribute("aria-pressed", String(button.dataset.language === language));
  });
  populateExamples();
  updateSpeedReference();
  groundTrackView.setAxisLabel(language === "en" ? "longitude / latitude" : "经度 / 纬度");
  moonView.setAxisLabels(
    text(language, "northPole"),
    text(language, "southPole"),
    text(language, "rotationDirection"),
  );
  if (latestResult) {
    renderResultSummary(latestResult);
    drawClearance(latestResult);
  }
  updatePlayback(latestPlayback);
};

const setNamedValue = (name: string, value: number): void => {
  (form.elements.namedItem(name) as HTMLInputElement).value = String(value);
};

const applyExample = (status: string): void => {
  const example = trajectoryExamples.find((item) => item.status === status);
  if (!example) return;
  setNamedValue("latitude", example.launch.latitudeDeg);
  setNamedValue("longitude", example.launch.longitudeDegEast);
  setNamedValue("height", example.launch.heightAboveTerrainM);
  setNamedValue("speed", example.launch.speedMS);
  setNamedValue("elevation", example.launch.elevationDeg);
  setNamedValue("azimuth", example.launch.azimuthDeg);
  setNamedValue("duration", example.config.durationS);
  setNamedValue("step", example.config.stepS);
  setNamedValue("tolerance", example.config.returnToleranceM);
  setNamedValue(
    "uncertainty",
    example.config.assumedTotalPositionErrorBoundM,
  );
  (form.elements.namedItem("earth") as HTMLInputElement).checked =
    example.config.includeEarthTide;
  (form.elements.namedItem("sun") as HTMLInputElement).checked =
    example.config.includeSolarTide;
  updateSpeedReference();
};

const readNumber = (data: FormData, key: string): number => {
  const value = Number(data.get(key));
  if (!Number.isFinite(value)) throw new Error(`Invalid ${key}`);
  return value;
};

form.addEventListener("input", (event) => {
  if (event.target !== exampleSelect) exampleSelect.value = "";
});
heightInput.addEventListener("input", updateSpeedReference);
exampleSelect.addEventListener("change", () => applyExample(exampleSelect.value));
playbackToggle.addEventListener("click", () => moonView.togglePlayback());
playbackRestart.addEventListener("click", () => moonView.restartPlayback());
playbackSpeed.addEventListener("change", () =>
  moonView.setPlaybackSpeed(Number(playbackSpeed.value)),
);
document.querySelectorAll<HTMLButtonElement>("[data-language]").forEach((button) => {
  button.addEventListener("click", () => applyLanguage(button.dataset.language as Language));
});

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
    assumedTotalPositionErrorBoundM: readNumber(data, "uncertainty"),
    minimumReturnTimeS: 1000,
    includeEarthTide: data.get("earth") === "on",
    includeSolarTide: data.get("sun") === "on",
  };
  progress.hidden = false;
  progressLabel.textContent = text(language, "loadingTiles");
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
    progressLabel.textContent =
      language === "en" ? event.data.label : `正在计算轨迹 ${value}%`;
    progressValue.textContent = `${value}%`;
    progressBar.style.width = `${value}%`;
  } else if (event.data.type === "result") {
    progress.hidden = true;
    runButton.disabled = false;
    updateResult(event.data.result);
  } else {
    progress.hidden = true;
    runButton.disabled = false;
    document.querySelector<HTMLElement>("#status")!.textContent = text(
      language,
      "modelError",
    );
    document.querySelector<HTMLElement>("#status-message")!.textContent =
      event.data.message;
    document.querySelector<HTMLElement>("#status-dot")!.className =
      "status-dot numerical_failure";
  }
};

document.querySelector<HTMLButtonElement>("#reset-button")!.addEventListener("click", () => {
  form.reset();
  exampleSelect.value = "";
  updateSpeedReference();
});

window.addEventListener("resize", () => {
  groundTrackView.resize();
  if (latestResult) drawClearance(latestResult);
});

applyLanguage(language);
groundTrackView.draw();
