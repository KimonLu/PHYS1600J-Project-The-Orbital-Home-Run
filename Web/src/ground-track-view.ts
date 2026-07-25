import type { SolverResult } from "./types";

export class GroundTrackView {
  private result: SolverResult | null = null;
  private scale = 1;
  private offsetX = 0;
  private offsetY = 0;
  private dragging = false;
  private pointerX = 0;
  private pointerY = 0;
  private axisLabel = "longitude / latitude";

  constructor(
    private readonly canvas: HTMLCanvasElement,
    private readonly image: HTMLImageElement,
  ) {
    canvas.tabIndex = 0;
    canvas.addEventListener("wheel", this.onWheel, { passive: false });
    canvas.addEventListener("pointerdown", this.onPointerDown);
    canvas.addEventListener("pointermove", this.onPointerMove);
    canvas.addEventListener("pointerup", this.onPointerUp);
    canvas.addEventListener("pointercancel", this.onPointerUp);
    canvas.addEventListener("dblclick", () => this.resetView());
    canvas.addEventListener("keydown", (event) => {
      if (event.key === "+" || event.key === "=") {
        event.preventDefault();
        this.zoomAt(1.25, this.canvas.clientWidth / 2, this.canvas.clientHeight / 2);
      } else if (event.key === "-") {
        event.preventDefault();
        this.zoomAt(0.8, this.canvas.clientWidth / 2, this.canvas.clientHeight / 2);
      } else if (event.key === "0") {
        event.preventDefault();
        this.resetView();
      }
    });
  }

  setAxisLabel(label: string): void {
    this.axisLabel = label;
    this.draw();
  }

  update(result: SolverResult): void {
    this.result = result;
    this.resetView();
  }

  resize(): void {
    this.clampOffsets();
    this.draw();
  }

  draw(): void {
    const context = this.resizeCanvas();
    const width = this.canvas.clientWidth;
    const height = this.canvas.clientHeight;
    context.fillStyle = "#111a20";
    context.fillRect(0, 0, width, height);
    context.save();
    context.translate(width / 2 + this.offsetX, height / 2 + this.offsetY);
    context.scale(this.scale, this.scale);
    context.translate(-width / 2, -height / 2);
    if (this.image.complete && this.image.naturalWidth > 0) {
      context.globalAlpha = 0.56;
      context.drawImage(this.image, 0, 0, width, height);
      context.globalAlpha = 1;
    }
    this.drawGrid(context, width, height);
    if (this.result) this.drawTrajectory(context, width, height, this.result);
    context.restore();
    context.fillStyle = "rgba(220,232,237,0.72)";
    context.font = '12px "IBM Plex Mono", monospace';
    context.textAlign = "left";
    context.fillText(`${this.scale.toFixed(2)}× · ${this.axisLabel}`, 12, height - 12);
  }

  private resizeCanvas(): CanvasRenderingContext2D {
    const ratio = Math.min(window.devicePixelRatio, 2);
    const width = Math.max(1, this.canvas.clientWidth);
    const height = Math.max(1, this.canvas.clientHeight);
    this.canvas.width = Math.round(width * ratio);
    this.canvas.height = Math.round(height * ratio);
    const context = this.canvas.getContext("2d")!;
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    return context;
  }

  private drawGrid(
    context: CanvasRenderingContext2D,
    width: number,
    height: number,
  ): void {
    context.strokeStyle = "rgba(255,255,255,0.13)";
    context.lineWidth = 1 / this.scale;
    for (let longitude = 0; longitude <= 360; longitude += 60) {
      const x = (longitude / 360) * width;
      context.beginPath();
      context.moveTo(x, 0);
      context.lineTo(x, height);
      context.stroke();
    }
    for (let latitude = -60; latitude <= 60; latitude += 30) {
      const y = ((90 - latitude) / 180) * height;
      context.beginPath();
      context.moveTo(0, y);
      context.lineTo(width, y);
      context.stroke();
    }
  }

  private drawTrajectory(
    context: CanvasRenderingContext2D,
    width: number,
    height: number,
    result: SolverResult,
  ): void {
    context.strokeStyle = result.status === "IMPACT" ? "#ff9c64" : "#79e4f5";
    context.lineWidth = 2 / this.scale;
    context.shadowBlur = 8 / this.scale;
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
      4 / this.scale,
      0,
      Math.PI * 2,
    );
    context.fill();
  }

  private resetView(): void {
    this.scale = 1;
    this.offsetX = 0;
    this.offsetY = 0;
    this.draw();
  }

  private zoomAt(factor: number, x: number, y: number): void {
    const oldScale = this.scale;
    const newScale = Math.min(8, Math.max(1, oldScale * factor));
    if (newScale === oldScale) return;
    const centerX = this.canvas.clientWidth / 2;
    const centerY = this.canvas.clientHeight / 2;
    const worldX = (x - centerX - this.offsetX) / oldScale;
    const worldY = (y - centerY - this.offsetY) / oldScale;
    this.scale = newScale;
    this.offsetX = x - centerX - worldX * newScale;
    this.offsetY = y - centerY - worldY * newScale;
    this.clampOffsets();
    this.draw();
  }

  private clampOffsets(): void {
    const maxX = (this.canvas.clientWidth * (this.scale - 1)) / 2;
    const maxY = (this.canvas.clientHeight * (this.scale - 1)) / 2;
    this.offsetX = Math.min(maxX, Math.max(-maxX, this.offsetX));
    this.offsetY = Math.min(maxY, Math.max(-maxY, this.offsetY));
  }

  private onWheel = (event: WheelEvent): void => {
    event.preventDefault();
    const bounds = this.canvas.getBoundingClientRect();
    this.zoomAt(
      Math.exp(-event.deltaY * 0.0015),
      event.clientX - bounds.left,
      event.clientY - bounds.top,
    );
  };

  private onPointerDown = (event: PointerEvent): void => {
    if (event.button !== 0) return;
    this.dragging = true;
    this.pointerX = event.clientX;
    this.pointerY = event.clientY;
    this.canvas.setPointerCapture(event.pointerId);
    this.canvas.classList.add("is-dragging");
  };

  private onPointerMove = (event: PointerEvent): void => {
    if (!this.dragging) return;
    this.offsetX += event.clientX - this.pointerX;
    this.offsetY += event.clientY - this.pointerY;
    this.pointerX = event.clientX;
    this.pointerY = event.clientY;
    this.clampOffsets();
    this.draw();
  };

  private onPointerUp = (event: PointerEvent): void => {
    if (!this.dragging) return;
    this.dragging = false;
    if (this.canvas.hasPointerCapture(event.pointerId)) {
      this.canvas.releasePointerCapture(event.pointerId);
    }
    this.canvas.classList.remove("is-dragging");
  };
}
