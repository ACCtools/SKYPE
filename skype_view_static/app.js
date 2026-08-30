"use strict";

const canvas = document.getElementById("depth-canvas");
const context = canvas.getContext("2d");
const canvasWrap = document.getElementById("canvas-wrap");
const tooltip = document.getElementById("tooltip");
const loading = document.getElementById("loading");
const sampleName = document.getElementById("sample-name");
const datasetSummary = document.getElementById("dataset-summary");
const zoomValue = document.getElementById("zoom-value");
const resetButton = document.getElementById("reset-view");
const lineWidthInput = document.getElementById("line-width");
const lineWidthValue = document.getElementById("line-width-value");

const TWO_PI = Math.PI * 2;
const SECTOR_GAP = 3 * Math.PI / 180;
const MIN_ZOOM = 1;
const MAX_ZOOM = 64;

const state = {
  data: null,
  sectors: [],
  width: 0,
  height: 0,
  radius: 0,
  lineWidthScale: Number(lineWidthInput.value),
  camera: { scale: 1, offsetX: 0, offsetY: 0 },
  hover: null,
  dragging: false,
  dragX: 0,
  dragY: 0,
  frameRequested: false,
};

function formatInteger(value) {
  return Math.round(value).toLocaleString("en-US");
}

function formatDepth(value) {
  if (!Number.isFinite(value)) return "N/A";
  return value.toLocaleString("en-US", {
    minimumFractionDigits: 2,
    maximumFractionDigits: 2,
  });
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function buildSectors() {
  const totalLength = state.data.chromosomes.reduce(
    (sum, chromosome) => sum + chromosome.length,
    0,
  );
  const availableAngle = TWO_PI - SECTOR_GAP * state.data.chromosomes.length;
  let angle = -Math.PI / 2 + SECTOR_GAP / 2;
  state.sectors = state.data.chromosomes.map((chromosome) => {
    const span = availableAngle * chromosome.length / totalLength;
    const sector = {
      ...chromosome,
      startAngle: angle,
      endAngle: angle + span,
      span,
      indexByStart: new Map(
        chromosome.starts.map((start, index) => [start, index]),
      ),
    };
    angle += span + SECTOR_GAP;
    return sector;
  });
}

function requestRender() {
  if (state.frameRequested) return;
  state.frameRequested = true;
  requestAnimationFrame(() => {
    state.frameRequested = false;
    render();
  });
}

function resizeCanvas() {
  const rect = canvasWrap.getBoundingClientRect();
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  state.width = Math.max(1, rect.width);
  state.height = Math.max(1, rect.height);
  state.radius = Math.min(state.width, state.height) * 0.43;
  canvas.width = Math.round(state.width * dpr);
  canvas.height = Math.round(state.height * dpr);
  canvas.style.width = `${state.width}px`;
  canvas.style.height = `${state.height}px`;
  context.setTransform(dpr, 0, 0, dpr, 0, 0);
  requestRender();
}

function pointOnCircle(angle, radius) {
  return [Math.cos(angle) * radius, Math.sin(angle) * radius];
}

function uprightTangentRotation(angle) {
  let rotation = angle + Math.PI / 2;
  while (rotation > Math.PI) rotation -= TWO_PI;
  while (rotation <= -Math.PI) rotation += TWO_PI;
  if (rotation > Math.PI / 2) rotation -= Math.PI;
  if (rotation < -Math.PI / 2) rotation += Math.PI;
  return rotation;
}

function depthRadius(value, innerRadius, outerRadius) {
  const fraction = Math.max(0, Math.min(value, state.data.scale_max))
    / state.data.scale_max;
  return innerRadius + fraction * (outerRadius - innerRadius);
}

function drawSectorFrame(sector, innerRadius, outerRadius, scale) {
  context.strokeStyle = "#6f7887";
  context.lineWidth = 0.85 / scale;
  context.beginPath();
  context.arc(0, 0, innerRadius, sector.startAngle, sector.endAngle);
  context.stroke();
  context.beginPath();
  context.arc(0, 0, outerRadius, sector.startAngle, sector.endAngle);
  context.stroke();
  const innerStart = pointOnCircle(sector.startAngle, innerRadius);
  const outerStart = pointOnCircle(sector.startAngle, outerRadius);
  const innerEnd = pointOnCircle(sector.endAngle, innerRadius);
  const outerEnd = pointOnCircle(sector.endAngle, outerRadius);
  context.beginPath();
  context.moveTo(innerStart[0], innerStart[1]);
  context.lineTo(outerStart[0], outerStart[1]);
  context.moveTo(innerEnd[0], innerEnd[1]);
  context.lineTo(outerEnd[0], outerEnd[1]);
  context.stroke();

  const thresholds = [
    [0.5 * state.data.median_depth, "#3977d4"],
    [2 * state.data.median_depth, "#dc4b56"],
  ];
  for (const [depth, color] of thresholds) {
    context.strokeStyle = color;
    context.globalAlpha = 0.72;
    context.lineWidth = 0.75 / scale;
    context.beginPath();
    context.arc(
      0,
      0,
      depthRadius(depth, innerRadius, outerRadius),
      sector.startAngle,
      sector.endAngle,
    );
    context.stroke();
  }
  context.globalAlpha = 1;
}

function drawDepthLine(sector, values, color, width, innerRadius, outerRadius, scale) {
  context.strokeStyle = color;
  context.lineWidth = width * state.lineWidthScale / scale;
  context.lineJoin = "round";
  context.lineCap = "round";
  context.beginPath();
  let previousStart = null;
  for (let index = 0; index < sector.starts.length; index += 1) {
    const start = sector.starts[index];
    const angle = sector.startAngle
      + sector.span * ((start - 1) / sector.length);
    const radius = depthRadius(values[index], innerRadius, outerRadius);
    const [x, y] = pointOnCircle(angle, radius);
    if (previousStart === null || start - previousStart !== state.data.bin_size) {
      context.moveTo(x, y);
    } else {
      context.lineTo(x, y);
    }
    previousStart = start;
  }
  context.stroke();
}

function drawTicksAndLabel(sector, outerRadius, scale) {
  const major = 25_000_000;
  const minor = 5_000_000;
  for (let position = 0; position <= sector.length; position += minor) {
    const angle = sector.startAngle + sector.span * position / sector.length;
    const isMajor = position % major === 0;
    const tickEnd = outerRadius + (isMajor ? 6 : 3) / scale;
    const [x1, y1] = pointOnCircle(angle, outerRadius + 1 / scale);
    const [x2, y2] = pointOnCircle(angle, tickEnd);
    context.strokeStyle = "#596273";
    context.lineWidth = 0.7 / scale;
    context.beginPath();
    context.moveTo(x1, y1);
    context.lineTo(x2, y2);
    context.stroke();

    if (isMajor && scale <= 8) {
      const [tx, ty] = pointOnCircle(angle, outerRadius + 10 / scale);
      context.save();
      context.translate(tx, ty);
      context.rotate(uprightTangentRotation(angle));
      context.fillStyle = "#657086";
      context.font = `${8.5 / scale}px ui-sans-serif, system-ui, sans-serif`;
      context.textAlign = "center";
      context.textBaseline = "middle";
      context.fillText(String(position / 1_000_000), 0, 0);
      context.restore();
    }
  }

  const middle = (sector.startAngle + sector.endAngle) / 2;
  const [labelX, labelY] = pointOnCircle(middle, outerRadius + 29 / scale);
  context.save();
  context.translate(labelX, labelY);
  context.rotate(uprightTangentRotation(middle));
  context.fillStyle = "#202a3d";
  context.font = `600 ${11 / scale}px ui-sans-serif, system-ui, sans-serif`;
  context.textAlign = "center";
  context.textBaseline = "middle";
  context.fillText(sector.name, 0, 0);
  context.restore();
}

function drawHover(innerRadius, outerRadius, scale) {
  if (!state.hover) return;
  const { sector, binStart, binCenter, index, available } = state.hover;
  const binEnd = Math.min(binStart + state.data.bin_size - 1, sector.length);
  const angleStart = sector.startAngle
    + sector.span * ((binStart - 1) / sector.length);
  const angleEnd = sector.startAngle
    + sector.span * (binEnd / sector.length);
  context.fillStyle = available
    ? "rgba(246, 180, 65, 0.27)"
    : "rgba(121, 132, 151, 0.18)";
  context.beginPath();
  context.arc(0, 0, outerRadius, angleStart, angleEnd);
  context.arc(0, 0, innerRadius, angleEnd, angleStart, true);
  context.closePath();
  context.fill();

  const cursorAngle = sector.startAngle
    + sector.span * ((binCenter - 1) / sector.length);
  const [x1, y1] = pointOnCircle(cursorAngle, innerRadius - 4 / scale);
  const [x2, y2] = pointOnCircle(cursorAngle, outerRadius + 4 / scale);
  context.strokeStyle = available ? "#e59b17" : "#7b8799";
  context.lineWidth = 1.2 / scale;
  context.beginPath();
  context.moveTo(x1, y1);
  context.lineTo(x2, y2);
  context.stroke();

  if (available) {
    const markerDepths = [
      [sector.reference_depth[index], "#111827"],
      [sector.predicted_depth[index], "#159447"],
    ];
    for (const [depth, color] of markerDepths) {
      const markerRadius = depthRadius(depth, innerRadius, outerRadius);
      const [markerX, markerY] = pointOnCircle(cursorAngle, markerRadius);
      context.fillStyle = color;
      context.strokeStyle = "#ffffff";
      context.lineWidth = 1 / scale;
      context.beginPath();
      context.arc(markerX, markerY, 3.2 / scale, 0, TWO_PI);
      context.fill();
      context.stroke();
    }
  }
}

function render() {
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  context.setTransform(dpr, 0, 0, dpr, 0, 0);
  context.clearRect(0, 0, state.width, state.height);
  if (!state.data || !state.width || !state.height) return;

  const { scale, offsetX, offsetY } = state.camera;
  const innerRadius = state.radius * 0.61;
  const outerRadius = state.radius * 0.94;
  context.save();
  context.translate(state.width / 2 + offsetX, state.height / 2 + offsetY);
  context.scale(scale, scale);

  for (const sector of state.sectors) {
    drawSectorFrame(sector, innerRadius, outerRadius, scale);
    drawDepthLine(
      sector,
      sector.reference_depth,
      "#111827",
      0.85,
      innerRadius,
      outerRadius,
      scale,
    );
    drawDepthLine(
      sector,
      sector.predicted_depth,
      "#159447",
      1.05,
      innerRadius,
      outerRadius,
      scale,
    );
    drawTicksAndLabel(sector, outerRadius, scale);
  }
  drawHover(innerRadius, outerRadius, scale);
  context.restore();
}

function normalizedAngle(angle, minimum) {
  while (angle < minimum) angle += TWO_PI;
  while (angle >= minimum + TWO_PI) angle -= TWO_PI;
  return angle;
}

function hitTest(screenX, screenY) {
  if (!state.data) return null;
  const { scale, offsetX, offsetY } = state.camera;
  const worldX = (screenX - state.width / 2 - offsetX) / scale;
  const worldY = (screenY - state.height / 2 - offsetY) / scale;
  const radius = Math.hypot(worldX, worldY);
  const innerRadius = state.radius * 0.61;
  const outerRadius = state.radius * 0.94;
  if (radius < innerRadius || radius > outerRadius) return null;

  let angle = Math.atan2(worldY, worldX);
  angle = normalizedAngle(angle, state.sectors[0].startAngle);
  const sector = state.sectors.find(
    (candidate) => angle >= candidate.startAngle && angle <= candidate.endAngle,
  );
  if (!sector) return null;

  const fraction = Math.max(
    0,
    Math.min(1, (angle - sector.startAngle) / sector.span),
  );
  const position = Math.min(
    sector.length,
    Math.floor(fraction * sector.length) + 1,
  );
  const binStart = Math.floor((position - 1) / state.data.bin_size)
    * state.data.bin_size + 1;
  const binEnd = Math.min(binStart + state.data.bin_size - 1, sector.length);
  const binCenter = Math.floor((binStart + binEnd) / 2);
  const index = sector.indexByStart.get(binStart);
  return {
    sector,
    binStart,
    binCenter,
    index,
    available: index !== undefined,
  };
}

function tooltipRows(rows) {
  return rows.map(([label, value]) => (
    `<div class="row"><span class="label">${escapeHtml(label)}</span>`
    + `<span>${escapeHtml(value)}</span></div>`
  )).join("");
}

function updateTooltip(hit, screenX, screenY) {
  if (!hit) {
    tooltip.hidden = true;
    return;
  }
  const { sector, binStart, binCenter, index, available } = hit;
  const binEnd = Math.min(binStart + state.data.bin_size - 1, sector.length);
  let html = `<strong>${escapeHtml(sector.name)}:`
    + `${formatInteger(binStart)}–${formatInteger(binEnd)}</strong>`;
  html += tooltipRows([
    ["Selected bin", `${formatInteger(binStart)}–${formatInteger(binEnd)}`],
    ["Bin center", `${sector.name}:${formatInteger(binCenter)}`],
  ]);
  if (available) {
    const reference = sector.reference_depth[index];
    const predicted = sector.predicted_depth[index];
    const predictedRaw = sector.predicted_depth_raw[index];
    html += tooltipRows([
      ["Reference depth", formatDepth(reference)],
      ["Predicted depth", formatDepth(predicted)],
      ["Difference", formatDepth(predicted - reference)],
      ["Reference CN", formatDepth(reference / state.data.median_depth * 2)],
      ["Predicted CN", formatDepth(predicted / state.data.median_depth * 2)],
    ]);
    if (predictedRaw !== predicted) {
      html += tooltipRows([["Raw prediction", formatDepth(predictedRaw)]]);
    }
  } else {
    html += '<div class="missing">Ignored bin — depth unavailable in clean coordinates</div>';
  }
  tooltip.innerHTML = html;
  tooltip.hidden = false;

  const margin = 12;
  const preferredLeft = screenX + 16;
  const preferredTop = screenY + 16;
  const width = tooltip.offsetWidth;
  const height = tooltip.offsetHeight;
  tooltip.style.left = `${Math.max(
    margin,
    Math.min(preferredLeft, state.width - width - margin),
  )}px`;
  tooltip.style.top = `${Math.max(
    margin,
    Math.min(preferredTop, state.height - height - margin),
  )}px`;
}

function pointerPosition(event) {
  const rect = canvas.getBoundingClientRect();
  return [event.clientX - rect.left, event.clientY - rect.top];
}

function updateHover(event) {
  const [x, y] = pointerPosition(event);
  state.hover = hitTest(x, y);
  canvas.classList.toggle("snapped", Boolean(state.hover));
  updateTooltip(state.hover, x, y);
  requestRender();
}

function setZoom(newScale, screenX, screenY) {
  const oldScale = state.camera.scale;
  const clamped = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, newScale));
  if (clamped === oldScale) return;
  const centerX = state.width / 2;
  const centerY = state.height / 2;
  const worldX = (screenX - centerX - state.camera.offsetX) / oldScale;
  const worldY = (screenY - centerY - state.camera.offsetY) / oldScale;
  state.camera.scale = clamped;
  state.camera.offsetX = screenX - centerX - worldX * clamped;
  state.camera.offsetY = screenY - centerY - worldY * clamped;
  if (clamped === MIN_ZOOM) {
    state.camera.offsetX = 0;
    state.camera.offsetY = 0;
  }
  zoomValue.textContent = `${Math.round(clamped * 100)}%`;
  requestRender();
}

function resetView() {
  state.camera.scale = 1;
  state.camera.offsetX = 0;
  state.camera.offsetY = 0;
  state.hover = null;
  canvas.classList.remove("snapped");
  tooltip.hidden = true;
  zoomValue.textContent = "100%";
  requestRender();
}

canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  const [x, y] = pointerPosition(event);
  const factor = Math.exp(-event.deltaY * 0.0015);
  setZoom(state.camera.scale * factor, x, y);
  state.hover = hitTest(x, y);
  canvas.classList.toggle("snapped", Boolean(state.hover));
  updateTooltip(state.hover, x, y);
}, { passive: false });

canvas.addEventListener("pointerdown", (event) => {
  const [x, y] = pointerPosition(event);
  state.dragging = true;
  state.dragX = x;
  state.dragY = y;
  canvas.classList.add("dragging");
  canvas.setPointerCapture(event.pointerId);
});

canvas.addEventListener("pointermove", (event) => {
  const [x, y] = pointerPosition(event);
  if (state.dragging) {
    state.camera.offsetX += x - state.dragX;
    state.camera.offsetY += y - state.dragY;
    state.dragX = x;
    state.dragY = y;
    state.hover = null;
    canvas.classList.remove("snapped");
    tooltip.hidden = true;
    requestRender();
    return;
  }
  updateHover(event);
});

function endDrag(event) {
  if (!state.dragging) return;
  state.dragging = false;
  canvas.classList.remove("dragging");
  if (canvas.hasPointerCapture(event.pointerId)) {
    canvas.releasePointerCapture(event.pointerId);
  }
  updateHover(event);
}

canvas.addEventListener("pointerup", endDrag);
canvas.addEventListener("pointercancel", endDrag);
canvas.addEventListener("pointerleave", () => {
  if (!state.dragging) {
    state.hover = null;
    canvas.classList.remove("snapped");
    tooltip.hidden = true;
    requestRender();
  }
});
canvas.addEventListener("dblclick", resetView);
resetButton.addEventListener("click", resetView);
lineWidthInput.addEventListener("input", () => {
  state.lineWidthScale = Number(lineWidthInput.value);
  lineWidthValue.value = `${state.lineWidthScale.toFixed(1)}×`;
  requestRender();
});

new ResizeObserver(resizeCanvas).observe(canvasWrap);

fetch("/api/depth", { cache: "no-store" })
  .then((response) => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  })
  .then((data) => {
    state.data = data;
    buildSectors();
    sampleName.textContent = data.sample;
    datasetSummary.textContent = `${data.clean_bins.toLocaleString("en-US")} / `
      + `${data.total_bins.toLocaleString("en-US")} clean bins displayed · `
      + `${data.ignored_bins.toLocaleString("en-US")} ignored · `
      + `${data.reference} · median depth ${formatDepth(data.median_depth)}`;
    loading.hidden = true;
    resizeCanvas();
  })
  .catch((error) => {
    loading.classList.add("error");
    loading.textContent = `Could not load SKYPE depth data: ${error.message}`;
    sampleName.textContent = "SKYPE depth viewer error";
    datasetSummary.textContent = "Check the server terminal for details.";
  });
