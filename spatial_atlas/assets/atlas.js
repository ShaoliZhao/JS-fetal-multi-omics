const geneInput = document.querySelector("#geneInput");
const geneList = document.querySelector("#geneList");
const plotButton = document.querySelector("#plotButton");
const pointSize = document.querySelector("#pointSize");
const capSlider = document.querySelector("#capSlider");
const currentGene = document.querySelector("#currentGene");
const dataMode = document.querySelector("#dataMode");
const cellCount = document.querySelector("#cellCount");
const ctrlCanvas = document.querySelector("#ctrlCanvas");
const jsCanvas = document.querySelector("#jsCanvas");
const ctrlMax = document.querySelector("#ctrlMax");
const jsMax = document.querySelector("#jsMax");

const fallbackGenes = ["ATOH1", "PCP4", "WNT7B", "SOX2", "MKI67", "FOXJ1"];
let manifest = { genes: fallbackGenes, demo: true };
let points = [];
let expression = [];

function seededNoise(index, seed) {
  const x = Math.sin(index * 12.9898 + seed * 78.233) * 43758.5453;
  return x - Math.floor(x);
}

function makeDemoPoints() {
  const celltypes = ["RL", "iGL", "PKC", "VZP", "BG/AST", "Endo"];
  const out = [];
  for (const sample of ["ctrl", "JS"]) {
    for (let i = 0; i < 900; i += 1) {
      const ring = i % 6;
      const angle = i * 0.28 + ring * 0.4;
      const radius = 0.18 + ring * 0.065 + seededNoise(i, sample === "ctrl" ? 1 : 2) * 0.035;
      const bend = sample === "ctrl" ? 0.03 : -0.04;
      const x = 0.5 + Math.cos(angle) * radius + bend * Math.sin(i * 0.013);
      const y = 0.5 + Math.sin(angle) * radius * 0.78 + (ring - 2.5) * 0.018;
      out.push({
        id: `${sample}_${i}`,
        sample,
        x: Math.max(0.02, Math.min(0.98, x)),
        y: Math.max(0.02, Math.min(0.98, y)),
        celltype: celltypes[ring]
      });
    }
  }
  return out;
}

function demoExpression(gene) {
  const seed = gene.split("").reduce((total, char) => total + char.charCodeAt(0), 0);
  return points.map((point, i) => {
    const ridge = Math.exp(-Math.pow(point.x - 0.46, 2) / 0.02) * Math.exp(-Math.pow(point.y - 0.42, 2) / 0.06);
    const outer = Math.exp(-Math.pow(point.y - 0.72, 2) / 0.025);
    const jsShift = point.sample === "JS" ? 0.68 : 1;
    const boost = gene.startsWith("WNT") && point.sample === "JS" ? 1.6 : 1;
    const pkc = gene === "PCP4" ? Math.exp(-Math.pow(point.y - 0.35, 2) / 0.018) * 1.8 : 0;
    const cycling = gene === "MKI67" ? Math.exp(-Math.pow(point.x - 0.68, 2) / 0.018) : 0;
    const noise = seededNoise(i, seed) * 0.2;
    return Math.max(0, (ridge * jsShift + outer * 0.55 + pkc + cycling + noise) * boost);
  });
}

function colorForValue(value, cap) {
  if (value <= 0) return "#dfe7e3";
  const t = Math.max(0, Math.min(1, value / cap));
  if (t < 0.34) return blend("#dfe7e3", "#6bb8aa", t / 0.34);
  if (t < 0.68) return blend("#6bb8aa", "#2e7a9b", (t - 0.34) / 0.34);
  return blend("#2e7a9b", "#f0bd4f", (t - 0.68) / 0.32);
}

function blend(a, b, t) {
  const ca = hexToRgb(a);
  const cb = hexToRgb(b);
  const r = Math.round(ca.r + (cb.r - ca.r) * t);
  const g = Math.round(ca.g + (cb.g - ca.g) * t);
  const bl = Math.round(ca.b + (cb.b - ca.b) * t);
  return `rgb(${r}, ${g}, ${bl})`;
}

function hexToRgb(hex) {
  const raw = hex.replace("#", "");
  return {
    r: parseInt(raw.slice(0, 2), 16),
    g: parseInt(raw.slice(2, 4), 16),
    b: parseInt(raw.slice(4, 6), 16)
  };
}

function percentile(values, p) {
  const sorted = values.filter((v) => Number.isFinite(v)).sort((a, b) => a - b);
  if (sorted.length === 0) return 1;
  const idx = Math.min(sorted.length - 1, Math.floor((p / 100) * sorted.length));
  return Math.max(sorted[idx], 0.001);
}

function drawSlice(canvas, sample, cap) {
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#f9faf6";
  ctx.fillRect(0, 0, width, height);

  const samplePoints = points
    .map((point, i) => ({ ...point, value: expression[i] || 0 }))
    .filter((point) => point.sample === sample);
  const radius = Number(pointSize.value);

  for (const point of samplePoints) {
    ctx.fillStyle = colorForValue(point.value, cap);
    ctx.globalAlpha = point.value > 0 ? 0.86 : 0.32;
    ctx.beginPath();
    ctx.arc(point.x * width, point.y * height, radius, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.globalAlpha = 1;
}

async function loadManifest() {
  try {
    const response = await fetch("data/manifest.json");
    if (response.ok) {
      manifest = await response.json();
    }
  } catch {
    manifest = { genes: fallbackGenes, demo: true };
  }
  geneList.innerHTML = "";
  for (const gene of manifest.genes || fallbackGenes) {
    const option = document.createElement("option");
    option.value = gene;
    geneList.appendChild(option);
  }
}

async function loadPoints() {
  try {
    const response = await fetch("data/points.json");
    if (response.ok) {
      const data = await response.json();
      points = data.points || [];
      dataMode.textContent = "exported";
      return;
    }
  } catch {
    points = [];
  }
  points = makeDemoPoints();
  dataMode.textContent = "demo";
}

async function loadGene(gene) {
  try {
    const response = await fetch(`data/genes/${encodeURIComponent(gene)}.json`);
    if (response.ok) {
      const data = await response.json();
      expression = data.values || [];
      dataMode.textContent = "exported";
      return;
    }
  } catch {
    expression = [];
  }
  expression = demoExpression(gene);
  dataMode.textContent = "demo";
}

async function updatePlot() {
  const gene = geneInput.value.trim().toUpperCase() || "ATOH1";
  geneInput.value = gene;
  currentGene.textContent = gene;
  await loadGene(gene);
  const cap = percentile(expression, Number(capSlider.value));
  drawSlice(ctrlCanvas, "ctrl", cap);
  drawSlice(jsCanvas, "JS", cap);

  const ctrlValues = points.map((point, i) => point.sample === "ctrl" ? expression[i] || 0 : null).filter((v) => v !== null);
  const jsValues = points.map((point, i) => point.sample === "JS" ? expression[i] || 0 : null).filter((v) => v !== null);
  ctrlMax.textContent = `max ${Math.max(...ctrlValues).toFixed(2)}`;
  jsMax.textContent = `max ${Math.max(...jsValues).toFixed(2)}`;
  cellCount.textContent = String(points.length);
}

plotButton.addEventListener("click", updatePlot);
geneInput.addEventListener("keydown", (event) => {
  if (event.key === "Enter") updatePlot();
});
pointSize.addEventListener("input", updatePlot);
capSlider.addEventListener("input", updatePlot);

loadManifest()
  .then(loadPoints)
  .then(updatePlot);
