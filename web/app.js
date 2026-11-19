// pypedigree in the browser: loads Pyodide and the pypedigree wheel, then
// draws the pedigree typed, pasted, dropped or opened by the user.
// Everything runs locally; the only downloads are Pyodide (from its CDN)
// and the files of this site.

const $ = (id) => document.getElementById(id);
const text = $("text");
const form = $("options");
const status = $("status");

let config = null;
let py = null;          // the pypedigree_web module, once loaded
let result = null;      // the last successful drawing
let baseName = "pedigree";
let timer = null;

function setStatus(message, isError = false) {
  status.textContent = message;
  status.classList.toggle("error", isError);
}

async function loadProgram() {
  config = await (await fetch("config.json")).json();
  $("version").textContent = config.version;
  showExamples(config.examples);

  const cdn = `https://cdn.jsdelivr.net/pyodide/v${config.pyodide}/full/`;
  const { loadPyodide } = await import(cdn + "pyodide.mjs");
  const pyodide = await loadPyodide({ indexURL: cdn });
  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");
  // Without its dependencies: matplotlib, numpy and click are only needed
  // for PNG/PDF output from Python and for the command line.
  const wheel = new URL(config.wheel, location.href).href;
  await micropip.install(wheel, { deps: false });
  const glue = await (await fetch("pypedigree_web.py")).text();
  pyodide.FS.writeFile("/home/pyodide/pypedigree_web.py", glue);
  py = pyodide.pyimport("pypedigree_web");
}

function showExamples(examples) {
  const list = $("examples");
  for (const example of examples) {
    const item = document.createElement("li");
    const button = document.createElement("button");
    button.type = "button";
    const title = document.createElement("strong");
    title.textContent = example.title;
    const description = document.createElement("span");
    description.textContent = example.description;
    button.append(title, description);
    button.addEventListener("click", () => openExample(example));
    item.append(button);
    list.append(item);
  }
}

// An example is named in the page's address by its file name without the
// extension, e.g. ?example=quadrants, so a link opens the page on it.
function exampleName(example) {
  return example.file.split("/").pop().replace(/\.[^.]+$/, "");
}

function setExampleInAddress(name) {
  const url = new URL(location.href);
  if (name) url.searchParams.set("example", name);
  else url.searchParams.delete("example");
  history.replaceState(null, "", url);
}

async function openExample(example) {
  setExampleInAddress(exampleName(example));
  const response = await fetch(example.file);
  text.value = await response.text();
  baseName = example.file.split("/").pop().replace(/\.[^.]+$/, "");
  $("file-name").textContent = `Example: ${example.title}`;
  form.reset();
  for (const [name, value] of Object.entries(example.options || {})) {
    form.elements[name].value = value;
  }
  $("style").value = example.style ? JSON.stringify(example.style, null, 2) : "";
  draw({ newFile: true });
}

function readFile(file) {
  const reader = new FileReader();
  reader.onload = () => {
    setExampleInAddress(null);
    text.value = reader.result;
    baseName = file.name.replace(/\.[^.]+$/, "") || "pedigree";
    $("file-name").textContent = file.name;
    draw({ newFile: true });
  };
  reader.onerror = () => setStatus(`Could not read ${file.name}.`, true);
  reader.readAsText(file);
}

function options() {
  const f = form.elements;
  return {
    family: f.family.value,
    phenotype_layout: f.phenotype_layout.value,
    pgs_scale: f.pgs_scale.value,
    carrier_notation: f.carrier_notation.value,
    min_font_size: f.min_font_size.value,
    show_labels: f.show_labels.checked,
    show_legend: f.show_legend.checked,
    show_na_labels: f.show_na_labels.checked,
    show_pgs_values: f.show_pgs_values.checked,
    style: f.style.value,
  };
}

function draw({ newFile = false } = {}) {
  if (!py) return;  // drawn once the program is loaded
  if (!text.value.trim()) {
    setStatus("Paste, drop or open a pedigree file, or try an example.");
    showDrawing(null);
    return;
  }
  const chosen = options();
  if (newFile) chosen.family = "";
  const out = JSON.parse(py.render(text.value, JSON.stringify(chosen)));
  showWarnings(out.warnings);
  if (!out.ok) {
    setStatus(out.error, true);
    showDrawing(null);
    return;
  }
  showFamilies(out.families, out.family);
  const families = out.families.length;
  const scope = out.family ? `family ${out.family}` :
    families > 1 ? `${families} families` : `family ${out.families[0]}`;
  setStatus(`${out.individuals} individuals, ${scope}.`);
  showDrawing(out);
}

function scheduleDraw() {
  clearTimeout(timer);
  timer = setTimeout(() => draw(), 300);
}

function showFamilies(families, chosen) {
  const select = form.elements.family;
  select.replaceChildren(new Option("All families", ""));
  for (const id of families) select.append(new Option(id, id));
  select.value = chosen;
  select.disabled = families.length < 2;
}

function showWarnings(warnings) {
  const list = $("warnings");
  list.replaceChildren();
  for (const message of warnings) {
    const item = document.createElement("li");
    item.textContent = message;
    list.append(item);
  }
}

function svgSize(svg) {
  const doc = new DOMParser().parseFromString(svg, "image/svg+xml");
  const root = doc.documentElement;
  return { width: parseFloat(root.getAttribute("width")), height: parseFloat(root.getAttribute("height")) };
}

function showDrawing(out) {
  result = out;
  $("viewer").hidden = $("downloads").hidden = $("hint").hidden = $("fit-label").hidden = !out;
  const frame = $("frame");
  if (!out) {
    frame.srcdoc = "";
    return;
  }
  let svg = out.svg;
  if (!$("fit").checked) {
    // Full size: the size of the static drawing, scrolled inside the frame.
    const { width, height } = svgSize(out.static_svg);
    svg = svg.replace('width="100%"', `width="${width}" height="${height}"`);
  }
  const page = "<!doctype html><html><head><meta charset='utf-8'>" +
    "<style>html,body{margin:0;background:white}</style></head><body>" +
    svg + "</body></html>";
  frame.srcdoc = page;
  fitFrame();
}

function fitFrame() {
  if (!result) return;
  const { width, height } = svgSize(result.static_svg);
  const frame = $("frame");
  const wanted = $("fit").checked ? frame.clientWidth * height / width : height + 20;
  frame.style.height = `${Math.max(12 * 16, Math.min(wanted + 4, window.innerHeight * 0.85))}px`;
}

function save(blob, name) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function svgBlob(svg) {
  return new Blob([svg], { type: "image/svg+xml" });
}

async function savePng() {
  const { width, height } = svgSize(result.static_svg);
  const scale = Number($("dpi").value) / 96;  // the SVG's size is in CSS pixels (96 per inch)
  const url = URL.createObjectURL(svgBlob(result.static_svg));
  try {
    const image = new Image();
    image.src = url;
    await image.decode();
    const canvas = document.createElement("canvas");
    canvas.width = Math.round(width * scale);
    canvas.height = Math.round(height * scale);
    const context = canvas.getContext("2d");
    context.fillStyle = "white";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0, canvas.width, canvas.height);
    const png = await new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
    if (!png) throw new Error("the drawing is too large for a PNG at this resolution");
    save(png, `${baseName}.png`);
  } catch (error) {
    setStatus(`Could not make the PNG: ${error.message}. Try a lower resolution, or the SVG.`, true);
  } finally {
    URL.revokeObjectURL(url);
  }
}

function wire() {
  text.addEventListener("input", scheduleDraw);
  form.addEventListener("change", () => draw());
  form.addEventListener("input", (event) => { if (event.target.name === "style") scheduleDraw(); });
  form.addEventListener("submit", (event) => event.preventDefault());
  $("file").addEventListener("change", (event) => {
    if (event.target.files.length) readFile(event.target.files[0]);
    event.target.value = "";
  });
  const drop = $("drop");
  drop.addEventListener("dragover", (event) => { event.preventDefault(); drop.classList.add("over"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("over"));
  drop.addEventListener("drop", (event) => {
    event.preventDefault();
    drop.classList.remove("over");
    if (event.dataTransfer.files.length) readFile(event.dataTransfer.files[0]);
  });
  $("dl-svg").addEventListener("click", () => save(svgBlob(result.svg), `${baseName}.svg`));
  $("dl-static").addEventListener("click", () => save(svgBlob(result.static_svg), `${baseName}_print.svg`));
  $("dl-png").addEventListener("click", savePng);
  $("fit").addEventListener("change", () => showDrawing(result));
  window.addEventListener("resize", fitFrame);
}

wire();
loadProgram()
  .then(() => {
    setStatus("Ready. Paste, drop or open a pedigree file, or try an example.");
    const wanted = new URLSearchParams(location.search).get("example");
    const example = wanted && config.examples.find((e) => exampleName(e) === wanted);
    if (example) openExample(example);
    else if (wanted) setStatus(`No example named "${wanted}". Choose one below.`, true);
    else if (text.value.trim()) draw({ newFile: true });
  })
  .catch((error) => {
    console.error(error);
    setStatus("The program could not be loaded. Check your internet connection " +
      "(it is downloaded once from cdn.jsdelivr.net) and reload the page. " +
      `Details: ${error.message}`, true);
  });
