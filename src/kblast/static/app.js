"use strict";

const $ = (id) => document.getElementById(id);
const token = document.querySelector('meta[name="kblast-token"]').content;
const state = { run: null, runId: null, stamp: null, tab: "hits", selected: new Set(),
  seen: new Set(), current: null, cache: new Map(), hitPage: 0, occurrencePage: 0, desktop: false };
const labels = { success: "Success", partial_success: "Finished with exceptions", running: "Running",
  pending: "Pending", no_mapping: "No mapping found", unresolved_mapping: "Unresolved mapping",
  retrieval_failure: "Retrieval failure", parse_failure: "Response parse failure", failure: "Run stopped", cancelled: "Cancelled",
  interrupted: "Interrupted", no_hits: "No hits", integrity_failure: "Target CDS mismatch" };
const esc = (value) => String(value ?? "").replace(/[&<>"']/g, (c) =>
  ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const number = (value) => value == null ? "—" : Number(value).toLocaleString("en-US");
const pct = (value) => Number(value).toFixed(1) + "%";
const badge = (status) => `<span class="status-badge ${esc(status)}">${esc(labels[status] || status)}</span>`;
const errorNote = (error) => error ? `<span class="row-error">${esc(error)}</span>` : "";
const mappingLabel = (o) => o.mapping_method === "ELink + GenBank CDS" ? "ELink / GenBank CDS" : `IPG row ${o.ipg_row}`;
const PAGE_SIZE = 100;
let toastTimer;
let ridTarget = null;
let importingRid = false;
// Two observed standard browser cycles (~138 s each), rounded to five minutes.
const WAITING_RECOVERY_SECONDS = 300;

async function api(path, data) {
  const options = { headers: { "X-Kblast-Token": token } };
  if (data !== undefined) {
    options.method = "POST";
    options.headers["Content-Type"] = "application/json";
    options.body = JSON.stringify(data);
  }
  const response = await fetch("/api/" + path, options);
  if (!response.ok) {
    let message = `Local request failed (${response.status}).`;
    try { message = (await response.json()).error || message; } catch (_) { /* Non-JSON error */ }
    throw new Error(message);
  }
  return response;
}

function toast(message) {
  $("toast").textContent = message;
  $("toast").hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { $("toast").hidden = true; }, 9000);
}

function empty(title, message) {
  return `<div class="empty-state"><svg viewBox="0 0 190 65" aria-hidden="true">
    <path class="empty-line" d="M8 31H182"/><path class="empty-gene" d="M14 21h25l9 10-9 10H14z"/>
    <path class="empty-target" d="M61 21h42l10 10-10 10H61z"/><path class="empty-gene" d="M163 21h-22l-10 10 10 10h22z"/>
    </svg><h3>${esc(title)}</h3><p>${esc(message)}</p></div>`;
}

function sequenceCount() {
  const lines = $("sequence").value.trim().split(/\r?\n/);
  const text = lines.filter((line) => !line.trimStart().startsWith(">")).join("");
  $("residue-count").textContent = number(text.replace(/\s/g, "").replace(/\*$/, "").length) + " aa";
}

function switchTab(tab) {
  state.tab = tab;
  for (const name of ["hits", "occurrences", "contexts"]) {
    $("tab-" + name).setAttribute("aria-selected", String(name === tab));
    $("tab-" + name).tabIndex = name === tab ? 0 : -1;
    $("panel-" + name).hidden = name !== tab;
  }
  if (tab === "contexts") renderContexts();
}

function updateControls() {
  const run = state.run;
  $("run-button").disabled = !!run?.active;
  $("run-button").innerHTML = run?.active ? 'Search in progress <span aria-hidden="true">…</span>' : 'Run search <span aria-hidden="true">→</span>';
  $("cancel-button").hidden = !run?.active;
  $("cancel-button").disabled = !!run?.cancelling;
  $("cancel-button").textContent = run?.cancelling ? "Cancelling…" : "Cancel run";
  $("resume-button").hidden = !run || run.active || ["success", "no_hits"].includes(run.status);
  $("export-all").disabled = !run;
  $("export-selected").disabled = !run || state.selected.size === 0;
  $("use-rid").disabled = !!run?.active || importingRid;
}

function renderRun() {
  const run = state.run;
  updateControls();
  if (!run) {
    $("hits-content").innerHTML = empty("Start with a protein.", "Your ranked BLAST hits and their genomic neighborhoods will appear here.");
    $("occurrences-content").innerHTML = empty("One hit can lead to many places.", "Nucleotide-backed IPG occurrences will appear after BLAST completes.");
    $("contexts-content").innerHTML = empty("Inspect the neighborhood.", "Annotated gene arrows appear as genomic contexts are retrieved.");
    return;
  }
  for (const [id, value] of [["hit", run.counts.hits], ["occurrence", run.counts.occurrences], ["context", run.counts.contexts]]) {
    $(id + "-count").textContent = number(value);
    $("tab-" + id + "-count").textContent = number(value);
  }
  const databaseLabel = run.settings.database === "nr_cluster_seq" ? "ClusteredNR" : run.settings.database;
  $("result-subtitle").textContent = `${number(run.query.sequence.length)} aa · top ${number(run.settings.top_n)} ${databaseLabel} hits · ±${number(run.settings.flank)} bp`;
  $("run-status").hidden = false;
  $("run-status-label").className = "status-badge " + run.status;
  $("run-status-label").textContent = run.cancelling ? "Cancelling" : labels[run.status] || run.status;
  $("run-id").textContent = run.rid ? "RID " + run.rid : "Run " + run.id.slice(0, 8);
  if (run.rid) {
    $("run-id").href = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=HTML&PROGRAM=blastp&PAGE=Proteins&RID=" + encodeURIComponent(run.rid);
    $("run-id").title = "Inspect this RID at NCBI Protein BLAST";
  } else { $("run-id").removeAttribute("href"); $("run-id").removeAttribute("title"); }
  $("run-message").textContent = run.message;
  renderBlastLifecycle();
  $("run-error").hidden = !run.error;
  $("run-error").textContent = run.error;
  $("progress-line").classList.toggle("done", !run.active);
  $("progress-line").hidden = !run.active && run.status !== "success" && run.status !== "partial_success";
  $("recovery").hidden = !!run.blast_metadata.program || !!run.hits.length;
  $("resume-original-rid").hidden = !run.original_rid || run.original_rid === run.rid;
  const external = ["external", "browser"].includes(run.rid_source);
  $("rid-provenance").hidden = !external && !run.original_rid;
  const metadata = run.blast_metadata;
  const available = metadata.returned_hit_count;
  const actual = Object.entries(metadata.parameters || {}).map(([k,v]) => `${k}=${v}`).join(" · ");
  $("rid-provenance").textContent = `${external ? "Externally supplied BLAST result ("+run.rid_source+")" : "Kblast URL API submission"}${run.original_rid && run.original_rid !== run.rid ? " · original RID "+run.original_rid : ""}${run.rid_adopted_at ? " · adopted "+new Date(run.rid_adopted_at).toLocaleString() : ""}${metadata.program ? " · returned "+metadata.program+" / "+metadata.database+" / "+metadata.version : ""}${available !== undefined ? " · "+available+" hits returned, consuming "+run.hits.length+" of requested "+run.settings.top_n : ""}${actual ? " · "+actual : ""}${external && run.rid_validation.status !== "validated" ? " · validation "+run.rid_validation.status : ""}${metadata.query_validation?.residual_risk ? " · Query identity is only partially verified; confirm the original protein." : ""}`;
  $("input-warnings").hidden = !run.query.warnings.length;
  $("input-warnings").textContent = run.query.warnings.join(" ");
  for (const o of run.occurrences) {
    if (o.status === "success" && !state.seen.has(o.id)) {
      state.seen.add(o.id);
      state.selected.add(o.id);
    }
  }
  renderHits();
  renderOccurrences();
  if (state.tab === "contexts") renderContexts();
  updateControls();
}

function renderBlastLifecycle() {
  const run = state.run;
  if (!run) return;
  const phases = {submitting:'Submitted', waiting:'Waiting at NCBI', ready:'Results ready',
    retrieving:'Retrieving', parsing:'Parsing', complete:'Parsed', cancelled:'Cancelled', failed:'Stopped'};
  const started = run.rid_adopted_at || run.submitted_at;
  const elapsed = started ? Math.max(0, Math.floor(((run.active ? Date.now() : Date.parse(run.updated_at)) - Date.parse(started)) / 1000)) : 0;
  const duration = `${Math.floor(elapsed / 60)}m ${elapsed % 60}s`;
  const poll = run.last_blast_get_at ? ` · last poll ${new Date(run.last_blast_get_at).toLocaleTimeString()}` : '';
  $('blast-lifecycle').textContent = `${phases[run.blast_phase] || 'Saved BLAST job'}${run.rid ? ' · RID received' : ''}${started ? ' · '+duration+' elapsed locally' : ''}${run.ncbi_status ? ' · NCBI '+run.ncbi_status : ''}${poll}`;
  $('waiting-recovery').hidden = !(run.rid_source === 'kblast_url_api' && run.ncbi_status === 'WAITING' && !run.blast_metadata.program && elapsed >= WAITING_RECOVERY_SECONDS);
}
setInterval(renderBlastLifecycle, 1000);

function pager(kind, total, page) {
  if (total <= PAGE_SIZE) return "";
  return `<div class="table-note">${number(page * PAGE_SIZE + 1)}–${number(Math.min(total, (page + 1) * PAGE_SIZE))} of ${number(total)} rows
    <button class="quiet" data-page="${kind}" data-delta="-1" ${page === 0 ? "disabled" : ""}>← Previous</button>
    <button class="quiet" data-page="${kind}" data-delta="1" ${(page + 1) * PAGE_SIZE >= total ? "disabled" : ""}>Next →</button></div>`;
}

function renderHits() {
  const hits = state.run.hits;
  if (!hits.length) {
    $("hits-content").innerHTML = empty(state.run.active ? "NCBI is processing the search." : "No protein hits returned.",
      state.run.active ? "The search runs in the background. Its RID and settings are saved locally." : "Check the run status above for details.");
    return;
  }
  const rows = hits.slice(state.hitPage * PAGE_SIZE, (state.hitPage + 1) * PAGE_SIZE).map((h) => `<tr>
    <td class="mono">${h.rank}</td><td class="mono accession">${esc(h.accession)}</td>
    <td class="description">${esc(h.title)}<small>${esc(h.organism || "Organism unavailable")}</small></td>
    <td>${pct(h.percent_identity)}</td><td>${pct(h.query_coverage)}</td><td class="mono">${Number(h.evalue).toExponential(2)}</td><td>${Number(h.bit_score).toFixed(1)}</td>
    <td>${badge(h.mapping_status)}${errorNote(h.mapping_error)}</td><td>${badge(h.sequence_status)}${errorNote(h.sequence_error)}</td></tr>`).join("");
  $("hits-content").innerHTML = `<div class="table-scroll"><table class="data-table"><thead><tr><th>Rank</th><th>Accession</th><th>Description / organism</th><th>Max identity</th><th>Query coverage (union)</th><th>Best E-value</th><th>Max bit score</th><th>IPG mapping</th><th>Protein FASTA</th></tr></thead><tbody>${rows}</tbody></table></div>${pager("hit", hits.length, state.hitPage)}`;
}

function renderOccurrences() {
  const rows = state.run.occurrences;
  if (!rows.length) {
    $("occurrences-content").innerHTML = empty(state.run.active ? "Waiting for occurrence expansion." : "No genomic occurrences.", state.run.status === "no_hits" ? "BLAST returned no subjects; occurrence expansion was not needed." : "Occurrences require a nucleotide accession and usable CDS location from IPG or NCBI nucleotide links. Hit mapping status is shown in Preview 1.");
    return;
  }
  const content = rows.slice(state.occurrencePage * PAGE_SIZE, (state.occurrencePage + 1) * PAGE_SIZE).map((o) => {
    const clip = (o.clipped_left ? '<span class="clip-badge">Left clipped</span>' : "") + (o.clipped_right ? '<span class="clip-badge">Right clipped</span>' : "");
    return `<tr ${o.status === "success" ? `class="clickable-row" data-context="${esc(o.id)}" tabindex="0" aria-label="Preview ${esc(o.nucleotide_accession)} context"` : ""}>
      <td class="mono">${o.hit_rank}</td><td class="mono accession">${esc(o.hit_accession)}</td><td class="mono accession">${esc(o.nucleotide_accession) || "—"}</td><td>${badge(o.status)}${errorNote(o.error)}</td>
      <td class="mono">${number(o.cds_start)}</td><td class="mono">${number(o.cds_end)}</td><td>${esc(o.strand) || "—"}</td>
      <td class="description">${esc(o.organism) || "—"}<small>${esc(o.strain)}</small></td><td class="mono">${esc(o.assembly) || "—"}</td>
      <td>${esc(o.source) || "—"}<small class="subtle"> · ${esc(mappingLabel(o))}</small></td><td class="mono">${number(o.requested_start)} … ${number(o.requested_end)}</td>
      <td class="mono">${number(o.context_start)} … ${number(o.context_end)}</td><td>${clip || (o.clipped_left === false ? "None" : "—")}</td>
      <td class="mono">${esc(o.protein_accession) || "—"}</td></tr>`;
  }).join("");
  $("occurrences-content").innerHTML = `<div class="table-scroll"><table class="data-table"><thead><tr><th>Rank</th><th>BLAST subject</th><th>Nucleotide</th><th>Status</th><th>CDS low</th><th>CDS high</th><th>Strand</th><th>Organism / strain</th><th>Assembly</th><th>Source / row</th><th>Requested context</th><th>Context bounds</th><th>Clipping</th><th>Protein accession</th></tr></thead><tbody>${content}</tbody></table></div>${pager("occurrence", rows.length, state.occurrencePage)}<div class="table-note">${number(state.run.counts.hits)} subject hits → ${number(rows.length)} genomic occurrences · ${number(state.run.counts.unresolved)} unresolved · ${number(state.run.counts.failed)} retrieval failures · ${number(state.run.counts.integrity_failures)} target CDS mismatches. Repeated rows are kept.</div>`;
}

function renderContexts() {
  const successful = state.run?.occurrences.filter((o) => o.status === "success") || [];
  $("selection-label").textContent = successful.length ? `${number(state.selected.size)} of ${number(successful.length)} selected · native genomic orientation` : "All resolved contexts are selected by default.";
  updateControls();
  if (!successful.length) {
    $("contexts-content").innerHTML = empty(state.run?.active ? "Contexts appear as retrieval completes." : "No successful genomic contexts.", "Successful contexts require a CDS matching the mapped coordinates and strand. Inspect hit and occurrence statuses for details.");
    return;
  }
  if (!successful.some((o) => o.id === state.current)) state.current = successful[0].id;
  const list = successful.map((o) => `<div class="context-item ${state.current === o.id ? "selected-preview" : ""}">
    <input type="checkbox" data-select="${esc(o.id)}" aria-label="Export ${esc(o.id)} ${esc(o.nucleotide_accession)}" ${state.selected.has(o.id) ? "checked" : ""}>
    <button data-preview="${esc(o.id)}"><span class="mono accession">${esc(o.nucleotide_accession)}</span><span class="item-organism">${esc(o.organism)}</span><span class="item-coordinates">${number(o.context_start)}–${number(o.context_end)} · hit ${o.hit_rank}, ${esc(mappingLabel(o))}</span></button></div>`).join("");
  $("contexts-content").innerHTML = `<div class="context-layout"><div class="context-list"><div class="selection-actions"><button class="text-button" data-selection="all">Select all</button><button class="text-button secondary-text" data-selection="none">Clear</button></div>${list}</div><div id="context-view" class="context-view"><div class="loading">Loading annotation…</div></div></div>`;
  loadContext(state.current);
}

async function loadContext(id) {
  const runId = state.runId;
  const occurrence = state.run.occurrences.find((o) => o.id === id);
  const key = runId + ":" + id + ":" + occurrence.context_sha256;
  try {
    let data = state.cache.get(key);
    if (!data) {
      data = await (await api(`runs/${runId}/contexts/${id}`)).json();
      state.cache.set(key, data);
    }
    if (state.current !== id || state.runId !== runId || !$("context-view")) return;
    const clip = (occurrence.clipped_left ? '<span class="clip-badge">Left clipped</span>' : "") + (occurrence.clipped_right ? '<span class="clip-badge">Right clipped</span>' : "");
    $("context-view").innerHTML = `<div class="context-heading"><div><h3>${esc(data.organism || occurrence.organism)}</h3><p>Hit ${occurrence.hit_rank} · ${esc(occurrence.hit_accession)} · ${esc(mappingLabel(occurrence))}</p></div><span class="status-badge success">Retrieved</span></div>
      <div class="context-region"><span class="mono">${esc(occurrence.record_accession)}:${number(data.start)}–${number(data.end)}</span><span class="context-native">${number(data.length)} bp</span>${clip}</div>
      <div class="legend"><span><i class="target"></i> Mapped CDS</span><span><i></i> Other CDS</span><span><i class="rna"></i> RNA feature</span></div>
      <div id="arrow-container" class="arrow-container"></div><div id="feature-detail" class="feature-detail"><p>Click a gene arrow to inspect its annotation.</p></div>
      <p class="context-footnote">Native coordinates increase left → right. Arrow direction follows the annotated strand. Joined features retain their separate parts.</p>`;
    if (data.omitted_features) $("context-view").insertAdjacentHTML("beforeend", `<p class="context-footnote">${number(data.omitted_features)} features with unavailable or remote coordinates are retained in GenBank and omitted from the arrow view.</p>`);
    drawArrows(data);
  } catch (error) {
    if (state.current === id && $("context-view")) $("context-view").innerHTML = `<div class="notice error">${esc(error.message)}</div>`;
  }
}

function drawArrows(data) {
  const container = $("arrow-container");
  const width = Math.max(620, Math.floor(container.clientWidth));
  const left = 36, right = width - 22, plotWidth = right - left;
  const scale = (coordinate) => left + (coordinate - data.start) / data.length * plotWidth;
  const occupied = [];
  const placements = [];
  const sorted = [...data.features].sort((a, b) => a.start - b.start || a.end - b.end);
  for (const feature of sorted) {
    const start = scale(feature.start), end = scale(feature.end + 1);
    const label = String(feature.label).slice(0, 12);
    const labelWidth = label.length * 6;
    const min = Math.min(start, (start + end) / 2 - labelWidth / 2);
    const max = Math.max(end, (start + end) / 2 + labelWidth / 2) + 10;
    let lane = occupied.findIndex((previous) => previous < min);
    if (lane < 0) { lane = occupied.length; occupied.push(0); }
    occupied[lane] = max;
    placements.push({ feature, lane, start, end, label });
  }
  const laneCount = Math.max(2, occupied.length);
  const axisY = 30 + laneCount * 45;
  let drawing = "";
  for (let lane = 0; lane < laneCount; lane++) drawing += `<line x1="${left}" x2="${right}" y1="${31 + lane * 45}" y2="${31 + lane * 45}" stroke="#e0e7e1" stroke-width="1"/>`;
  for (const { feature, lane, start, end, label } of placements) {
    const y = 31 + lane * 45;
    const color = feature.target ? "#7960aa" : feature.type === "CDS" || feature.type === "gene" ? "#779ba9" : "#b08438";
    let paths = "";
    if (feature.parts.length > 1) paths += `<line x1="${start}" x2="${end}" y1="${y}" y2="${y}" stroke="${color}" stroke-dasharray="3 3"/>`;
    for (const part of feature.parts) {
      const x1 = Math.max(left, scale(part.start));
      const x2 = Math.min(right, Math.max(x1 + 2, scale(part.end + 1)));
      const tip = Math.min(7, (x2 - x1) / 2), top = y - 7, bottom = y + 7;
      const points = part.strand === -1 ? `${x2},${top} ${x1 + tip},${top} ${x1},${y} ${x1 + tip},${bottom} ${x2},${bottom}`
        : part.strand === 1 ? `${x1},${top} ${x2 - tip},${top} ${x2},${y} ${x2 - tip},${bottom} ${x1},${bottom}`
        : `${x1},${top} ${x2},${top} ${x2},${bottom} ${x1},${bottom}`;
      paths += `<polygon points="${points}" fill="${color}" stroke="${color}" stroke-width="1" ${feature.partial ? 'stroke-dasharray="2 2"' : ""}/>`;
    }
    const title = `${feature.label} · ${feature.product || feature.type} · ${number(feature.start)}–${number(feature.end)} (${feature.strand === -1 ? "−" : feature.strand === 1 ? "+" : "unknown strand"})`;
    drawing += `<g data-feature="${feature.index}" tabindex="0" role="button" aria-label="${esc(title)}"><title>${esc(title)}</title>${paths}<text x="${(start + end) / 2}" y="${y + 23}" font-size="9" fill="${feature.target ? "#7960aa" : "#6a7f73"}" text-anchor="middle">${esc(label)}</text></g>`;
  }
  drawing += `<line x1="${left}" x2="${right}" y1="${axisY}" y2="${axisY}" stroke="#bdccc1"/>`;
  for (let i = 0; i <= 4; i++) {
    const x = left + plotWidth * i / 4;
    const coordinate = Math.round(data.start + (data.end - data.start) * i / 4);
    drawing += `<line x1="${x}" x2="${x}" y1="${axisY}" y2="${axisY + 5}" stroke="#bdccc1"/><text x="${x}" y="${axisY + 21}" font-size="9" fill="#85968a" text-anchor="${i === 0 ? "start" : i === 4 ? "end" : "middle"}">${number(coordinate)}</text>`;
  }
  if (!data.features.length) drawing += `<text x="${width / 2}" y="60" text-anchor="middle" fill="#6d7e79" font-size="12">No drawable CDS or RNA features in this annotation.</text>`;
  container.innerHTML = `<svg style="width:${width}px" viewBox="0 0 ${width} ${axisY + 35}" aria-label="Annotated genomic neighborhood" role="img">${drawing}</svg>`;
  const inspect = (event) => {
    const arrow = event.target.closest("[data-feature]");
    if (!arrow) return;
    const f = data.features.find((feature) => feature.index === Number(arrow.dataset.feature));
    $("feature-detail").innerHTML = `<strong>${esc(f.gene || f.locus_tag || f.type)}</strong>${f.target ? ' <span class="status-badge">Mapped CDS</span>' : ""}<p class="detail-product">${esc(f.product || "No product annotation")}</p><p class="mono">${number(f.start)}–${number(f.end)} · strand ${f.strand === -1 ? "−" : f.strand === 1 ? "+" : "unknown"} · ${esc(f.type)}${f.partial ? " · partial feature" : ""}</p><p>${esc([f.locus_tag, f.protein_id].filter(Boolean).join(" · "))}</p>`;
  };
  container.addEventListener("click", inspect);
  container.addEventListener("keydown", (event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); inspect(event); } });
}

async function refresh() {
  if (!state.runId) return;
  const id = state.runId;
  try {
    const run = await (await api("runs/" + id)).json();
    if (state.runId !== id) return;
    state.run = run;
    const stamp = run.updated_at + ":" + run.active + ":" + run.cancelling;
    if (stamp !== state.stamp) { state.stamp = stamp; renderRun(); }
    else updateControls();
  } catch (error) { toast(error.message); }
}

async function loadRun(id) {
  state.run = null;
  state.runId = id;
  state.stamp = null;
  state.selected.clear(); state.seen.clear(); state.cache.clear();
  state.current = null; state.hitPage = 0; state.occurrencePage = 0;
  await refresh();
  if (state.run) {
    $("sequence").value = ">" + state.run.query.description + "\n" + state.run.query.sequence;
    $("top-n").value = state.run.settings.top_n;
    $("flank").value = state.run.settings.flank;
    sequenceCount();
  }
}

async function exportResults(selected) {
  const id = state.runId;
  if (!id) return;
  try {
    if (state.desktop) {
      const result = await (await api(`runs/${id}/export-desktop`, { selected })).json();
      if (result.error) throw new Error(result.error);
      if (result.path) toast("Exported to " + result.path);
      return;
    }
    const response = await api(`runs/${id}/export`, { selected });
    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = `Kblast_${state.run.created_at.slice(0, 10)}_${id.slice(0, 8)}.zip`;
    document.body.appendChild(anchor); anchor.click(); anchor.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
    toast("Export bundle prepared. Check your downloads.");
  } catch (error) { toast(error.message); }
}

$("run-button").addEventListener("click", async () => {
  $("input-error").hidden = true;
  if (!$("top-n").checkValidity() || !$("flank").checkValidity() || !$("email").checkValidity()) {
    $("input-error").textContent = "Check the hit count, flank size, and contact email.";
    $("input-error").hidden = false; return;
  }
  $("run-button").disabled = true;
  try {
    const result = await (await api("runs", { sequence: $("sequence").value,
      top_n: $("top-n").valueAsNumber, flank: $("flank").valueAsNumber, email: $("email").value })).json();
    switchTab("hits"); await loadRun(result.id);
  } catch (error) { $("input-error").textContent = error.message; $("input-error").hidden = false; updateControls(); }
});
$("cancel-button").addEventListener("click", async () => {
  try { await api(`runs/${state.runId}/cancel`, {}); await refresh(); } catch (error) { toast(error.message); }
});
$("resume-button").addEventListener("click", async () => {
  try { await api(`runs/${state.runId}/resume`, { email: $("email").value }); await refresh(); }
  catch (error) { toast(error.message); }
});

function openRidDialog(target) {
  ridTarget = target;
  $("rid-dialog-title").textContent = target ? "Switch to existing NCBI RID" : "Use existing NCBI RID";
  $("rid-dialog-scope").textContent = target ? `Use the saved ${state.run.query.sequence.length} aa query, top ${state.run.settings.top_n}, ±${number(state.run.settings.flank)} bp.` : "Create a fresh run from the protein, Top N and flank entered on the left.";
  $("existing-rid").value = "";
  $("rid-error").hidden = true;
  $("rid-dialog").showModal();
  $("existing-rid").focus();
}
$("use-rid").addEventListener("click", () => openRidDialog(null));
$("switch-rid").addEventListener("click", () => openRidDialog(state.runId));
$("resume-original-rid").addEventListener("click", () => {
  openRidDialog(state.runId);
  $("existing-rid").value = state.run.original_rid;
  $("rid-source").value = "external";
});
$("close-rid").addEventListener("click", () => { if (!importingRid) $("rid-dialog").close(); });
$("rid-dialog").addEventListener("cancel", (event) => { if (importingRid) event.preventDefault(); });
$("rid-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  if (importingRid) return;
  const rid = $("existing-rid").value.trim();
  if (!/^[A-Z0-9][A-Z0-9-]{5,63}$/.test(rid)) return;
  const id = ridTarget;
  importingRid = true;
  $("adopt-rid").disabled = true;
  $("close-rid").disabled = true;
  $("rid-error").hidden = true;
  try {
    const data = { rid, rid_source: $("rid-source").value, email: $("email").value };
    if (id) {
      let run = await (await api(`runs/${id}`)).json();
      if (run.active) {
        $("adopt-rid").textContent = "Stopping local work…";
        await api(`runs/${id}/cancel`, {});
        while (run.active) {
          await new Promise(resolve => setTimeout(resolve, 1000));
          run = await (await api(`runs/${id}`)).json();
        }
      }
      await api(`runs/${id}/resume`, data);
      await loadRun(id);
    } else {
      if (!$("top-n").checkValidity() || !$("flank").checkValidity() || !$("email").checkValidity()) throw new Error("Check the hit count, flank and contact email.");
      const result = await (await api("runs", {...data, sequence: $("sequence").value, top_n: $("top-n").valueAsNumber, flank: $("flank").valueAsNumber})).json();
      await loadRun(result.id);
    }
    switchTab("hits");
    $("rid-dialog").close();
    toast("Existing RID saved. Retrieving and validating its NCBI result.");
  } catch (error) { $("rid-error").textContent = error.message; $("rid-error").hidden = false; }
  finally { importingRid = false; $("adopt-rid").disabled = false; $("close-rid").disabled = false; $("adopt-rid").textContent = "Retrieve existing result"; updateControls(); }
});
async function openNcbi() {
  try {
    if (state.desktop) await api("open-protein-blast", {});
    else window.open("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins", "_blank", "noopener,noreferrer");
  } catch (error) { toast(error.message); }
}
$("open-ncbi").addEventListener("click", openNcbi);
$("waiting-open-ncbi").addEventListener("click", openNcbi);
$("copy-query").addEventListener("click", async () => {
  try { await navigator.clipboard.writeText($("sequence").value); toast("Query copied. Paste it into NCBI and select true nr and your intended algorithm."); }
  catch (_) { $("sequence").focus(); $("sequence").select(); toast("Query selected. Press ⌘C to copy it."); }
});
$("sequence").addEventListener("input", () => { sequenceCount(); $("input-error").hidden = true; });
$("import-button").addEventListener("click", () => $("file-input").click());
$("file-input").addEventListener("change", async (event) => {
  const file = event.target.files[0];
  if (!file) return;
  if (file.size > 2 * 1024 * 1024) { toast("Choose a protein FASTA smaller than 2 MB."); return; }
  $("sequence").value = await file.text(); sequenceCount(); event.target.value = ""; $("input-error").hidden = true;
});
$("example-button").addEventListener("click", async () => {
  try { const data = await (await api("example")).json(); $("sequence").value = data.sequence; sequenceCount(); $("input-error").hidden = true; toast("Loaded public actinorhodin cyclase CAC44203.1. Ready to search."); }
  catch (error) { toast(error.message); }
});
$("export-all").addEventListener("click", () => exportResults(null));
$("export-selected").addEventListener("click", () => exportResults([...state.selected]));
document.querySelector(".tabs").addEventListener("click", (event) => {
  const tab = event.target.closest("[data-tab]"); if (tab) switchTab(tab.dataset.tab);
});
document.querySelector(".tabs").addEventListener("keydown", (event) => {
  if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
  const names = ["hits", "occurrences", "contexts"];
  const index = (names.indexOf(state.tab) + (event.key === "ArrowRight" ? 1 : 2)) % 3;
  switchTab(names[index]); $("tab-" + names[index]).focus(); event.preventDefault();
});
document.querySelector(".preview-card").addEventListener("click", (event) => {
  const context = event.target.closest("[data-context]");
  const preview = event.target.closest("[data-preview]");
  const selection = event.target.closest("[data-selection]");
  const page = event.target.closest("[data-page]");
  if (context) { state.current = context.dataset.context; switchTab("contexts"); }
  if (preview) { state.current = preview.dataset.preview; renderContexts(); }
  if (selection) {
    state.selected = new Set(selection.dataset.selection === "all" ? state.run.occurrences.filter((o) => o.status === "success").map((o) => o.id) : []);
    renderContexts();
  }
  if (page && !page.disabled) {
    state[page.dataset.page + "Page"] += Number(page.dataset.delta);
    page.dataset.page === "hit" ? renderHits() : renderOccurrences();
  }
});
$("occurrences-content").addEventListener("keydown", (event) => {
  const row = event.target.closest("[data-context]");
  if (row && (event.key === "Enter" || event.key === " ")) { state.current = row.dataset.context; switchTab("contexts"); event.preventDefault(); }
});
$("contexts-content").addEventListener("change", (event) => {
  const checkbox = event.target.closest("[data-select]"); if (!checkbox) return;
  checkbox.checked ? state.selected.add(checkbox.dataset.select) : state.selected.delete(checkbox.dataset.select);
  $("selection-label").textContent = `${number(state.selected.size)} of ${number(state.run.counts.contexts)} selected · native genomic orientation`;
  updateControls();
});
$("history-button").addEventListener("click", async () => {
  try {
    const data = await (await api("runs")).json();
    $("history-content").innerHTML = data.runs.length ? data.runs.map((r) => `<button class="history-run" data-run="${esc(r.id)}"><span><strong>${esc(r.query)}</strong><small>${esc(new Date(r.created_at).toLocaleString())} · ${r.counts.hits} hits → ${r.counts.occurrences} occurrences → ${r.counts.contexts} contexts</small></span>${badge(r.status)}</button>`).join("") : '<p>No saved runs yet. Start with a protein query.</p>';
    $("history-dialog").showModal();
  } catch (error) { toast(error.message); }
});
$("close-history").addEventListener("click", () => $("history-dialog").close());
$("history-content").addEventListener("click", async (event) => {
  const button = event.target.closest("[data-run]"); if (!button) return;
  $("history-dialog").close(); await loadRun(button.dataset.run);
});

async function initialize() {
  renderRun();
  try {
    const config = await (await api("config")).json();
    state.desktop = config.desktop;
    $("key-status").innerHTML = `<i></i> ${config.api_key_configured ? "E-utilities API key configured" : "E-utilities · standard rate"}`;
    $("email").value = config.email;
    $("data-directory").textContent = config.data_directory;
    const runs = await (await api("runs")).json();
    const selected = new URLSearchParams(location.search).get("run");
    if (selected || runs.active_id || runs.runs.length) await loadRun(selected || runs.active_id || runs.runs[0].id);
  } catch (error) { toast(error.message); }
  setInterval(() => { if (state.run?.active) refresh(); }, 1500);
}
initialize();
