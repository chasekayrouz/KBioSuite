'use strict';
const $ = (id) => document.getElementById(id);
const token = document.querySelector('meta[name=kmatch-token]').content;
let job = null, desktop = false, dependencies = [], viewer = null, loadedRun = null;
let customMsa = null, validationTimer = null, started = null, previousBusy = false;
let defaultsLoaded = false, backendDefinitions = new Map();
let requestedDetails = null;
let hitOrder = [];
let selectedHit = null, alignmentSupport = null, exportDialogOpen = false;
async function api(path, body) {
  const response = await fetch('/api/' + path, {method: body === undefined ? 'GET' : 'POST',
    headers: {'X-Kmatch-Token': token, 'Content-Type': 'application/json'},
    body: body === undefined ? undefined : JSON.stringify(body)});
  if (!response.ok) {const error = await response.json(); throw new Error(error.error || 'Kmatch could not complete the operation.');}
  return response;
}
function error(message) {$('error').textContent = message; $('error').hidden = !message;}
function percent(value) {return (100 * value).toFixed(1) + '%';}
function number(value, digits=3) {return value == null ? 'Unavailable' : Number(value).toFixed(digits);}
async function refreshDependencies() {
  const response = await api('dependencies'); const body = await response.json();
  dependencies = body.items;
  $('dependencies').replaceChildren();
  const labels = ['Boltz-2', 'AlphaFold 3', 'Foldseek', 'PDB100'];
  const components = ['boltz', 'af3-code', 'foldseek', 'pdb100'];
  dependencies.forEach((state, index) => {
    const card = document.createElement('div'); card.className = 'dependency';
    const title = document.createElement('strong'); title.textContent = labels[index];
    const readiness = document.createElement('p'); readiness.className = state.ready ? 'ready' : 'missing';
    readiness.textContent = state.message;
    const button = document.createElement('button'); button.textContent = index === 1 ? 'Install AF3 code' : state.ready ? (index === 3 ? 'Check / update' : 'Verify / setup') : 'Install';
    button.dataset.component = components[index]; button.disabled = job?.busy || false;
    button.addEventListener('click', async () => {
      try {error(''); await api('setup', {component: components[index]}); started = Date.now(); await poll();} catch(e) {error(e.message);}
    });
    card.append(title, readiness, button);
    if (index === 1) {
      const help = document.createElement('p'); help.textContent = 'Weights remain external. Obtain af3.bin.zst directly from Google under its separate terms.';card.append(help);
    }
    if (state.version) {const version = document.createElement('small'); version.textContent = 'Version: ' + state.version; card.append(version);}
    $('dependencies').append(card);
  });
  const ready = dependencies[0]?.ready && dependencies[2]?.ready && dependencies[3]?.ready;
  $('dependency-summary').textContent = ready ? 'Boltz-2 & PDB100 ready' : 'Setup required';
  if (!ready) $('setup-panel').open = true;
}
function nativeApi() {return window.kmatchNative();}
function nativeReady() {return desktop && window.kmatchBridge.ready;}
async function loadStructure(prediction, runId) {
  if (runId === loadedRun || !prediction) return;
  const structure = await (await api('structure')).text();
  try {
    if (!viewer) viewer = $3Dmol.createViewer($('viewer'), {backgroundColor: '#f9fbfa', antialias: true});
    viewer.removeAllModels(); viewer.addModel(structure, prediction.structure_path.endsWith('.cif') ? 'cif' : 'pdb');
    loadedRun = runId; $('viewer-empty').hidden = true; styleStructure(); viewer.zoomTo(); viewer.render();
    $('confidence-legend').hidden = !prediction.per_residue_confidence?.length;
  } catch(e) {error('3D view could not initialize: ' + e.message + '. The structure can still be saved.');}
}
function styleStructure() {
  if (!viewer) return;
  const confidence = $('color-mode').value === 'confidence';
  viewer.setStyle({}, {cartoon: {colorfunc: (atom) => confidence ? (atom.b >= 90 ? '#0053d6' : atom.b >= 70 ? '#65cbf3' : atom.b >= 50 ? '#ffdb13' : '#ff7d45') : '#377b57'}});
  viewer.render();
}
function selectHit(index) {
  const hit = job.hits[index]; if (!hit) return;
  selectedHit = index;
  document.querySelectorAll('#hits tr').forEach(row => row.classList.toggle('selected', Number(row.dataset.hitIndex) === index));
  $('hit-detail').hidden = false; $('hit-title').textContent = `#${hitOrder.indexOf(index)+1} · ${hit.target}`;
  const fields = [['Query TM-score', number(hit.query_tm)], ['Target TM-score',number(hit.target_tm)], ['Alignment TM-score',number(hit.alignment_tm)], ['LDDT',number(hit.lddt)],
    ['Homology probability',percent(hit.probability)], ['E-value',hit.evalue.toExponential(2)], ['Bit score',number(hit.bits,1)], ['Sequence identity',percent(hit.sequence_identity)],
    ['Query coverage',percent(hit.query_coverage)], ['Target coverage',percent(hit.target_coverage)], ['Query coordinates',`${hit.query_start}–${hit.query_end} / ${hit.query_length}`], ['Target coordinates',`${hit.target_start}–${hit.target_end} / ${hit.target_length}`]];
  $('metrics').replaceChildren(); fields.forEach(([label,value]) => {const node = document.createElement('div');node.className='metric';const title=document.createElement('span');title.textContent=label;node.append(title, document.createTextNode(value));$('metrics').append(node);});
  $('alignment').textContent = 'Query  ' + hit.query_alignment + '\nTarget ' + hit.target_alignment;
  $('transform').textContent = 'Foldseek rotation u (row-major; target → query):\n' + hit.rotation.map(row => row.join(', ')).join('\n') + '\nTranslation t: ' + hit.translation.join(', ') + '\nUsed directly for 3D alignment export; no refitting.';
  renderAlignmentExport(); renderPdbDetails();
}
let renderedHits = null;
function renderHits() {
  const signature = job.run_id + ':' + job.hits.length + ':' + $('hit-sort').value;
  if (signature === renderedHits) return;
  const previousSelection = renderedHits?.startsWith(job.run_id + ':') ? selectedHit : null;
  renderedHits = signature;
  selectedHit = null;
  hitOrder = job.hit_rankings?.[$('hit-sort').value] || [];
  $('hits').replaceChildren();$('hit-detail').hidden = true;
  hitOrder.forEach((index,rank) => {
    const hit = job.hits[index];
    const row = document.createElement('tr'); row.dataset.hitIndex = index;
    const values = [rank+1,hit.target,'Loading…',number(hit.query_tm),percent(hit.probability),hit.evalue.toExponential(1),percent(hit.query_coverage),hit.alignment_length+' aa',percent(hit.sequence_identity)];
    values.forEach((value,column) => {const cell = document.createElement('td');if (column === 1) {const button=document.createElement('button');button.textContent=value;button.addEventListener('click',()=>selectHit(index));cell.append(button);
      const match=hit.target.match(/^([0-9][a-z0-9]{3})(?:-assembly[0-9]+)?_/i);
      if (match) {const link=document.createElement('a');link.href='https://www.rcsb.org/structure/'+match[1].toUpperCase();link.textContent='RCSB ↗';link.target='_blank';link.rel='noopener noreferrer';link.className='pdb-link';link.addEventListener('click',async event=>{event.stopPropagation();if(desktop){event.preventDefault();try{await nativeApi().open_pdb(job.run_id,index);}catch(e){error(e.message);}}});cell.append(link);}
    }else {cell.textContent=value;if(column===2){cell.className='protein-name';cell.dataset.target=hit.target;}}row.append(cell);});
    row.addEventListener('click',()=>selectHit(index));$('hits').append(row);
  });
  if (previousSelection !== null && hitOrder.includes(previousSelection)) selectHit(previousSelection);
  $('results').hidden = !job.hits.length; $('results-empty').hidden = !!job.hits.length;
  if (job.phase === 'complete' && job.prediction && !job.hits.length) $('results-empty').textContent='No hits met the search criteria. The prediction and search diagnostics are retained in the run folder.';
  else $('results-empty').textContent='Homologs will appear after prediction and the local Foldseek search.';
  $('hit-count').textContent = job.hits.length ? `· ${hitOrder.length} of ${job.hits.length} hits` : '';
  if(job.hits.length && requestedDetails!==job.run_id){requestedDetails=job.run_id; requestPdbDetails();}
}
async function requestPdbDetails() {
  if(!job?.hits.length) return;
  try {await api('pdb-details',{run_id:job.run_id,allow_network:$('pdb-network').checked});}
  catch(e){$('pdb-status').textContent=e.message;}
}
function renderPdbDetails() {
  $('reload-pdb').disabled=!job?.hits.length;
  $('pdb-status').textContent=job.pdb_details_status || '';
  document.querySelectorAll('.protein-name').forEach(cell=>{
    const data=job.pdb_details?.[cell.dataset.target];
    cell.textContent=data?.protein_name || (data?.status==='ready'?'Chain name unavailable':data?.status==='unavailable'?'Unavailable':'Loading…');
    cell.title=data?.protein_name || data?.message || '';
  });
  if(selectedHit===null)return;
  const data=job.pdb_details?.[job.hits[selectedHit]?.target];
  $('ligands').replaceChildren();
  $('ligand-status').textContent=data?.status==='ready' ?
    (data.ligands.length?`Entry-wide annotation · cached ${data.retrieved_utc} · not every ligand necessarily contacts this hit chain.`:'No non-water small-molecule ligands reported in this entry.') :
    (data?.message || 'Loading public PDB ligand details…');
  (data?.ligands || []).forEach(ligand=>{const item=document.createElement('li');item.textContent=`${ligand.id} · ${ligand.name}`+(ligand.formula?` · ${ligand.formula}`:'')+` · ${ligand.copies} instance(s)`;$('ligands').append(item);});
}
function renderAlignmentExport() {
  const busy = job.busy || exportDialogOpen;
  const exporting = job.busy && job.operation === 'alignment_export';
  const option = $('alignment-format').querySelector('option[value="pse"]');
  option.disabled = !alignmentSupport?.pse;
  if (option.disabled && $('alignment-format').value === 'pse') $('alignment-format').value = 'cif';
  $('alignment-format').disabled = busy; $('full-atom').disabled=busy;
  $('export-alignment').disabled = busy || selectedHit === null || (desktop && !nativeReady());
  $('export-alignment').textContent = exporting ? 'Exporting…' : 'Export 3D alignment…';
  $('cancel-alignment').hidden = !exporting;
  $('alignment-progress').hidden = !exporting && !(job.operation === 'alignment_export' && job.error);
  $('alignment-progress').textContent = job.error || job.message;
  $('alignment-note').textContent = ($('full-atom').checked ? 'Full-atom PDB hit: matching ribbons, contacting ligands and aligned binding-site residues. '+($('pdb-network').checked?'Missing coordinates download from RCSB and are cached.':'Cached coordinates only.') : 'Local Cα trace only: no ligands or hit side chains; target remains a trace.') + ' ' + (alignmentSupport?.message || '');
  const result = job.alignment_export;
  $('alignment-result').hidden = !result;
  if (result) {
    $('alignment-result-text').textContent = `${result.target} · ${result.format.toUpperCase()} · aligned Cα RMSD ${number(result.metrics.aligned_ca_rmsd_after_angstrom,2)} Å` + (result.saved_to ? ` · Saved to ${result.saved_to}` : ' · Ready to download');
    $('download-alignment').hidden = desktop;
    $('download-alignment').disabled = busy;
  }
}
async function poll() {
  try {
    const response = await api('status'); const data = await response.json();
    job = data.job; desktop = data.desktop; alignmentSupport = data.alignment_support; window.setNativeState(desktop);
    if (!defaultsLoaded) {backendDefinitions=new Map(data.backends.map(item=>[item.id,item])); $('backend').replaceChildren(...data.backends.map(item=>{const option=document.createElement('option');option.value=item.id;option.textContent=item.label;return option;})); $('backend').value=backendDefinitions.has(data.defaults.predictor)?data.defaults.predictor:data.backends[0].id; backendNote(); defaultsLoaded=true;}
    $('asset-root').textContent = data.runs === data.home + '/runs' ? 'Assets & run outputs: ' + data.home : 'Assets: ' + data.home + ' · Run outputs: ' + data.runs;
    $('status').textContent = job.message; $('status-dot').classList.toggle('busy',job.busy);
    // A failed setup (e.g. missing Command Line Tools) is also shown inside the setup panel,
    // next to the button that started it; the status card is further down the page.
    const setupFailed = !job.busy && job.operation === 'setup' && Boolean(job.error);
    $('setup-error').textContent = setupFailed ? job.error : ''; $('setup-error').hidden = !setupFailed;
    $('run').disabled = job.busy || exportDialogOpen; $('cancel').hidden = !job.busy;
    ['sequence','backend','msa-mode','choose-msa'].forEach(id=>$(id).disabled=job.busy);
    document.querySelectorAll('[data-component]').forEach(button=>button.disabled=job.busy);
    $('export').disabled = !job.prediction || (desktop && !nativeReady()); $('reveal').disabled = !job.directory || !nativeReady();
    $('choose-msa').disabled = job.busy || !nativeReady();
    if (job.error) error(job.error);
    $('elapsed').textContent = job.busy && started ? Math.floor((Date.now()-started)/1000)+' s elapsed' : '';
    document.querySelectorAll('[data-phase]').forEach(element=>element.classList.toggle('active',element.dataset.phase===job.phase));
    const warnings = job.prediction?.warnings || [];
    $('prediction-warning').textContent = warnings.map(item=>item.message).join(' ');
    $('prediction-warning').hidden = !warnings.length;
    if (job.prediction) {
      const p=job.prediction, c=p.confidence;
      $('prediction-summary').textContent = `${backendDefinitions.get(p.predictor)?.name || p.predictor} ${p.predictor_version} · ${p.sequence_length} residues · mean pLDDT ${number(c.mean_plddt,1)} · pTM ${number(c.ptm)} · ${Math.round(p.runtime_seconds)} s`;
      await loadStructure(p,job.run_id);
    } else {
      if (viewer && loadedRun) {viewer.removeAllModels();viewer.render();loadedRun=null;}
      $('viewer-empty').hidden=false;$('confidence-legend').hidden=true;
      $('prediction-summary').textContent='Local inference uses the Apple GPU. Confidence is a model estimate.';
    }
    renderHits();
    renderAlignmentExport(); renderPdbDetails();
    if (previousBusy && !job.busy) await refreshDependencies(); previousBusy=job.busy;
  } catch(e) {error(e.message);}
}
$('sequence').addEventListener('input',()=>{clearTimeout(validationTimer);validationTimer=setTimeout(async()=>{try {const result=await (await api('validate',{sequence:$('sequence').value})).json();$('length').textContent=result.length+' residues';$('validation').textContent='Valid single protein sequence';$('validation').classList.remove('error');}catch(e){$('length').textContent='Check sequence';$('validation').textContent=e.message;$('validation').classList.add('error');}},250);});
$('run').addEventListener('click',async()=>{try {error('');const mode=$('msa-mode').value;if(mode==='custom'&&!customMsa)throw new Error('Choose a local A3M alignment first.');await api('run',{sequence:$('sequence').value,backend:$('backend').value,network_msa:mode==='network',msa_path:mode==='custom'?customMsa:null});started=Date.now();await poll();}catch(e){error(e.message);}});
$('cancel').addEventListener('click',async()=>{try{await api('cancel',{});}catch(e){error(e.message);}});
$('setup-toggle').addEventListener('click',()=>{$('setup-panel').open=!$('setup-panel').open;});
$('msa-mode').addEventListener('change',()=>{const mode=$('msa-mode').value;$('custom-msa').hidden=mode!=='custom';$('msa-note').textContent=mode==='network'?'This mode sends your sequence to api.colabfold.com for sequence searching. The MSA is cached; structure prediction and Foldseek run locally.':mode==='custom'?'Choose an existing local A3M whose first sequence matches this query.':'MSA-free accuracy varies. An alignment can substantially improve difficult folds.';});
function backendNote() {$('backend-note').textContent=backendDefinitions.get($('backend').value)?.note || '';}
$('backend').addEventListener('change',backendNote);
$('choose-msa').addEventListener('click',async()=>{try{customMsa=await nativeApi().choose_msa();$('msa-name').textContent=customMsa?customMsa.split('/').pop():'No alignment selected';}catch(e){error(e.message);}});
$('color-mode').addEventListener('change',styleStructure);
$('reset-view').addEventListener('click',()=>{if(viewer){viewer.zoomTo();viewer.render();}});
$('export').addEventListener('click',async()=>{try{if(desktop){await nativeApi().export_structure();}else{const response=await api('export');const url=URL.createObjectURL(await response.blob());const link=document.createElement('a');link.href=url;link.download=job.prediction.structure_path.split('/').pop()+(job.prediction.predictor==='af3'?'.zip':'');link.click();URL.revokeObjectURL(url);}}catch(e){error(e.message);}});
$('reveal').addEventListener('click',async()=>{try{await nativeApi().reveal_run();}catch(e){error(e.message);}});
$('export-alignment').addEventListener('click',async()=>{
  if (selectedHit === null || job.busy || exportDialogOpen) return;
  const runId = job.run_id, hitIndex = selectedHit, format = $('alignment-format').value;
  const fullAtom=$('full-atom').checked, allowNetwork=$('pdb-network').checked;
  try {
    error(''); exportDialogOpen = true; renderAlignmentExport(); $('run').disabled = true;
    const began = desktop ? await nativeApi().export_alignment(runId, hitIndex, format, fullAtom, allowNetwork) :
      await api('alignment-export', {run_id:runId, hit_index:hitIndex, format, full_atom:fullAtom, allow_network:allowNetwork});
    if (began) started = Date.now();
  } catch(e) {error(e.message);}
  finally {exportDialogOpen = false; await poll();}
});
$('download-alignment').addEventListener('click',async()=>{
  try {
    const result=job.alignment_export; if (!result) return;
    const response=await api('alignment-export/'+encodeURIComponent(result.export_id));
    const url=URL.createObjectURL(await response.blob()); const link=document.createElement('a');
    link.href=url; link.download=`Kmatch_${result.target.replace(/[^A-Za-z0-9_]/g,'_')}_alignment.${result.format}`;
    link.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
  } catch(e) {error(e.message);}
});
$('cancel-alignment').addEventListener('click',async()=>{try{await api('cancel',{});}catch(e){error(e.message);}});
new ResizeObserver(()=>{if(viewer){viewer.resize();viewer.render();}}).observe($('viewer'));
(async()=>{await poll();try{await refreshDependencies();}catch(e){error(e.message);}setInterval(poll,800);})();

$('reload-pdb').addEventListener('click',requestPdbDetails);
$('pdb-network').addEventListener('change',()=>{requestPdbDetails();renderAlignmentExport();});
$('full-atom').addEventListener('change',renderAlignmentExport);

$('hit-sort').addEventListener('change',()=>{renderHits();renderPdbDetails();});
