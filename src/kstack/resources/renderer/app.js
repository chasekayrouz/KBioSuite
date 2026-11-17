/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
/* Live rendering + interaction only. Exports consume the project through separate modules. */
(function () {
  'use strict';
  const M = window.CSModel, $ = id => document.getElementById(id), NS = 'http://www.w3.org/2000/svg';
  let project, geo, native = false, selected = null, raf = 0, dirtyRows = new Set(), rebuild = true;
  let redrawConnectors = true, toastTimer, viewport = {zoom: 1, x: 20, y: 20}, drag = null;
  const rowNodes = new Map(), metrics = {firstRenderMs: 0, lastEditMs: 0, connectorMs: 0, frames: 0};
  let viewportSize = {width: $('viewport').clientWidth, height: $('viewport').clientHeight};
  const svgEl = (tag, attrs = {}, text) => {
    const element = document.createElementNS(NS, tag);
    Object.entries(attrs).forEach(([key,value]) => element.setAttribute(key, value));
    if (text !== undefined) element.textContent = text;
    return element;
  };
  const el = (tag, cls, text) => {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (text !== undefined) node.textContent = text;
    return node;
  };
  function bridge(message) {
    if (native && window.webkit?.messageHandlers?.kstack) window.webkit.messageHandlers.kstack.postMessage(message);
  }
  function changed(reason) {
    project.metadata.modified = new Date().toISOString();
    bridge({action: 'changed', reason, project});
  }
  function toast(text) {
    $('toast').textContent = text; $('toast').hidden = false;
    clearTimeout(toastTimer); toastTimer = setTimeout(() => { $('toast').hidden = true; }, 2300);
  }
  function download(kind, content, filename, mime) {
    if (!content) { toast('Select at least one protein first.'); return; }
    if (native) { bridge({action: 'export', kind, content, filename}); return; }
    const url = URL.createObjectURL(new Blob([content], {type: mime}));
    const a = el('a'); a.href = url; a.download = filename; document.body.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 3000);
  }
  function invalidate(rows = [], full = false, connections = true) {
    rows.forEach(id => dirtyRows.add(id)); rebuild ||= full; redrawConnectors ||= connections;
    if (!raf) raf = requestAnimationFrame(render);
  }
  function transform() {
    const z = viewport.zoom;
    // The browser paints SVG primitives at the current viewport resolution.
    // No CSS-scaled scene, promoted layer, or bitmap connector intermediate.
    $('genes-svg').setAttribute('viewBox', `${-viewport.x/z} ${-viewport.y/z} ${viewportSize.width/z} ${viewportSize.height/z}`);
    for (const [cid, row] of geo?.rows || []) {
      const actions = $('row-controls').querySelector(`[data-cluster="${cid}"]`);
      if (!actions) continue;
      const controlZoom=z*Math.min(geo.sx,geo.sy);
      actions.style.zoom = controlZoom;
      actions.style.left = `${(10*geo.sy*z + viewport.x)/controlZoom}px`;
      actions.style.top = `${((row.y + 30*geo.sy)*z + viewport.y)/controlZoom}px`;
    }
    $('zoom-value').textContent = `${Math.round(viewport.zoom * 100)}%`;
  }
  function zoom(factor, center) {
    const box = $('viewport').getBoundingClientRect();
    const cx = center ? center.x - box.left : box.width / 2, cy = center ? center.y - box.top : box.height / 2;
    const next = Math.max(Math.min(.08, viewport.zoom), Math.min(12, viewport.zoom * factor));
    viewport.x = cx - (cx - viewport.x) * next / viewport.zoom;
    viewport.y = cy - (cy - viewport.y) * next / viewport.zoom;
    viewport.zoom = next; transform();
  }
  function fit() {
    if (!geo) return;
    const box = $('viewport').getBoundingClientRect();
    // Read the chosen geometry even when a scale edit is awaiting its frame.
    // Very wide stacks may need less than the usual 8% interactive zoom floor.
    viewport.zoom = Math.min(1.5, Math.max(1,box.width - 40) / M.geometry(project).width);
    viewport.x = 20; viewport.y = 18; transform();
  }
  function button(text, title, fn, cls = '') {
    const node = el('button', cls, text); node.title = title; node.setAttribute('aria-label',title);
    node.addEventListener('click', fn); return node;
  }
  function makeActions(cluster, row) {
    const actions = el('div','row-actions'); actions.dataset.cluster = cluster.id;
    if (native) {
      const handle = el('span','handle','⋮⋮'); handle.title = 'Drag to reorder cluster'; handle.setAttribute('role','button');
      handle.setAttribute('aria-label',`Reorder ${cluster.name}`); handle.tabIndex = 0;
      handle.addEventListener('pointerdown', event => startReorder(event, cluster.id));
      actions.append(handle);
      const up = button('↑', `Move ${cluster.name} up`, () => reorder(cluster.id,-1));
      const down = button('↓', `Move ${cluster.name} down`, () => reorder(cluster.id,1));
      const i=project.layout.order.indexOf(cluster.id);
      up.disabled=i===0; down.disabled=i===project.clusters.length-1;
      actions.append(up,down,button('↔',`Reverse ${cluster.name}`,()=>{
        project.layout.rows[cluster.id].orientation *= -1; changed('Reverse cluster'); invalidate([cluster.id]);
        if(selected) inspect(selected);
      }));
    }
    const proteins = cluster.cds.filter(g => g.protein);
    const all = proteins.length > 0 && proteins.every(g => project.fasta_selection.includes(g.id));
    const selection = button(all ? '− Proteins' : '+ Proteins', `Select all proteins in ${cluster.name}`, () => {
      M.selectCluster(project,cluster.id,!all); changed('Cluster FASTA selection'); invalidate([cluster.id]); updateSelection();
      if (selected) inspect(selected);
    }, 'select-cluster');
    selection.setAttribute('aria-pressed',all); selection.disabled=!proteins.length;
    actions.append(selection);
    if(cluster.warnings.length || cluster.cds.some(g=>g.warnings.length)) {
      actions.append(button('!', `Import notes for ${cluster.name}`,()=>inspectNotes(cluster.id)));
    }
    return actions;
  }
  function drawRow(cid) {
    const cluster = geo.index.clusters.get(cid), row = geo.rows.get(cid);
    let group=rowNodes.get(cid);
    if (!group) { group=svgEl('g',{'data-cluster':cid}); rowNodes.set(cid,group); $('gene-rows').append(group); }
    group.replaceChildren();
    const {sx,sy} = geo;
    const text = (x,y,size,attrs,value) => svgEl('text',{...M.textAttributes(x,y,size,geo),...attrs},value);
    if (native) {
      // Dense rows must not place a later row's transparent drag target over
      // the preceding row's gene glyphs.
      const hitHeight=Math.min(54*sy,geo.rowHeight);
      const zone=svgEl('rect',{x:row.x,y:row.y-hitHeight/2,width:row.end-row.x,height:hitHeight,class:'row-hitzone'});
      zone.addEventListener('pointerdown',event=>startRowDrag(event,cid)); group.append(zone);
    }
    group.append(svgEl('line',{x1:row.x,x2:row.end,y1:row.y,y2:row.y,stroke:'#c4ccd6','stroke-width':sy}));
    group.append(text(12*sy,row.y+3*sy,12,{class:'cluster-name'},cluster.name.length>28?cluster.name.slice(0,27)+'…':cluster.name));
    for(const gene of cluster.cds) {
      const shape=geo.genes.get(gene.id);
      shape.segments.forEach((segment,i)=>{
        if (i && !gene.origin_spanning) group.append(svgEl('line',{x1:shape.segments[i-1].x2,x2:segment.x1,y1:row.y,y2:row.y,
          stroke:M.color(project,gene,geo.index),'stroke-dasharray':`${3*sx} ${2*sx}`,'stroke-width':1.5*sy}));
        const path=svgEl('path',{d:M.arrowPath(segment),fill:M.color(project,gene,geo.index),stroke:'#485368','stroke-width':.7*sy,'stroke-linejoin':'round',
          class:`gene${selected===gene.id?' selected':''}${project.fasta_selection.includes(gene.id)?' in-fasta':''}`,
          'data-gene':gene.id,tabindex:0,role:'button','aria-haspopup':native?'menu':'false','aria-label':`${cluster.name}: ${M.label(gene)}`});
        path.addEventListener('click',event=>{event.stopPropagation();inspect(gene.id);});
        path.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();inspect(gene.id);}});
        path.addEventListener('pointerenter',event=>{
          $('tooltip').textContent=`${M.label(gene)}${gene.product?' · '+gene.product:''}${gene.origin_spanning?' · spans origin':''}`;
          $('tooltip').hidden=false; moveTooltip(event);
        });
        path.addEventListener('pointermove',moveTooltip);
        path.addEventListener('pointerleave',()=>{$('tooltip').hidden=true;});
        group.append(path);
      });
      if(gene.origin_spanning) {
        group.append(text(row.x-10*sy,row.y-15*sy,10,{fill:M.color(project,gene,geo.index)},'↔'));
        group.append(text(row.end+3*sy,row.y-15*sy,10,{fill:M.color(project,gene,geo.index)},'↔'));
      }
      if(project.visualization.labels&&shape.segments.length) {
        const segment=shape.segments.reduce((a,b)=>b.x2-b.x1>a.x2-a.x1?b:a);
        const label=M.label(gene);
        group.append(text((segment.x1+segment.x2)/2,row.y-18*sy,8,{'text-anchor':'middle',fill:'#536074'},label.length>26?label.slice(0,25)+'…':label));
      }
    }
    const previous=$('row-controls').querySelector(`[data-cluster="${cid}"]`);
    const actions=makeActions(cluster,row); if(previous)previous.replaceWith(actions);else $('row-controls').append(actions);
  }
  function moveTooltip(e) {
    $('tooltip').style.left=`${Math.min(e.clientX+13,window.innerWidth-310)}px`;
    $('tooltip').style.top=`${Math.max(8,e.clientY-35)}px`;
  }
  function drawConnectors() {
    const start=performance.now(), layer=$('connectors-svg'), connections=M.connectors(project,geo);
    // Reuse polygons during row dragging; no duplicate export tree or per-zoom rebuild.
    connections.forEach((connection,i)=>{
      let polygon=layer.children[i];
      if(!polygon){polygon=svgEl('polygon');layer.append(polygon);}
      polygon.setAttribute('points',connection.points.map(p=>p.join(',')).join(' '));
      polygon.setAttribute('fill',connection.color);
      polygon.dataset.query=connection.query;polygon.dataset.target=connection.target;
    });
    while(layer.children.length>connections.length)layer.lastElementChild.remove();
    metrics.connectorMs=performance.now()-start;
  }
  function render() {
    raf=0;const start=performance.now();const oldGeo=geo;geo=M.geometry(project);
    if(!oldGeo||geo.minimum!==oldGeo.minimum) {
      if(oldGeo) { viewport.x+=M.cameraShift(oldGeo,geo,viewport.zoom); }
      rebuild=true;
    }
    if(!oldGeo||geo.sx!==oldGeo.sx||geo.sy!==oldGeo.sy) {
      const style=$('genes-svg').style;
      style.setProperty('--gene-stroke',1.5*geo.sy);style.setProperty('--selected-stroke',2.3*geo.sy);
      style.setProperty('--fasta-dash',`${3*geo.sx} ${geo.sx}`);
    }
    if(rebuild) {
      $('gene-rows').replaceChildren();$('scale-svg').replaceChildren();$('row-controls').replaceChildren();rowNodes.clear();
      project.layout.order.forEach(drawRow);
      const {bp, x, y, end}=M.scaleBar(project,geo), sy=geo.sy;
      $('scale-svg').append(svgEl('line',{x1:x,x2:end,y1:y,y2:y,stroke:'#30405a','stroke-width':2*sy}),
        svgEl('text',{...M.textAttributes((x+end)/2,y+18*sy,10,geo),'text-anchor':'middle',fill:'#647187'},`${bp/1000} kb`));
      [x,end].forEach(x=>$('scale-svg').append(svgEl('line',{x1:x,x2:x,y1:y-4*sy,y2:y+4*sy,stroke:'#30405a','stroke-width':sy})));
    } else dirtyRows.forEach(drawRow);
    if(redrawConnectors)drawConnectors();
    transform();
    dirtyRows.clear();rebuild=false;redrawConnectors=false;metrics.lastEditMs=performance.now()-start;metrics.frames++;
    if(!metrics.firstRenderMs)metrics.firstRenderMs=metrics.lastEditMs;
    $('render-time').textContent=`${project.clusters.length} clusters · ${geo.index.genes.size} CDS · shared bp scale`;
  }
  function updateSelection() {
    $('selection-count').textContent=project.fasta_selection.length;
    $('export-fasta').disabled=!project.fasta_selection.length;
  }
  function inspect(gid) {
    const gene=geo.index.genes.get(gid);if(!gene)return;
    const previous=selected;selected=gid;
    if(previous)dirtyRows.add(geo.index.genes.get(previous).cluster_id);invalidate([gene.cluster_id]);
    const c=geo.index.clusters.get(gene.cluster_id),group=geo.index.groups.get(gene.group_id), container=$('inspector-content');
    container.replaceChildren();$('inspector').classList.remove('empty');
    const title=el('h2','gene-title',M.label(gene));title.style.color=M.color(project,gene,geo.index);container.append(title,el('p','product',gene.product||'No product annotation'));
    const actions=el('div','protein-actions');
    const copy=button('Copy Protein Sequence','Copy Protein Sequence',async()=>{
      if(native){bridge({action:'copy',content:gene.protein});toast('Protein copied');return;}
      try{await navigator.clipboard.writeText(gene.protein);}catch(error){
        const input=el('textarea');input.value=gene.protein;input.style.position='fixed';input.style.left='-9999px';document.body.append(input);input.select();
        const copied=document.execCommand('copy');input.remove();if(!copied){toast('Copy unavailable; select the sequence below.');return;}
      }toast('Protein copied');
    },'primary');copy.disabled=!gene.protein;
    const inSelection=project.fasta_selection.includes(gid);
    const select=button(inSelection?'Remove from FASTA Selection':'Add to FASTA Selection','Toggle gene FASTA selection',()=>{
      const ids=new Set(project.fasta_selection);inSelection?ids.delete(gid):ids.add(gid);project.fasta_selection=[...ids].sort();
      changed('Gene FASTA selection');updateSelection();inspect(gid);
    });select.disabled=!gene.protein;actions.append(copy,select);container.append(actions);
    const facts=el('dl','facts');
    const pairs=[['Cluster',c.name],['Gene',gene.gene||'—'],['Locus tag',gene.locus_tag||'—'],['Coordinates',M.coordinates(gene)],
      ['Strand',gene.strand===1?'+':gene.strand===-1?'−':'unknown'],['CDS length',`${gene.cds_length.toLocaleString()} bp`],
      ['Protein length',`${gene.protein_length.toLocaleString()} aa`],['Homolog group',group?.members.length>1?`${project.homolog_groups.indexOf(group)+1} · ${group.members.length} members`:'singleton'],
      ['Translation',gene.translation_source],['Internal ID',gene.id]];
    pairs.forEach(([key,value])=>{const row=el('div');row.append(el('dt','',key),el('dd','',value));facts.append(row);});container.append(facts);
    if(gene.warnings.length)container.append(el('p','warnings',gene.warnings.join(' ')));
    if(native) {
      container.append(el('div','section-label','GENE COLORS'));
      const controls=el('div','color-tools'),label=el('label','', 'Override color '),input=el('input');input.type='color';input.value=M.color(project,gene,geo.index);input.setAttribute('aria-label','Gene override color');label.append(input);controls.append(label);
      controls.append(button('Apply to gene','Apply color to this gene',()=>{M.setGeneColor(project,gid,input.value);changed('Gene color');invalidate([gene.cluster_id]);inspect(gid);}),
        button('Apply to group','Apply color to homolog group',()=>{M.setGroupColor(project,gene.group_id,input.value);changed('Group color');invalidate(project.clusters.map(c=>c.id));inspect(gid);}),
        button('Reset gene','Reset this gene override',()=>{M.resetColor(project,gene);changed('Reset gene color');invalidate([gene.cluster_id]);inspect(gid);}),
        button('Reset group','Reset group to generated colors',()=>{M.resetColor(project,gene,true);changed('Reset group color');invalidate(project.clusters.map(c=>c.id));inspect(gid);}));container.append(controls);
    }
    container.append(el('div','section-label','COMPLETE PROTEIN SEQUENCE'),el('pre','sequence',gene.protein||'No interpretable protein'));
    const details=el('details','qualifiers'),summary=el('summary','','GenBank qualifiers');details.append(summary);
    const dl=el('dl');for(const [key,values]of Object.entries(gene.qualifiers)){if(key==='translation')continue;dl.append(el('dt','',key),el('dd','',values.join('; ')));}details.append(dl);container.append(details);
  }
  function reorder(cid,delta) {
    if(!native)return;const before=[...project.layout.order];if(!M.moveRow(project,cid,delta))return;
    changed('Reorder clusters');invalidate(before.filter((id,i)=>id!==project.layout.order[i]));
  }
  function inspectNotes(cid) {
    const cluster=geo.index.clusters.get(cid), previous=selected, container=$('inspector-content');
    selected=null;container.replaceChildren();$('inspector').classList.remove('empty');
    if(previous)invalidate([geo.index.genes.get(previous).cluster_id]);
    container.append(el('h2','gene-title',cluster.name),el('p','product','Import notes — every CDS remains in this project.'));
    for(const note of cluster.warnings)container.append(el('p','warnings',note));
    for(const gene of cluster.cds.filter(g=>g.warnings.length)) {
      container.append(button(M.label(gene),`Inspect retained CDS ${gene.id}`,()=>inspect(gene.id)),el('p','warnings',gene.warnings.join(' ')));
    }
  }
  function startRowDrag(e,cid) {
    if(!native||e.button!==0)return;e.preventDefault();e.stopPropagation();
    $('viewport').setPointerCapture(e.pointerId);drag={kind:'row',cid,start:e.clientX,offset:project.layout.rows[cid].offset,moved:false};
  }
  function startReorder(e,cid) {
    if(!native||e.button!==0)return;e.preventDefault();e.stopPropagation();
    $('viewport').setPointerCapture(e.pointerId);drag={kind:'reorder',cid,start:e.clientY,index:project.layout.order.indexOf(cid),moved:false};
  }
  $('viewport').addEventListener('pointerdown',e=>{
    if(e.button!==0||e.target.closest('button,.gene,.handle'))return;
    e.preventDefault();$('viewport').setPointerCapture(e.pointerId);drag={kind:'pan',x:e.clientX,y:e.clientY,vx:viewport.x,vy:viewport.y};
  });
  window.addEventListener('pointermove',e=>{
    if(!drag)return;
    if(drag.kind==='pan'){viewport.x=drag.vx+e.clientX-drag.x;viewport.y=drag.vy+e.clientY-drag.y;transform();}
    if(drag.kind==='row'){
      const delta=(e.clientX-drag.start)/(viewport.zoom*geo.scale);
      if(Math.abs(delta)>1){drag.moved=true;project.layout.rows[drag.cid].offset=Math.round((drag.offset+delta)*1000)/1000;invalidate([drag.cid]);}
    }
    if(drag.kind==='reorder'){
      const next=Math.max(0,Math.min(project.clusters.length-1,Math.round(drag.index+(e.clientY-drag.start)/(viewport.zoom*geo.rowHeight))));
      const current=project.layout.order.indexOf(drag.cid);
      if(next!==current){const before=[...project.layout.order];project.layout.order.splice(current,1);project.layout.order.splice(next,0,drag.cid);drag.moved=true;invalidate(before.filter((id,i)=>id!==project.layout.order[i]));}
    }
  });
  window.addEventListener('pointerup',()=>{if(drag?.moved)changed(drag.kind==='row'?'Horizontal cluster offset':'Reorder clusters');drag=null;});
  window.addEventListener('pointercancel',()=>{if(drag?.moved)changed('Cluster layout');drag=null;});
  $('viewport').addEventListener('wheel',e=>{e.preventDefault();if(e.shiftKey){viewport.x-=e.deltaY;transform();}else zoom(Math.exp(-e.deltaY*.0015),{x:e.clientX,y:e.clientY});},{passive:false});
  $('viewport').addEventListener('keydown',e=>{if(['ArrowUp','ArrowDown','ArrowLeft','ArrowRight'].includes(e.key)){e.preventDefault();viewport.x+=e.key==='ArrowLeft'?35:e.key==='ArrowRight'?-35:0;viewport.y+=e.key==='ArrowUp'?35:e.key==='ArrowDown'?-35:0;transform();}});
  $('zoom-in').onclick=()=>zoom(1.2);$('zoom-out').onclick=()=>zoom(1/1.2);$('fit').onclick=fit;
  const presentationKeys = {width:'width_scale',height:'height_scale',spacing:'row_spacing_scale'};
  function setScale(axis, percent) {
    if (!Number.isFinite(percent)) return;
    const value=Math.max(25,Math.min(400,percent))/100, key=presentationKeys[axis];
    $(`${axis}-scale`).value=value*100;
    if ((project.visualization[key] ?? 1) === value) return;
    project.visualization[key]=value;changed(`${axis} scale`);invalidate([],true);
  }
  for (const axis of Object.keys(presentationKeys)) {
    const input=$(`${axis}-scale`);
    input.oninput=()=>{if(input.value!=='' && input.validity.valid)setScale(axis,Number(input.value));};
    // Incomplete/out-of-range typing must not leave a displayed value that
    // disagrees with the saved/exported factor after leaving the field.
    input.onchange=()=>{
      if(input.value!=='' && Number.isFinite(input.valueAsNumber))setScale(axis,Math.round(input.valueAsNumber/5)*5);
      else input.value=(project.visualization[presentationKeys[axis]] ?? 1)*100;
    };
  }
  $('reset-scale').onclick=()=>{for(const axis of Object.keys(presentationKeys))setScale(axis,100);};
  new ResizeObserver(([entry])=>{
    viewportSize={width:entry.contentRect.width,height:entry.contentRect.height};transform();
  }).observe($('viewport'));
  $('connectors').onchange=()=>{project.visualization.connectors=$('connectors').checked;changed('Connectors');invalidate();};
  $('labels').onchange=()=>{project.visualization.labels=$('labels').checked;changed('Gene labels');invalidate(project.layout.order,false,false);};
  $('display-identity').onchange=()=>{const value=Number($('display-identity').value);if(!Number.isFinite(value))return;project.visualization.connector_identity=Math.max(0,Math.min(100,value))/100;$('display-identity').value=project.visualization.connector_identity*100;changed('Connector display threshold');invalidate();};
  const exportPreview = window.CSExportPreview.create({getProject:()=>project,
    onExport:content=>download('svg',content,'Kstack.svg','image/svg+xml')});
  const geneMenu = window.CSGeneMenu.attach({surface:$('genes-svg'),host:document.body,enabled:()=>native,
    onAlign:gid=>{
      const offsets=M.alignmentOffsets(project,gid), rows=Object.keys(offsets);
      if(!rows.length)return;
      for(const cid of rows)project.layout.rows[cid].offset=offsets[cid];
      changed('Auto-align gene group');invalidate(rows);
    }});
  function openExportPreview() { geneMenu.close();exportPreview.open(); }
  $('export-svg').onclick=openExportPreview;
  $('export-fasta').onclick=()=>download('faa',window.CSFASTA.generate(project),'Kstack.faa','text/plain');
  $('close-inspector').onclick=()=>{const previous=selected;selected=null;$('inspector').classList.add('empty');$('inspector-content').replaceChildren(el('p','empty-inspector','Select a gene to inspect its protein and annotations.'));if(previous)invalidate([geo.index.genes.get(previous).cluster_id]);};
  function setProject(p) {
    if(p.format!=='Kstack'||p.schema_version!==1)throw new Error('Unsupported Kstack project schema.');
    geneMenu.close();exportPreview.close();
    project=JSON.parse(JSON.stringify(p));selected=null;geo=null;metrics.firstRenderMs=0;
    $('project-title').textContent=project.metadata.title||'Your stack';
    $('summary').textContent=`${project.clusters.length} clusters · ${Math.round(project.analysis.identity*100)}% identity · ${Math.round(project.analysis.coverage*100)}% coverage`;
    $('connectors').checked=project.visualization.connectors;$('labels').checked=project.visualization.labels;
    for(const [axis,key] of Object.entries(presentationKeys))$(`${axis}-scale`).value=(project.visualization[key] ?? 1)*100;
    $('display-identity').value=Math.round(project.visualization.connector_identity*100);updateSelection();
    $('inspector').classList.add('empty');$('inspector-content').replaceChildren(el('p','empty-inspector','Select any gene to inspect its annotations and complete protein sequence.'));
    invalidate([],true);requestAnimationFrame(()=>{fit();bridge({action:'ready',metrics});});
  }
  window.Kstack={setProject,getProject:()=>JSON.parse(JSON.stringify(project)),
    setNativeMode(value){native=!!value;document.body.classList.toggle('native',native);$('interaction-hint').textContent=native?'Drag a row to align · Handles reorder · Right-click a gene to auto-align':'Click a gene to inspect · Drag background to pan · Scroll to zoom';invalidate([],true);},
    exportSVG:()=>window.CSSVG.generate(project),exportFASTA:()=>window.CSFASTA.generate(project),
    openExportPreview,
    inspect,inspectNotes,fit,zoom,metrics,getViewport:()=>({...viewport}),
    selectCluster(cid,value=true){M.selectCluster(project,cid,value);changed('Cluster FASTA selection');updateSelection();invalidate([cid]);},
    selectGene(gid,value=true){const gene=M.index(project).genes.get(gid);if(!gene?.protein)return;const selection=new Set(project.fasta_selection);value?selection.add(gid):selection.delete(gid);project.fasta_selection=[...selection].sort();changed('Gene FASTA selection');updateSelection();invalidate([gene.cluster_id]);},
    edit(action,cid,value){if(!native)throw new Error('Structural editing is available in the desktop app.');if(action==='reverse')project.layout.rows[cid].orientation*=-1;if(action==='offset')project.layout.rows[cid].offset=Number(value);if(action==='move')M.moveRow(project,cid,Number(value));changed(action);invalidate([],true);}
  };
  try{setProject(JSON.parse($('kstack-project').textContent));}catch(error){$('project-title').textContent='Unable to open this stack';$('inspector-content').textContent=error.message;}
})();
