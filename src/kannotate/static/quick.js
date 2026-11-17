'use strict';
(() => {
  const el=id=>document.getElementById(id), model=QuickModel;
  let result=null, revision=0, pendingRevision=0, ready=false, busy=false, valid=false, importing=false;
  let validationTimer, validationGeneration=0, source='All', weaker=false, sort='score', ascending=false;
  let selection=new model.Selection(), memory=new model.ViewMemory(), layout=null;
  let inputOpen=true, resourcesOpen=false, detailsOpen=false, alignmentIndex=0, detailKey='';
  const number=value=>value===0?'0':value.toExponential(2);
  function frameMode(){
    const quick=!el('quick-workspace').hidden;
    document.body.classList.toggle('quick-focus',!!result&&quick&&!inputOpen&&!resourcesOpen);
    el('quick-input-panel').hidden=!inputOpen;
    el('quick-edit').hidden=el('quick-resources').hidden=!result;
    el('quick-edit').textContent=inputOpen?'Back to results':'Edit FASTA';
    el('quick-resources').textContent=resourcesOpen?'Back to results':'Local databases';
  }
  function switchWorkspace(quick){
    el('quick-workspace').hidden=!quick;el('file-workspace').hidden=quick;
    for(const [id,active] of [['open-quick',quick],['open-files',!quick]]){
      el(id).classList.toggle('active',active);el(id).setAttribute('aria-pressed',String(active));
    }
    frameMode();
    if(quick){if(inputOpen)el('quick-input').focus();else renderArchitecture();}
  }
  el('open-quick').onclick=()=>switchWorkspace(true);
  el('open-files').onclick=()=>switchWorkspace(false);
  el('quick-edit').onclick=()=>{inputOpen=!inputOpen;resourcesOpen=false;frameMode();if(inputOpen)el('quick-input').focus();};
  el('quick-resources').onclick=()=>{resourcesOpen=!resourcesOpen;frameMode();};
  function controls(){
    el('quick-search').disabled=busy||importing||!ready||!valid;el('quick-clear').disabled=busy||importing;
    el('quick-file').disabled=el('quick-input').disabled=busy||importing;
  }
  el('quick-file').addEventListener('change',async()=>{
    const file=el('quick-file').files[0];if(!file||busy||importing)return;
    importing=true;++validationGeneration;clearTimeout(validationTimer);controls();
    el('quick-error').textContent='';el('quick-import-status').textContent=`Loading ${file.name}…`;
    el('quick-import-status').classList.remove('failure');
    try{
      const data=await api('quick/import',{text:await file.text(),filename:file.name});
      el('quick-input').value=data.fasta;
      el('quick-import-status').textContent=[`Loaded ${file.name}.`,...data.warnings].join(' ');
    }catch(error){el('quick-import-status').textContent=error.message;el('quick-import-status').classList.add('failure');}
    finally{importing=false;el('quick-file').value='';await validate();controls();}
  });
  async function validate(){
    const generation=++validationGeneration, text=el('quick-input').value;
    valid=false;controls();
    if(!text.trim()){el('quick-validation').textContent='Paste 1–20 proteins. Raw sequence gets the name query_1.';el('quick-validation').classList.remove('failure');return;}
    try{
      const data=await api('quick/validate',{text});
      if(generation!==validationGeneration)return;
      valid=true;el('quick-validation').textContent=`${data.count} protein${data.count===1?'':'s'} · ${data.lengths.reduce((a,b)=>a+b,0).toLocaleString()} residues`;
      el('quick-validation').classList.remove('failure');
    }catch(error){if(generation!==validationGeneration)return;el('quick-validation').textContent=error.message;el('quick-validation').classList.add('failure');}
    controls();
  }
  el('quick-input').addEventListener('input',()=>{
    el('quick-import-status').textContent='';
    ++validationGeneration;valid=false;controls();clearTimeout(validationTimer);
    validationTimer=setTimeout(validate,120);
  });
  async function start(){
    if(busy||importing)return;
    clearTimeout(validationTimer);await validate();if(!valid||!ready||busy||importing)return;
    try{
      el('quick-error').textContent='';busy=true;controls();
      el('quick-running').textContent='Starting local Pfam + CDD search…';
      await api('quick/search',{text:el('quick-input').value});await refresh();
    }catch(error){busy=false;controls();el('quick-error').textContent=error.message;}
  }
  el('quick-search').onclick=start;
  el('quick-input').addEventListener('keydown',event=>{
    if((event.metaKey||event.ctrlKey)&&event.key==='Enter'){event.preventDefault();start();}
  });
  el('quick-clear').onclick=async()=>{
    try{const response=await api('quick/clear',{});revision=response.revision;}
    catch(error){el('quick-error').textContent=error.message;return;}
    ++validationGeneration;clearTimeout(validationTimer);el('quick-input').value='';
    el('quick-file').value='';el('quick-import-status').textContent='';
    result=null;selection=new model.Selection();memory=new model.ViewMemory();layout=null;
    valid=false;inputOpen=true;resourcesOpen=false;el('quick-results').hidden=true;frameMode();
    el('quick-error').textContent='';el('quick-running').textContent='';hideTooltip();
    el('quick-validation').classList.remove('failure');validate();el('quick-input').focus();
  };
  window.quickState=async(state, databasesReady)=>{
    ready=databasesReady;busy=state.job.busy;controls();
    el('quick-running').textContent=busy||state.job.stopped?state.job.message:'';
    el('quick-running').classList.toggle('running',busy);
    if(state.job.error)el('quick-error').textContent=state.job.error;
    if(state.job.quick_revision!==revision&&state.job.quick_revision!==pendingRevision){
      pendingRevision=state.job.quick_revision;
      try{
        const response=await api('quick/results');revision=response.revision;
        result=response.result;selection=new model.Selection();memory=new model.ViewMemory();detailKey="";alignmentIndex=0;
        if(result)renderResult();else{el('quick-results').hidden=true;inputOpen=true;frameMode();}
      }catch(error){el('quick-error').textContent=error.message;}
      finally{pendingRevision=0;}
    }
  };
  function renderResult(){
    showDetails(false);
    inputOpen=false;resourcesOpen=false;frameMode();el('quick-results').hidden=false;
    const entries=Object.values(result.provenance).flatMap(p=>p.queries), cached=entries.filter(p=>p.cached).length;
    el('quick-summary').textContent=`${result.queries.length} proteins · Pfam ${result.database_versions.pfam} + CDD ${result.database_versions.cdd} · ${cached===entries.length&&entries.length?'cached · ':''}${result.timings.total.toFixed(2)} s engines`;
    el('quick-provenance').textContent=JSON.stringify({threads:result.threads,timings:result.timings,
      engines:Object.fromEntries(Object.entries(result.provenance).map(([name,p])=>[name,{
        database:p.identity?.database,executable:p.identity?.executable?.[0],version:p.identity?.executable_version,
        parameters:p.identity?.parameters,output_format:p.identity?.output_format,queries:p.queries,
      }])),
    },null,2);
    el('quick-warnings').textContent=result.warnings.join(' · ');
    el('quick-queries').replaceChildren(...result.queries.map(query=>{
      const button=document.createElement('button');button.dataset.query=query.id;button.title=query.name;
      const name=document.createElement('strong');name.textContent=query.name;
      const meta=document.createElement('span');meta.textContent=`${query.length.toLocaleString()} aa · ${query.hits.filter(h=>h.significant).length} hits`;
      button.append(name,meta);button.onclick=()=>selectQuery(query);return button;
    }));
    selectQuery(result.queries[0]);window.scrollTo(0,0);
  }
  function saveView(){memory.save(selection.query,source,weaker,{graph:el('architecture-scroll').scrollTop,table:el('hits-panel').scrollTop,details:el('details-panel').scrollTop,alignment:alignmentIndex});}
  function restoreView(){
    const state=memory.restore(selection.query,source,weaker);
    el('architecture-scroll').scrollTop=state.graph;el('hits-panel').scrollTop=state.table;el('details-panel').scrollTop=state.details;
  }
  function selectQuery(query){
    saveView();selection.selectQuery(query);alignmentIndex=memory.restore(query,source,weaker).alignment||0;
    el('query-name').textContent=query.name;el('query-name').title=query.name;
    const index=result.queries.indexOf(query);
    el('query-length').textContent=`${query.length.toLocaleString()} amino acids`;
    el('query-position').textContent=`${index+1} / ${result.queries.length} proteins`;
    el('previous-query').disabled=index===0;el('next-query').disabled=index===result.queries.length-1;
    for(const button of el('quick-queries').children){
      const selected=button.dataset.query===query.id;button.classList.toggle('active',selected);
      button.setAttribute('aria-pressed',String(selected));button.tabIndex=selected?0:-1;
    }
    el('copy-status').textContent='';hideTooltip();renderHits();restoreView();
  }
  function moveQuery(delta){
    if(!result)return;
    const index=Math.max(0,Math.min(result.queries.length-1,result.queries.indexOf(selection.query)+delta));
    selectQuery(result.queries[index]);
    const button=el('quick-queries').children[index];button.focus({preventScroll:true});reveal(el('quick-queries'),button);
  }
  el('previous-query').onclick=()=>moveQuery(-1);el('next-query').onclick=()=>moveQuery(1);
  el('quick-queries').onkeydown=event=>{
    if(['ArrowDown','ArrowUp','Home','End'].includes(event.key)){
      event.preventDefault();moveQuery({ArrowDown:1,ArrowUp:-1,Home:-20,End:20}[event.key]);
    }
  };
  function hits(){return selection.query?model.visibleHits(selection.query,source,weaker):[];}
  function currentHit(){return hits().find(h=>h.id===selection.hit)||null;}
  // Scroll only the intended pane. scrollIntoView also moves the outer workspace.
  function reveal(container,node,inset=0){
    if(!node)return;
    const box=node.getBoundingClientRect(),frame=container.getBoundingClientRect();
    if(box.top<frame.top+inset)container.scrollTop+=box.top-frame.top-inset;
    else if(box.bottom>frame.bottom)container.scrollTop+=box.bottom-frame.bottom;
  }
  function activate(id,fromGraph=false,span=null){
    const hit=selection.selectHit(id);if(!hit)return;
    alignmentIndex=span?Math.max(0,(hit.alignments||[]).findIndex(a=>a.query_start<=span[1]&&a.query_end>=span[0])):0;
    updateHighlight();
    if(fromGraph)reveal(el('hits-panel'),el('domain-rows').querySelector(`[data-hit="${id}"]`),32);
    else reveal(el('architecture-scroll'),el('domain-architecture').querySelector(`[data-hit="${id}"]`),6);
  }
  function describe(hit,span=null){
    return `${hit.source} · ${hit.name} · ${hit.accession}\n${hit.description.length>360?hit.description.slice(0,357)+'…':hit.description}\nResidues ${span?span.join('–'):hit.spans.map(s=>s.join('–')).join(', ')} · E-value ${number(hit.evalue)} · bitscore ${hit.bitscore.toFixed(1)}`;
  }
  function updateHighlight(){
    const hit=currentHit();
    for(const item of document.querySelectorAll('[data-hit]')){
      const selected=item.dataset.hit===hit?.id;item.classList.toggle('selected',selected);item.setAttribute('aria-selected',String(selected));
    }
    el('copy-hit').disabled=el('inspect-hit').disabled=!hit;
    el('hit-detail').textContent=hit?`${hit.source} · ${hit.name} · ${hit.accession} · ${hit.spans.map(s=>s.join('–')).join(', ')}`:'Select a domain or row to inspect its evidence.';
    el('hit-detail').title=hit?hit.description:'';
    if(detailsOpen)renderDetails();
  }
  function renderHits(){
    hideTooltip();layout=model.architecture(hits(),selection.query.length);
    renderTable();renderArchitecture();updateHighlight();
  }
  function renderTable(){
    const visible=hits().slice().sort((a,b)=>{
      const delta=typeof a[sort]==='number'?a[sort]-b[sort]:a[sort].localeCompare(b[sort]);
      return (ascending?delta:-delta)||a.evalue-b.evalue||a.id.localeCompare(b.id);
    });
    el('hit-count').textContent=`${visible.length} matches`;
    el('hit-count').title=['Pfam','CDD'].map(name=>`${visible.filter(h=>h.source===name).length} ${name}`).join(' · ')+(weaker?'':' · score ≥6');
    el('domain-empty').hidden=visible.length>0;
    el('domain-empty').textContent=selection.query.hits.length?'No domains match these filters. Try All or Show weaker hits.':'No significant domains detected. This does not imply that the protein lacks a function.';
    el('domain-rows').replaceChildren(...visible.map(hit=>{
      const row=document.createElement('tr');row.dataset.hit=hit.id;row.tabIndex=0;
      row.setAttribute('aria-label',`${hit.source} ${hit.name}, residues ${hit.start} to ${hit.end}`);
      for(const value of [hit.source,hit.name,hit.accession,`${hit.start}–${hit.end}`,number(hit.evalue),hit.bitscore.toFixed(1),hit.description]){
        const cell=document.createElement('td');cell.textContent=String(value);cell.title=String(value);row.append(cell);
      }
      row.onclick=()=>activate(hit.id);
      row.onkeydown=event=>{
        if(event.key==='Enter'||event.key===' '){event.preventDefault();activate(hit.id);}
        else if(event.key==='ArrowDown'||event.key==='ArrowUp'){
          const next=event.key==='ArrowDown'?row.nextElementSibling:row.previousElementSibling;
          if(next){event.preventDefault();next.focus({preventScroll:true});activate(next.dataset.hit);reveal(el('hits-panel'),next,32);}
        }
      };
      return row;
    }));
    for(const button of document.querySelectorAll('[data-sort]'))button.parentElement.setAttribute('aria-sort',button.dataset.sort===sort?(ascending?'ascending':'descending'):'none');
  }
  const svgNS='http://www.w3.org/2000/svg';
  function svg(tag,attributes={},text=null){
    const node=document.createElementNS(svgNS,tag);
    for(const [key,value] of Object.entries(attributes))node.setAttribute(key,String(value));
    if(text!==null)node.textContent=text;return node;
  }
  function hideTooltip(){el('domain-tooltip').hidden=true;}
  function renderArchitecture(){
    if(!selection.query||!layout||el('quick-workspace').hidden)return;
    hideTooltip();
    const container=el('domain-architecture'),width=container.clientWidth;
    if(!width)return;
    const left=48,right=width-16,span=right-left,height=Math.max(12,layout.tracks*28+8);
    container.setAttribute('viewBox',`0 0 ${width} ${height}`);container.setAttribute('height',height);
    const scale=el('domain-scale');scale.setAttribute('viewBox',`0 0 ${width} 46`);scale.style.width=width+'px';
    const ticks=[svg('line',{x1:left,y1:14,x2:right,y2:14,class:'backbone'})],length=selection.query.length;
    const tickCount=Math.min(6,Math.max(2,Math.floor(span/80)),length);
    for(let i=0;i<tickCount;i++){
      const residue=Math.round(1+(length-1)*i/Math.max(1,tickCount-1)),x=left+(residue-1)/length*span;
      ticks.push(svg('line',{x1:x,y1:9,x2:x,y2:20,class:'tick'}),svg('text',{x,y:35,'text-anchor':i===0?'start':i===tickCount-1?'end':'middle',class:'coordinate'},residue));
    }
    scale.replaceChildren(...ticks);
    const nodes=[];
    for(const group of layout.groups)nodes.push(svg('text',{x:3,y:18+group.track*28,class:'source-label'},group.source));
    for(const block of layout.blocks){
      const x=left+block.x*span,w=block.width*span,y=2+block.track*28;
      const node=svg('g',{'data-hit':block.id,tabindex:0,role:'button','aria-label':`${block.hit.source} ${block.hit.name}, residues ${block.start} to ${block.end}`,class:'domain-block'});
      node.append(svg('rect',{x,y,width:w,height:22,rx:4,fill:block.color,class:block.hit.source==='Pfam'?'pfam-block':'cdd-block'}));
      const chars=Math.floor((w-10)/6.5);
      if(chars>=4){const name=block.hit.name;node.append(svg('text',{x:x+w/2,y:y+15,'text-anchor':'middle',class:'domain-label'},name.length>chars?name.slice(0,chars-1)+'…':name));}
      node.onclick=()=>activate(block.id,true,[block.start,block.end]);
      node.onkeydown=event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();activate(block.id,true,[block.start,block.end]);}};
      const tooltip=event=>{
        const tip=el('domain-tooltip');tip.textContent=describe(block.hit,[block.start,block.end]);tip.hidden=false;
        const box=node.getBoundingClientRect();tip.style.left=Math.max(8,Math.min(event.clientX||box.x,innerWidth-tip.offsetWidth-12))+'px';
        tip.style.top=Math.max(8,Math.min((event.clientY||box.y)+18,innerHeight-tip.offsetHeight-12))+'px';
      };
      node.onmouseenter=tooltip;node.onfocus=tooltip;node.onmouseleave=node.onblur=hideTooltip;nodes.push(node);
    }
    container.replaceChildren(...nodes);
    el('track-count').textContent=layout.tracks>1?`${layout.tracks} tracks · scroll for overlapping hits`:'1-based, inclusive residues';
    updateHighlight();
  }
  let lastWidth=0,resizeFrame;
  new ResizeObserver(entries=>{
    const width=entries[0].contentRect.width;if(Math.abs(width-lastWidth)<1)return;
    lastWidth=width;cancelAnimationFrame(resizeFrame);resizeFrame=requestAnimationFrame(()=>{renderArchitecture();detailKey='';if(detailsOpen)renderDetails();});
  }).observe(el('architecture-scroll'));
  el('architecture-scroll').addEventListener('scroll',hideTooltip,{passive:true});
  for(const button of el('source-filters').children)button.onclick=()=>{
    saveView();source=button.dataset.source;
    for(const item of el('source-filters').children){const selected=item===button;item.classList.toggle('active',selected);item.setAttribute('aria-pressed',String(selected));}
    if(selection.query){renderHits();restoreView();}
  };
  el('show-weaker').onchange=()=>{saveView();weaker=el('show-weaker').checked;if(selection.query){renderHits();restoreView();}};
  for(const button of document.querySelectorAll('[data-sort]'))button.onclick=()=>{
    ascending=sort===button.dataset.sort?!ascending:true;sort=button.dataset.sort;
    if(selection.query){renderTable();el('hits-panel').scrollTop=0;updateHighlight();}
  };
  function showDetails(open){
    detailsOpen=open;el('details-panel').hidden=!open;el('hits-panel').hidden=open;
    for(const [id,active] of [['show-hits',!open],['show-details',open]]){el(id).classList.toggle('active',active);el(id).setAttribute('aria-selected',String(active));}
    if(open)renderDetails();
  }
  el('show-hits').onclick=()=>showDetails(false);el('show-details').onclick=el('inspect-hit').onclick=()=>showDetails(true);
  function renderDetails(){
    const hit=currentHit(),key=`${hit?.id}:${alignmentIndex}`;
    if(key===detailKey)return;detailKey=key;el('details-panel').scrollTop=0;
    el('detail-empty').hidden=!!hit;el('selected-evidence').hidden=!hit;if(!hit)return;
    el('evidence-name').textContent=`${hit.source} · ${hit.name} · ${hit.accession}`;
    el('evidence-stats').textContent=`Residues ${hit.spans.map(s=>s.join('–')).join(', ')} · E-value ${number(hit.evalue)} · bitscore ${hit.bitscore.toFixed(1)}`;
    el('evidence-description-text').textContent=hit.description;
    const alignments=hit.alignments||[];alignmentIndex=Math.min(alignmentIndex,Math.max(0,alignments.length-1));
    el('alignment-segment').hidden=alignments.length<2;
    el('alignment-segment').replaceChildren(...alignments.map((a,i)=>{
      const option=document.createElement('option');option.value=i;option.textContent=`Segment ${i+1}: query ${a.query_start}–${a.query_end}`;option.selected=i===alignmentIndex;return option;
    }));
    const alignment=alignments[alignmentIndex];el('copy-alignment').hidden=!alignment;
    if(alignment)el('evidence-stats').textContent=`Query ${alignment.query_start}–${alignment.query_end} · model ${alignment.model_start}–${alignment.model_end} · E-value ${number(alignment.evalue)} · bitscore ${alignment.bitscore.toFixed(1)}`;
    el('alignment-title').textContent=alignment?'Query ↔ CDD model consensus':'Alignment';
    el('alignment-note').textContent=alignment?'Statistics above describe this aligned segment (HSP). | identical residue; – gap. The model is a profile consensus, not an individual protein.':
      hit.source==='Pfam'?'Pfam alignment strings are not retained by the current HMMER search. Coordinates and statistics remain available above.':'No alignment strings are available for this hit.';
    el('alignment-text').hidden=!alignment;
    if(alignment){const columns=Math.max(20,Math.min(70,Math.floor((el('details-panel').clientWidth-155)/7.3/10)*10));el('alignment-text').textContent=model.alignmentText(alignment,columns);}
  }
  el('alignment-segment').onchange=()=>{alignmentIndex=Number(el('alignment-segment').value);renderDetails();};
  async function copy(text){
    try{
      if(navigator.clipboard?.writeText)await navigator.clipboard.writeText(text);
      else{const buffer=el('copy-buffer'),focused=document.activeElement;buffer.value=text;buffer.focus();buffer.select();if(!document.execCommand('copy'))throw Error('Select the text and press ⌘C to copy.');focused?.focus({preventScroll:true});}
      el('copy-status').textContent='Copied';
    }catch(error){el('copy-status').textContent=error.message;}
  }
  function copyHit(){const h=currentHit();if(h)copy([h.source,h.name,h.accession,h.start,h.end,number(h.evalue),h.bitscore,h.description].join('\t'));}
  el('copy-hit').onclick=copyHit;
  el('copy-alignment').onclick=()=>{const hit=currentHit(),a=hit?.alignments?.[alignmentIndex];if(a)copy(`${selection.query.name}\nCDD ${hit.name} ${hit.accession}\n${model.alignmentText(a)}`);};
  el('copy-query').onclick=()=>{const q=selection.query;if(q)copy(`>${q.name}\n${q.sequence}\n`);};
  document.addEventListener('keydown',event=>{
    if((event.metaKey||event.ctrlKey)&&event.key.toLowerCase()==='c'&&document.activeElement?.closest('[data-hit]')&&!window.getSelection().toString()){
      event.preventDefault();copyHit();
    }
  });
})();
