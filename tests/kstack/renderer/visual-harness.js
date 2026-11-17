// Run in the existing browser/WKWebView harness against public or synthetic data.
// Returns artifacts to the caller; performs no filesystem or network operations.
(async function () {
  const C=window.Kstack, M=window.CSModel, $=id=>document.getElementById(id);
  const original=C.getProject(), wasNative=document.body.classList.contains('native');
  const checks=[], exports=[], timings=[], textMetrics=new Map();
  const frame=()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)));
  const check=(condition,name)=>{if(!condition)throw Error(name);checks.push(name);};
  const equal=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
  const setScale=(x,y,spacing=100)=>{for(const [axis,value]of [['width',x],['height',y],['spacing',spacing]]){const input=$(`${axis}-scale`);input.value=value;input.dispatchEvent(new Event('input'));}};
  const content=p=>JSON.stringify([p.clusters,p.layout,p.homolog_groups,p.pairwise_hits]);
  const previewMarkup=svg=>{
    const copy=svg.cloneNode(true);
    for(const gene of copy.querySelectorAll('path[data-gene]'))for(const attr of ['tabindex','role','aria-haspopup','aria-label'])gene.removeAttribute(attr);
    return copy.innerHTML;
  };
  const near=(a,b)=>Math.abs(a-b)<.001;
  const points=path=>[...path.matchAll(/[ML]([\d.e+-]+),([\d.e+-]+)/g)].map(m=>[+m[1],+m[2]]);
  function validateText(root,source,name) {
    check(!root.querySelector('[transform],[textLength],[lengthAdjust]'),`${name}: ${source} has no glyph-stretching attributes`);
    const sy=M.geometry(C.getProject()).sy;
    for(const [i,node]of [...root.querySelectorAll('text')].entries()) {
      const style=getComputedStyle(node),size=parseFloat(style.fontSize),width=node.getComputedTextLength();
      check(style.transform==='none',`${name}: ${source} text ${i} has no CSS transform`);
      const matrix=node.getScreenCTM();
      check(near(matrix.a,matrix.d)&&near(matrix.b,0)&&near(matrix.c,0)&&near(size,+node.getAttribute('font-size')),`${name}: ${source} text ${i} has a uniform screen projection and native font size`);
      // Native optical font sizing can change letter widths at different font
      // sizes. Compare Width/spacing changes at the same Height/font size.
      const key=`${source}/${i}/${node.textContent}/${sy}`,reference=textMetrics.get(key),metric={size,width};
      if(reference)check(near(metric.size,reference.size)&&Math.abs(metric.width-reference.width)<.15,`${name}: ${source} text ${i} retains natural font proportions (${metric.width} vs ${reference.width})`);
      else textMetrics.set(key,metric);
    }
  }
  function validateArrow(path,s,name) {
    const v=points(path),width=s.x2-s.x1;
    check(near(Math.min(...v.map(p=>p[0])),s.x1)&&near(Math.max(...v.map(p=>p[0])),s.x2),`${name}: arrow retains genomic span`);
    check(near(Math.max(...v.map(p=>p[1]))-Math.min(...v.map(p=>p[1])),20*s.sy),`${name}: arrow height follows Height only`);
    if(s.strand) {
      check(v.length===(width<=10*s.sy?3:7),`${name}: short-arrow clamp`);
      if(v.length===7)check(near(Math.abs(v[3][0]-v[2][0])/(v[4][1]-v[2][1]),.5)&&near(v[0][1],v[1][1])&&near(v[5][1],v[6][1]),`${name}: stable head ratio and rectangular shaft`);
      const mirror=points(M.arrowPath({...s,strand:-s.strand}));
      check(v.every(([x,y],i)=>near(x+mirror[i][0],s.x1+s.x2)&&near(y,mirror[i][1])),`${name}: exact strand mirror`);
    }
  }
  function validateExport(svg,name) {
    const doc=new DOMParser().parseFromString(svg,'image/svg+xml');
    check(!doc.querySelector('parsererror'),`${name}: valid SVG XML`);
    const root=document.importNode(doc.documentElement,true);
    // Test-only temporary measurement of actual text/primitive bounds, never a product export tree.
    root.style.position='absolute';root.style.visibility='hidden';root.style.pointerEvents='none';
    document.body.append(root);
    try {
      const [left,top,width,height]=root.getAttribute('viewBox').split(' ').map(Number);
      check(+root.getAttribute('width')===width&&+root.getAttribute('height')===height,`${name}: dimensions match viewBox`);
      const outside=[];
      for(const node of root.querySelectorAll('text,path,polygon,line')) {
        const b=node.getBBox();
        for(const [x,y]of [[b.x,b.y],[b.x+b.width,b.y+b.height]]) {
          if(x<left-.001||x>left+width+.001||y<top-.001||y>top+height+.001)outside.push(node.outerHTML.slice(0,120));
        }
      }
      check(!outside.length,`${name}: all vector/text bounds inside export (${outside.join(';')})`);
      check(!root.querySelector('image,foreignObject,use,canvas,g'),`${name}: editable root SVG primitives`);
      validateText(root,'export',name);
      const p=C.getProject(),geo=M.geometry(p);
      for(const path of root.querySelectorAll('[data-kind=gene]')) {
        const s=geo.genes.get(path.dataset.gene).segments[+path.dataset.segment];
        check(path.getAttribute('d')===M.arrowPath(s),`${name}: exported gene ${path.dataset.gene}/${path.dataset.segment}`);
        validateArrow(path.getAttribute('d'),s,`${name}: export ${path.dataset.gene}/${path.dataset.segment}`);
      }
    } finally { root.remove(); }
  }
  try {
    C.setNativeMode(true);await frame();
    // Both UI entry points share disposable export state. Exercise the actual
    // controls and pointer handlers in WKWebView, with capture stubbed only for
    // synthetic pointer events (real input is qualified separately).
    for (const native of [false, true]) {
      C.setNativeMode(native); await frame();
      const saved=C.getProject(),source=C.exportSVG(),options={order:[...saved.layout.order],hidden:[],offsets:{},visualization:{...saved.visualization}};
      const prefix=native?'Native export preview':'Offline export preview';
      const previewSVG=()=>$('export-preview-canvas').querySelector('svg');
      const compare=()=>{
        const expected=document.importNode(new DOMParser().parseFromString(window.CSSVG.generate(saved,options),'image/svg+xml').documentElement,true);
        check(previewMarkup(previewSVG())===expected.innerHTML&&previewSVG().getAttribute('viewBox')===expected.getAttribute('viewBox'),`${prefix}: preview matches model-generated export`);
        check(equal(C.getProject(),saved)&&C.exportSVG()===source,`${prefix}: main stack and dirty timestamp unchanged`);
      };
      $('export-svg').focus();$('export-svg').click(); await frame();
      check($('export-preview').open,`${prefix}: toolbar opens dialog`); compare();
      if (saved.clusters.length>1) {
        const hidden=options.order[1];options.hidden=[hidden];
        $('export-clusters').querySelector(`[data-cluster="${hidden}"] input`).click();await frame();compare();
        check(previewSVG().querySelectorAll('[data-kind="baseline"]').length===saved.clusters.length-1,`${prefix}: hidden row leaves no gap`);
        const id=options.order[0],handle=$('export-clusters').querySelector(`[data-cluster="${id}"] .export-handle`);
        handle.dispatchEvent(new KeyboardEvent('keydown',{key:'ArrowDown',bubbles:true}));
        options.order.splice(0,1);options.order.splice(1,0,id);await frame();compare();
        check(document.activeElement===handle,`${prefix}: keyboard reorder retains focus`);
      }
      for(const [id,key,value]of [['width','width_scale',65],['height','height_scale',175],['spacing','row_spacing_scale',150],['identity','connector_identity',100]]) {
        $(`export-${id}`).value=value;$(`export-${id}`).dispatchEvent(new Event('input'));options.visualization[key]=value/100;
      }
      $('export-labels').click();options.visualization.labels=!saved.visualization.labels;await frame();compare();
      $('export-connectors').click();options.visualization.connectors=!saved.visualization.connectors;await frame();compare();
      $('export-hide-all').click();await frame();
      check($('export-save').disabled&&!previewSVG()&&$('export-preview-canvas').textContent.includes('at least one'),`${prefix}: empty selection disables export`);
      $('export-show-all').click();options.hidden=[];await frame();compare();
      $('export-reset').click();options.order=[...saved.layout.order];options.visualization={...saved.visualization};await frame();compare();
      if (saved.clusters.length>1) {
        const list=$('export-clusters'),capture=list.setPointerCapture;list.setPointerCapture=()=>{};
        try {
          const first=list.children[0],second=list.children[1],a=first.getBoundingClientRect(),b=second.getBoundingClientRect();
          first.querySelector('.export-handle').dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:81,clientY:a.top+a.height/2,bubbles:true}));
          window.dispatchEvent(new PointerEvent('pointermove',{pointerId:81,clientY:b.top+b.height/2,bubbles:true}));
          window.dispatchEvent(new PointerEvent('pointerup',{pointerId:81,bubbles:true}));
          [options.order[0],options.order[1]]=[options.order[1],options.order[0]];await frame();compare();
          check(list.children[1].dataset.cluster===options.order[1],`${prefix}: list pointer drag reorders`);
        } finally {list.setPointerCapture=capture;}
        const canvas=$('export-preview-canvas'),captureCanvas=canvas.setPointerCapture;canvas.setPointerCapture=()=>{};
        try {
          const geo=M.geometry(window.CSSVG.projectView(saved,options)),matrix=previewSVG().getScreenCTM();
          const a=new DOMPoint(geo.left,geo.rows.get(options.order[0]).y).matrixTransform(matrix);
          const b=new DOMPoint(geo.left,geo.rows.get(options.order[1]).y).matrixTransform(matrix);
          previewSVG().dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:82,clientX:a.x,clientY:a.y,bubbles:true}));
          window.dispatchEvent(new PointerEvent('pointermove',{pointerId:82,clientX:b.x,clientY:b.y,bubbles:true}));
          [options.order[0],options.order[1]]=[options.order[1],options.order[0]];await frame();compare();
          window.dispatchEvent(new PointerEvent('pointercancel',{pointerId:82,bubbles:true}));
          [options.order[0],options.order[1]]=[options.order[1],options.order[0]];await frame();compare();
          check(list.children[0].dataset.cluster===options.order[0],`${prefix}: figure drag updates live and cancellation restores order`);
          // Horizontal quick drag: pixel motion must stay stable while the
          // normalized minimum and fitted SVG extent change under the pointer.
          for(const fit of [true,false]) {
            $('export-fit').checked=fit;$('export-fit').dispatchEvent(new Event('change'));await frame();
            const view=window.CSSVG.projectView(saved,options),geo=M.geometry(view);
            const cid=view.layout.order.reduce((a,b)=>view.layout.rows[a].offset<view.layout.rows[b].offset?a:b);
            const other=view.layout.order.find(id=>id!==cid),otherGene=geo.index.clusters.get(other).cds[0].id;
            const screen=id=>canvas.querySelector(`[data-gene="${id}"]`).getBoundingClientRect().left;
            const stationary=screen(otherGene),matrix=previewSVG().getScreenCTM();
            const a=new DOMPoint(geo.rows.get(cid).x,geo.rows.get(cid).y).matrixTransform(matrix);
            const delta=-40;
            previewSVG().dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:83,clientX:a.x,clientY:a.y,bubbles:true}));
            window.dispatchEvent(new PointerEvent('pointermove',{pointerId:83,clientX:a.x+delta,clientY:a.y+1,bubbles:true}));await frame();
            check(near(screen(otherGene),stationary),`${prefix}: horizontal drag keeps other rows stationary (fit ${fit})`);
            options.offsets[cid]=Math.round((view.layout.rows[cid].offset+delta/(matrix.a*geo.scale))*1e6)/1e6;
            window.dispatchEvent(new PointerEvent('pointerup',{pointerId:83,bubbles:true}));await frame();compare();
            const moved=canvas.querySelector(`[data-gene="${geo.index.clusters.get(cid).cds[0].id}"]`);
            moved.dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:84,clientX:a.x,clientY:a.y,bubbles:true}));
            window.dispatchEvent(new PointerEvent('pointermove',{pointerId:84,clientX:a.x+50,clientY:a.y,bubbles:true}));await frame();
            window.dispatchEvent(new PointerEvent('pointercancel',{pointerId:84,bubbles:true}));await frame();compare();
          }
        } finally {canvas.setPointerCapture=captureCanvas;}
        const reference=saved.clusters.find(c=>c.id===options.order[0]).cds.find(g=>g.segments.length);
        const path=canvas.querySelector(`[data-gene="${reference.id}"]`),rect=path.getBoundingClientRect();
        path.dispatchEvent(new MouseEvent('contextmenu',{bubbles:true,cancelable:true,clientX:rect.left,clientY:rect.top}));
        const menu=$('export-preview').querySelector('.gene-menu');
        check(!menu.hidden&&document.activeElement===menu.firstChild,`${prefix}: gene menu opens and focuses action`);
        Object.assign(options.offsets,M.alignmentOffsets(window.CSSVG.projectView(saved,options),reference.id));
        menu.firstChild.click();await frame();compare();
        check(menu.hidden,`${prefix}: auto-align closes menu`);
        const focused=canvas.querySelector(`[data-gene="${reference.id}"]`);focused.focus();
        focused.dispatchEvent(new KeyboardEvent('keydown',{key:'F10',shiftKey:true,bubbles:true,cancelable:true}));
        check(!menu.hidden,`${prefix}: keyboard gene actions available`);
        menu.firstChild.dispatchEvent(new KeyboardEvent('keydown',{key:'Escape',bubbles:true,cancelable:true}));
        check(menu.hidden&&document.activeElement===focused&&$('export-preview').open,`${prefix}: Escape restores gene focus without closing preview`);
        $('export-reset').click();options.offsets={};options.order=[...saved.layout.order];await frame();compare();
      }
      $('export-cancel').click();await frame();
      check(!$('export-preview').open&&document.activeElement===$('export-svg'),`${prefix}: cancel restores focus`);
      C.openExportPreview();await frame();
      check($('export-preview').open,`${prefix}: native menu API opens preview`);
      const expected=document.importNode(new DOMParser().parseFromString(source,'image/svg+xml').documentElement,true);
      check(previewMarkup(previewSVG())===expected.innerHTML,`${prefix}: reopening starts from the main stack`);
      $('export-preview').dispatchEvent(new Event('cancel',{cancelable:true}));await frame();
      check(!$('export-preview').open&&equal(C.getProject(),saved),`${prefix}: Escape closes without changes`);
    }
    C.setNativeMode(true);await frame();
    {
    const baseline=C.getProject(),moved=structuredClone(baseline),reference=moved.clusters[0].cds.find(g=>g.segments.length);
    for(const cid of moved.layout.order)moved.layout.rows[cid].offset+=moved.layout.order.indexOf(cid)*431;
    C.setProject(moved);await frame();
    const expectedOffsets=M.alignmentOffsets(moved,reference.id),gene=$('genes-svg').querySelector(`[data-gene="${reference.id}"]`);
    gene.dispatchEvent(new MouseEvent('contextmenu',{bubbles:true,cancelable:true,clientX:100,clientY:100}));
    const menu=document.body.querySelector(':scope > .gene-menu');
    check(!menu.hidden,'Main: gene context menu opens');menu.firstChild.click();await frame();
    const aligned=C.getProject();
    for(const cid of moved.layout.order)check(aligned.layout.rows[cid].offset===(expectedOffsets[cid]??moved.layout.rows[cid].offset),'Main: auto-align updates only calculated offsets');
    check(equal(aligned.clusters,moved.clusters)&&equal(aligned.pairwise_hits,moved.pairwise_hits)&&equal(aligned.homolog_groups,moved.homolog_groups),'Main: auto-align preserves biological data');
    check(!$('genes-svg').querySelector('.cluster-info'),'Main: cluster size/direction subtext removed');
    C.setProject(baseline);await frame();
    }
    check(!document.querySelector('canvas'), 'No bitmap connector intermediate');
    $('labels').checked=true;$('labels').dispatchEvent(new Event('change'));await frame();
    const cases=[[100,100,100],[50,100,100],[200,100,100],[100,50,100],[100,200,100],[65,175,100],
      [100,100,50],[100,100,150],[65,175,150],[25,400,100],[400,25,100],[25,400,25],[400,25,400]];
    for(const [x,y,spacing]of cases) {
      const start=performance.now();setScale(x,y,spacing);
      // Hold the camera constant for font measurements. Native optical sizing
      // may choose different letterforms at a different displayed zoom level.
      C.zoom(1/C.getViewport().zoom);await frame();
      const name=`${x}x${y}-spacing${spacing}`,p=C.getProject(),geo=M.geometry(p);
      check(content(p)===content(original),`${name}: biology and alignment unchanged`);
      check(getComputedStyle($('scene')).transform==='none'&&getComputedStyle($('scene')).willChange==='auto',`${name}: scene is not CSS magnified/promoted`);
      validateText($('genes-svg'),'live',name);
      const paths=[...document.querySelectorAll('.gene')];
      check(paths.length===[...geo.genes.values()].reduce((n,g)=>n+g.segments.length,0),`${name}: all gene segments rendered`);
      for(const path of paths) {
        const s=geo.genes.get(path.dataset.gene).segments.find(s=>M.arrowPath(s)===path.getAttribute('d'));
        check(!!s,`${name}: live gene ${path.dataset.gene}`);
        validateArrow(path.getAttribute('d'),s,`${name}: live ${path.dataset.gene}`);
      }
      const connectors=M.connectors(p,geo),nodes=[...$('connectors-svg').children];
      check(nodes.length===connectors.length&&nodes.every((n,i)=>n.getAttribute('points')===connectors[i].points.map(p=>p.join(',')).join(' ')),`${name}: connector alignment`);
      for(const c of connectors)check(near(c.points[0][1],geo.genes.get(c.query).row.y+11*geo.sy)&&near(c.points[2][1],geo.genes.get(c.target).row.y-11*geo.sy),`${name}: connector endpoints follow row centers`);
      const z=C.getViewport().zoom, viewportBox=$('viewport').getBoundingClientRect();
      for(const [i,cid]of p.layout.order.entries()) {
        const actions=document.querySelector(`#row-controls [data-cluster="${cid}"]`).getBoundingClientRect();
        const row=geo.rows.get(cid),expected=viewportBox.top+$('viewport').clientTop+C.getViewport().y+(row.y+30*geo.sy)*z;
        check(Math.abs(actions.top-expected)<1,`${name}: row controls follow their row anchor`);
        const zone=document.querySelector(`#gene-rows [data-cluster="${cid}"] .row-hitzone`);
        check(+zone.getAttribute('height')<=geo.rowHeight&&near(+zone.getAttribute('y')+ +zone.getAttribute('height')/2,row.y),`${name}: alignment hit target cannot cover adjacent genes`);
        if(i)check(near(row.y-geo.rows.get(p.layout.order[i-1]).y,90*y/100*spacing/100),`${name}: independent row-center distance`);
      }
      const svg=C.exportSVG();validateExport(svg,name);exports.push({name,svg});
      const prior=C.getProject().visualization;C.fit();await frame();
      check(equal(C.getProject().visualization,prior),`${name}: Fit preserves presentation scales`);
      check(geo.width*C.getViewport().zoom<=$('viewport').clientWidth-36,`${name}: Fit contains complete stack width`);
      C.zoom(12/C.getViewport().zoom);await frame();
      check(paths[0]===document.querySelector('.gene')&&nodes[0]===$('connectors-svg').firstElementChild,`${name}: 1200% zoom reuses vector nodes`);
      check(C.exportSVG()===svg,`${name}: camera does not alter export`);
      timings.push({name,renderMs:C.metrics.lastEditMs,connectorMs:C.metrics.connectorMs,validationMs:performance.now()-start});
    }
    setScale(100,100,25);C.fit();await frame();
    const figureBox=$('viewport').getBoundingClientRect();
    const visibleGene=[...document.querySelectorAll('.gene')].find(path=>{
      const b=path.getBoundingClientRect();return b.width>5&&b.left>figureBox.left&&b.right<figureBox.right&&b.top>figureBox.top&&b.bottom<figureBox.bottom;
    });
    check(!!visibleGene,'Dense spacing has a visible gene hit target');
    const geneBox=visibleGene.getBoundingClientRect(),hit=document.elementFromPoint(geneBox.x+geneBox.width/2,geneBox.y+geneBox.height/2);
    check(hit?.closest('.gene')?.dataset.gene===visibleGene.dataset.gene,'Dense spacing leaves gene hit testing unobstructed');
    hit.dispatchEvent(new MouseEvent('click',{bubbles:true}));await frame();
    check(document.querySelector('.gene.selected')?.dataset.gene===visibleGene.dataset.gene,'Gene click opens inspection at dense spacing');
    setScale(65,175,150);C.fit();await frame();
    const p=C.getProject(),cid=p.layout.order[1]||p.layout.order[0],gene=p.clusters.find(c=>c.id===cid).cds.find(g=>g.protein&&g.segments.length);
    C.inspect(gene.id);await frame();
    check(document.querySelector('.sequence').textContent===gene.protein,'Inspector shows the complete protein');
    const color=document.querySelector('[aria-label="Gene override color"]');color.value='#123456';
    document.querySelector('[aria-label="Apply color to this gene"]').click();await frame();
    check(C.getProject().color_overrides.genes[gene.id]==='#123456','Gene color interaction');
    document.querySelector('[aria-label="Apply color to homolog group"]').click();await frame();
    check(C.getProject().color_overrides.groups[gene.group_id]==='#123456','Group color interaction');
    document.querySelector('[aria-label="Reset group to generated colors"]').click();await frame();
    check(!C.getProject().color_overrides.groups[gene.group_id],'Group color reset');
    C.selectGene(gene.id);C.selectCluster(cid);await frame();
    check(C.exportFASTA().includes(gene.protein.slice(0,60)),'Gene/cluster FASTA selection');
    C.edit('reverse',cid);await frame();check(C.getProject().layout.rows[cid].orientation===-p.layout.rows[cid].orientation,'Cluster reversal');
    // Dispatch the same pointer handlers with an actual scaled coordinate delta.
    const viewport=$('viewport'),capture=viewport.setPointerCapture;viewport.setPointerCapture=()=>{};
    try {
      const z=C.getViewport().zoom,geometry=M.geometry(C.getProject()),before=C.getProject().layout.rows[cid].offset;
      document.querySelector(`[data-cluster="${cid}"] .row-hitzone`).dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:42,clientX:300,bubbles:true}));
      window.dispatchEvent(new PointerEvent('pointermove',{pointerId:42,clientX:300+1000*geometry.scale*z,bubbles:true}));
      window.dispatchEvent(new PointerEvent('pointerup',{pointerId:42,bubbles:true}));await frame();
      check(Math.abs(C.getProject().layout.rows[cid].offset-before-1000)<.002,'Horizontal drag maps scaled pixels back to bp');
      if(p.clusters.length>1) {
        const beforeIndex=C.getProject().layout.order.indexOf(cid),delta=beforeIndex? -1:1;
        document.querySelector(`[data-cluster="${cid}"] .handle`).dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:43,clientY:300,bubbles:true}));
        window.dispatchEvent(new PointerEvent('pointermove',{pointerId:43,clientY:300+delta*geometry.rowHeight*z,bubbles:true}));
        window.dispatchEvent(new PointerEvent('pointerup',{pointerId:43,bubbles:true}));await frame();
        check(C.getProject().layout.order.indexOf(cid)===beforeIndex+delta,'Reorder drag uses scaled row spacing');
      }
      const saved=C.getProject(),camera=C.getViewport(),svg=C.exportSVG();
      viewport.dispatchEvent(new PointerEvent('pointerdown',{button:0,pointerId:44,clientX:300,clientY:300,bubbles:true}));
      window.dispatchEvent(new PointerEvent('pointermove',{pointerId:44,clientX:345,clientY:330,bubbles:true}));
      window.dispatchEvent(new PointerEvent('pointerup',{pointerId:44,bubbles:true}));await frame();
      check(near(C.getViewport().x-camera.x,45)&&near(C.getViewport().y-camera.y,30)&&equal(C.getProject(),saved)&&C.exportSVG()===svg,'Pan preserves saved presentation, dirty timestamp, and export');
    } finally { viewport.setPointerCapture=capture; }
    const groups=C.getProject().homolog_groups,identity=$('display-identity');
    identity.value=100;identity.dispatchEvent(new Event('change'));await frame();
    check($('connectors-svg').children.length===M.connectors(C.getProject()).length&&equal(C.getProject().homolog_groups,groups),'Identity filtering preserves homolog groups');
    identity.value=original.visualization.connector_identity*100;identity.dispatchEvent(new Event('change'));
    $('labels').click();$('connectors').click();await frame();
    check(!document.querySelector('#gene-rows text[text-anchor="middle"]'),'Gene label toggle');
    check($('connectors-svg').children.length===0,'Connector toggle');
    $('connectors').click();await frame();
    $('reset-scale').click();await frame();
    check(['width_scale','height_scale','row_spacing_scale'].every(key=>C.getProject().visualization[key]===1),'Reset restores 100% / 100% / 100%');
    $('width-scale').value=900;$('width-scale').dispatchEvent(new Event('change'));await frame();
    check(C.getProject().visualization.width_scale===4&&$('width-scale').value==='400','Out-of-range scale input clamps on change');
    $('height-scale').value='';$('height-scale').dispatchEvent(new Event('change'));
    check($('height-scale').value==='100','Empty scale input restores the actual value');
    $('spacing-scale').value=5;$('spacing-scale').dispatchEvent(new Event('change'));await frame();
    check(C.getProject().visualization.row_spacing_scale===.25&&$('spacing-scale').value==='25','Row spacing validates with the same bounds');
    return {ok:true,checks,exports,timings};
  } finally {C.setProject(original);C.setNativeMode(wasNative);await frame();}
})
