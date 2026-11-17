const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('../../../src/kstack/resources/renderer/model.js');
const SVG = require('../../../src/kstack/resources/export/svg.js');
const fixture = require('../fixtures/synthetic/project.json');
const points = path => [...path.matchAll(/[ML]([\d.e+-]+),([\d.e+-]+)/g)].map(m => [+m[1],+m[2]]);
const close = (a,b) => assert.ok(Math.abs(a-b)<1e-8,`${a} ≠ ${b}`);

test('forward block arrow has a flat shaft, vertical tail and fixed triangular head', () => {
  assert.deepEqual(points(M.arrowPath({x1:20,x2:80,y:50,strand:1})),
    [[20,43],[70,43],[70,40],[80,50],[70,60],[70,57],[20,57]]);
  for (const width of [11,20,75,300]) {
    const shape=points(M.arrowPath({x1:0,x2:width,y:0,strand:1}));
    close(shape[3][0]-shape[2][0],10);
  }
});
test('reverse arrows exactly mirror the forward vertices', () => {
  for (const width of [.025,5,10,10.01,75]) {
    const s={x1:20,x2:20+width,y:50,strand:1}, forward=points(M.arrowPath(s));
    const reverse=points(M.arrowPath({...s,strand:-1}));
    reverse.forEach(([x,y],i)=>{close(x,s.x1+s.x2-forward[i][0]);close(y,forward[i][1]);});
  }
});
test('short genes clamp to nonintersecting triangles without inflating genomic width', () => {
  for (const width of [.025,.5,2,9.999,10]) for(const strand of [-1,1]) {
    const shape=points(M.arrowPath({x1:2,x2:2+width,y:30,strand}));
    assert.equal(shape.length,3);assert.ok(shape.flat().every(Number.isFinite));
    close(Math.max(...shape.map(p=>p[0]))-Math.min(...shape.map(p=>p[0])),width);
    assert.deepEqual(shape.map(p=>p[1]),[20,30,40]);
  }
  assert.equal(points(M.arrowPath({x1:0,x2:1,y:0,strand:0})).length,4);
});
test('optional omitted scales equal explicit 100% without mutating the project', () => {
  const p=structuredClone(fixture), before=structuredClone(p), base=M.geometry(p), svg=SVG.generate(p);
  assert.deepEqual(p,before);Object.assign(p.visualization,{width_scale:1,height_scale:1,row_spacing_scale:1});
  assert.deepEqual(M.geometry(p),{...base,index:M.index(p)});assert.equal(SVG.generate(p),svg);
  close(base.rows.get('a').y,M.TOP);close(base.height,M.TOP+2*M.ROW_HEIGHT+95);
});
const cases=[[1,1,1],[.5,1,1],[2,1,1],[1,.5,1],[1,2,1],[.65,1.75,1],
  [1,1,.5],[1,1,1.5],[.65,1.75,1.5],[.25,4,.25],[4,.25,4]];
for (const [sx,sy,spacing] of cases) {
  test(`intrinsic arrow proportions at ${sx} × ${sy}, spacing ${spacing}`,()=>{
    // A long ordinary CDS must retain a 1:2 head-length:height ratio even when
    // its genomic span is compressed. This fails for a head based on widthScale.
    const s={x1:25,x2:25+400*sx,y:60*sy,sx,sy,strand:1};
    const forward=points(M.arrowPath(s)),reverse=points(M.arrowPath({...s,strand:-1}));
    assert.equal(forward.length,7);
    close(forward[3][0]-forward[2][0],10*sy);
    close((forward[3][0]-forward[2][0])/(forward[4][1]-forward[2][1]),.5);
    close(forward[5][1]-forward[1][1],14*sy);
    close(forward[0][0],s.x1);close(forward[3][0],s.x2);
    close(forward[0][1],forward[1][1]);close(forward[5][1],forward[6][1]);
    reverse.forEach(([x,y],i)=>{close(x,s.x1+s.x2-forward[i][0]);close(y,forward[i][1]);});
    for(const fraction of [.001,.5,1,1.001])for(const strand of [-1,1]) {
      const width=10*sy*fraction,shape=points(M.arrowPath({...s,x2:s.x1+width,strand}));
      assert.equal(shape.length,fraction<=1?3:7);
      close(Math.max(...shape.map(p=>p[0]))-Math.min(...shape.map(p=>p[0])),width);
      close(Math.max(...shape.map(p=>p[1]))-Math.min(...shape.map(p=>p[1])),20*sy);
      if(fraction>1)close(Math.abs(shape[3][0]-shape[2][0]),10*sy);
      assert.ok(shape.every(([x])=>x>=s.x1&&x<=s.x1+width));
    }
  });
  test(`layout and editable SVG at ${sx} × ${sy}, spacing ${spacing}`,()=>{
    const p=structuredClone(fixture), base=M.geometry(p), connections=M.connectors(p,base), bar=M.scaleBar(p,base);
    const biology=structuredClone(p);
    Object.assign(p.visualization,{width_scale:sx,height_scale:sy,row_spacing_scale:spacing,labels:true});
    const before=structuredClone(p), geo=M.geometry(p), svg=SVG.generate(p);
    assert.deepEqual(p,before);
    for(const key of ['clusters','layout','homolog_groups','pairwise_hits'])assert.deepEqual(p[key],biology[key]);
    close(geo.width,(base.width-M.LEFT-65)*sx+(M.LEFT+65)*sy);
    close(geo.height,(M.TOP+95)*sy+(p.clusters.length-1)*M.ROW_HEIGHT*sy*spacing);
    for(const [gid,shape] of geo.genes) shape.segments.forEach((segment,i)=>{
      const original=base.genes.get(gid).segments[i],vertices=points(M.arrowPath(segment));
      close(segment.x1,geo.left+(original.x1-base.left)*sx);
      close(segment.x2-segment.x1,(original.x2-original.x1)*sx);
      close(segment.y,M.TOP*sy+(original.y-M.TOP)*sy*spacing);
      const width=segment.x2-segment.x1;
      if(segment.strand) {
        assert.equal(vertices.length,width<=10*sy?3:7);
        if(vertices.length===7)close(Math.abs(vertices[3][0]-vertices[2][0]),10*sy);
      }
      assert.ok(svg.includes(`d="${M.arrowPath(segment)}"`));
    });
    M.connectors(p,geo).forEach((c,i)=>{
      assert.equal(c.query,connections[i].query);assert.equal(c.target,connections[i].target);
      c.points.forEach(([x],j)=>close(x,geo.left+(connections[i].points[j][0]-base.left)*sx));
      const top=geo.genes.get(c.query).row,bottom=geo.genes.get(c.target).row;
      close(c.points[0][1],top.y+11*sy);close(c.points[1][1],top.y+11*sy);
      close(c.points[2][1],bottom.y-11*sy);close(c.points[3][1],bottom.y-11*sy);
    });
    const scaledBar=M.scaleBar(p,geo);assert.equal(scaledBar.bp,bar.bp);
    close(scaledBar.end-scaledBar.x,(bar.end-bar.x)*sx);close(scaledBar.y,geo.height-28*sy);
    const root=svg.match(/<svg[^>]+>/)[0], bounds=root.match(/viewBox="([^"]+)"/)[1].split(' ').map(Number);
    close(+root.match(/width="([^"]+)"/)[1],bounds[2]);close(+root.match(/height="([^"]+)"/)[1],bounds[3]);
    assert.ok(bounds[0]+bounds[2]>=geo.width-1e-8);close(bounds[3],geo.height);
    const inside=([x,y])=>{assert.ok(x>=bounds[0]&&x<=bounds[0]+bounds[2],`x=${x}`);assert.ok(y>=0&&y<=bounds[3],`y=${y}`);};
    for(const m of svg.matchAll(/<path[^>]+d="([^"]+)"/g))points(m[1]).forEach(inside);
    for(const m of svg.matchAll(/points="([^"]+)"/g))m[1].split(' ').map(v=>v.split(',').map(Number)).forEach(inside);
    for(const m of svg.matchAll(/<line[^>]+>/g)) {
      const attr=k=>+m[0].match(new RegExp(`${k}="([^"]+)"`))[1];inside([attr('x1'),attr('y1')]);inside([attr('x2'),attr('y2')]);
    }
    assert.doesNotMatch(svg,/<(?:image|g|use|foreignObject|canvas)\b|style=|transform=|textLength=|lengthAdjust=/);
    for(const text of svg.matchAll(/<text[^>]+>/g)) {
      const fontSize=+text[0].match(/font-size="([^"]+)"/)[1];
      assert.ok([9,10,12].some(size=>Math.abs(fontSize-size*sy)<1e-8));
    }
  });
}
test('text attributes use natural font sizing, independent of Width and Row spacing',()=>{
  for(const sx of [.25,.5,1,2,4])for(const sy of [.25,.5,1,1.75,4])for(const spacing of [.25,1,4]) {
    assert.deepEqual(M.textAttributes(15,30,12,{sx,sy,spacing}),{x:15,y:30,'font-size':12*sy});
  }
});
test('Row spacing changes only row Y positions and the footer, not glyphs or X geometry',()=>{
  const p=structuredClone(fixture);Object.assign(p.visualization,{width_scale:.65,height_scale:1.75});
  const base=M.geometry(p);
  for(const spacing of [.25,.5,1.5,4]) {
    p.visualization.row_spacing_scale=spacing;const geo=M.geometry(p);
    close(geo.width,base.width);close(geo.rowHeight,base.rowHeight*spacing);
    for(const [i,cid]of p.layout.order.entries()) {
      const delta=i*base.rowHeight*(spacing-1);
      assert.deepEqual(geo.rows.get(cid),{...base.rows.get(cid),y:base.rows.get(cid).y+delta});
      for(const gene of geo.index.clusters.get(cid).cds)geo.genes.get(gene.id).segments.forEach((s,j)=>{
        const original=points(M.arrowPath(base.genes.get(gene.id).segments[j]));
        points(M.arrowPath(s)).forEach(([x,y],k)=>{close(x,original[k][0]);close(y,original[k][1]+delta);});
      });
    }
    close(geo.height-base.height,(p.clusters.length-1)*base.rowHeight*(spacing-1));
  }
});
test('scaled left-bound compensation preserves stationary rows during drag',()=>{
  const p=structuredClone(fixture);Object.assign(p.visualization,{width_scale:2,height_scale:.5});
  const before=M.geometry(p),cid=p.layout.order.find(id=>p.layout.rows[id].offset===before.minimum);
  p.layout.rows[cid].offset-=1234;const after=M.geometry(p),zoom=3.4;
  for(const c of p.clusters) {
    const dx=(after.rows.get(c.id).x-before.rows.get(c.id).x)*zoom+M.cameraShift(before,after,zoom);
    close(dx,c.id===cid?-1234*M.SCALE*2*zoom:0);
  }
});
test('export bounds include unusually long labels at anisotropic scales',()=>{
  const p=structuredClone(fixture);Object.assign(p.visualization,{labels:true,width_scale:.25,height_scale:4});
  p.clusters[0].cds[0].gene='Wide label '.repeat(40);
  const svg=SVG.generate(p),[left,,width,height]=svg.match(/viewBox="([^"]+)"/)[1].split(' ').map(Number);
  assert.ok(left<0);assert.ok(left+width>M.geometry(p).width);assert.equal(height,M.geometry(p).height);
});
test('one through fifty rows preserve proportional row and export bounds',()=>{
  for(const count of [1,2,20,50]) {
    const p=structuredClone(fixture),source=p.clusters[0];p.clusters=[];p.layout={order:[],rows:{}};
    p.homolog_groups=[];p.pairwise_hits=[];
    for(let i=0;i<count;i++){const id=`r${i}`;p.clusters.push({...source,id,cds:[]});p.layout.order.push(id);p.layout.rows[id]={offset:-i*137.3,orientation:1};}
    for(const spacing of [.25,.5,1,1.5,4]) {
      Object.assign(p.visualization,{width_scale:4,height_scale:.25,row_spacing_scale:spacing});
      const geo=M.geometry(p);close(geo.rows.get(`r${count-1}`).y,(M.TOP+(count-1)*M.ROW_HEIGHT*spacing)*.25);
      assert.ok(SVG.generate(p).includes(`height="${geo.height}"`));
    }
  }
});
