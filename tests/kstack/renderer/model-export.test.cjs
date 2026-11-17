const test = require('node:test');
const assert = require('node:assert/strict');
const M = require('../../../src/kstack/resources/renderer/model.js');
const SVG = require('../../../src/kstack/resources/export/svg.js');
const FASTA = require('../../../src/kstack/resources/export/fasta.js');
const fixture = require('../fixtures/synthetic/project.json');
const p = () => structuredClone(fixture);

test('common bp scale preserves gene/intergenic distances', () => {
  const project=p(),g=M.geometry(project);
  const a=g.genes.get('a0').segments[0],b=g.genes.get('a1').segments[0];
  assert.equal(a.x2-a.x1,300*M.SCALE);assert.equal(b.x1-a.x1,4000*M.SCALE);
});
test('horizontal movement affects only figure offset', () => {
  const project=p(),before=structuredClone(project.clusters);project.layout.rows.a.offset+=500;
  const g=M.geometry(project);assert.deepEqual(project.clusters,before);
  assert.equal(g.rows.get('a').x-g.rows.get('b').x,1500*M.SCALE);
});
test('dragging the leftmost row keeps every other gene stationary on screen', () => {
  const project=p(),before=M.geometry(project),zoom=1.7;
  const cid=project.layout.order.find(id=>project.layout.rows[id].offset===before.minimum);
  project.layout.rows[cid].offset-=500;
  const after=M.geometry(project),camera=M.cameraShift(before,after,zoom);
  for(const c of project.clusters) {
    const old=before.genes.get(c.cds[0].id).segments[0].x1*zoom;
    const next=after.genes.get(c.cds[0].id).segments[0].x1*zoom+camera;
    assert.ok(Math.abs((next-old)-(c.id===cid ? -500*M.SCALE*zoom : 0))<1e-9);
  }
});
test('reversal preserves lengths and reverses arrows', () => {
  const project=p();project.layout.rows.a.orientation=-1;const g=M.geometry(project),s=g.genes.get('a0').segments[0];
  assert.equal(s.x1-g.rows.get('a').x,(10000-1300)*M.SCALE);assert.equal(s.strand,-1);assert.equal(s.x2-s.x1,7.5);
});
test('connector edges are direct rather than component combinations', () => {
  const project=p();project.layout.order=['a','c','b'];const connections=M.connectors(project);
  assert.equal(connections.length,1);assert.equal(connections[0].query,'c0');assert.equal(connections[0].target,'b0');
  assert.equal(project.clusters[0].cds[0].group_id,project.clusters[2].cds[0].group_id);
});
test('connectors require analysis identity AND shorter coverage', () => {
  const project=p();project.pairwise_hits[0].identity=.299;project.pairwise_hits[1].shorter_coverage=.69;
  assert.equal(M.connectors(project).length,1);
});
test('gapped alignment columns cannot create a failing-coverage connector', () => {
  const project=p(),hit=project.pairwise_hits[0];hit.alignment_length=150;hit.shorter_coverage=.699999;
  assert.ok(hit.alignment_length/Math.min(hit.query_length,hit.target_length)>1);
  assert.equal(M.passes(project,hit),false);assert.equal(M.connectors(project).length,2);
  assert.equal((SVG.generate(project).match(/data-kind="connector"/g)||[]).length,2);
});
test('connector coverage is inclusive at the configured boundary', () => {
  const project=p(),hit=project.pairwise_hits[0];hit.alignment_length=150;
  for(const [coverage,expected] of [[.699999,false],[.7-Number.EPSILON,false],[.7,true],[.700001,true],[1,true]]) {
    hit.shorter_coverage=coverage;assert.equal(M.passes(project,hit),expected);
  }
});
test('invalid stored coverage cannot pass a connector threshold', () => {
  const project=p(),hit=project.pairwise_hits[0];
  for(const coverage of [-.01,1.01,NaN,Infinity]) {
    hit.shorter_coverage=coverage;assert.equal(M.passes(project,hit),false);
  }
});
test('connector display filter leaves groups intact', () => {
  const project=p(),groups=structuredClone(project.homolog_groups);project.visualization.connector_identity=.65;
  assert.equal(M.connectors(project).length,0);assert.deepEqual(project.homolog_groups,groups);
});
test('connector toggle hides connectors', () => {const project=p();project.visualization.connectors=false;assert.equal(M.connectors(project).length,0);});
test('duplications retain separate direct connections', () => {
  const project=p(),copy=structuredClone(project.clusters[0].cds[0]);copy.id='a2';copy.segments[0].start+=500;copy.segments[0].end+=500;
  project.clusters[0].cds.push(copy);project.homolog_groups[0].members.push('a2');project.pairwise_hits.push({...project.pairwise_hits[0],query:'a2'});
  assert.equal(M.connectors(project).length,4);
});
test('generated group colors and singleton color', () => {
  const project=p();assert.equal(M.color(project,project.clusters[0].cds[0]),M.color(project,project.clusters[2].cds[0]));
});
test('gene override takes precedence and group recolor reaches all members', () => {
  const project=p(),gene=project.clusters[0].cds[0];M.setGeneColor(project,gene.id,'#123456');
  assert.equal(M.color(project,gene),'#123456');M.setGroupColor(project,gene.group_id,'#abcdef');
  assert.equal(M.color(project,gene),'#abcdef');assert.equal(M.color(project,project.clusters[2].cds[0]),'#abcdef');
});
test('reset group restores generated color', () => {
  const project=p(),gene=project.clusters[0].cds[0],base=M.color(project,gene);M.setGroupColor(project,gene.group_id,'#abcdef');M.setGeneColor(project,gene.id,'#123456');M.resetColor(project,gene,true);
  assert.equal(M.color(project,gene),base);
});
test('invalid color values rejected', () => {assert.throws(()=>M.setGeneColor(p(),'a0','javascript:bad'));});
test('vertical reordering preserves biological coordinates', () => {const project=p(),cs=structuredClone(project.clusters);assert.ok(M.moveRow(project,'a',1));assert.deepEqual(project.layout.order,['b','a','c']);assert.deepEqual(cs,project.clusters);});
test('out-of-bounds reordering is harmless', () => {assert.equal(M.moveRow(p(),'a',-1),false);});
test('selected-protein FASTA only includes selection', () => {
  const project=p();project.fasta_selection=['a0'];const text=FASTA.generate(project);assert.equal((text.match(/^>/gm)||[]).length,1);assert.match(text,/>a\|0\|a0 enzyme/);assert.ok(text.endsWith('\n'));
});
test('whole-cluster selection adds all proteins and can remove them', () => {
  const project=p();M.selectCluster(project,'a');assert.deepEqual(project.fasta_selection,['a0','a1']);
  assert.equal((FASTA.generate(project).match(/^>/gm)||[]).length,2);M.selectCluster(project,'a',false);assert.equal(FASTA.generate(project),'');
});
test('FASTA skips CDSs without interpretable protein', () => {const project=p();project.clusters[0].cds[0].protein='';M.selectCluster(project,'a');assert.deepEqual(project.fasta_selection,['a1']);});
test('FASTA headers unique with duplicate labels and safe multiline annotation', () => {
  const project=p();project.clusters[0].name='A weird / cluster';project.clusters[0].cds.forEach(g=>{g.gene='same';g.product='bad\n>header';});M.selectCluster(project,'a');
  const headers=FASTA.generate(project).split('\n').filter(x=>x.startsWith('>'));assert.equal(new Set(headers).size,2);assert.match(headers[0],/^>A_weird_cluster\|same\|a0 /);assert.equal(headers.length,2);
});
test('FASTA lines wrapped to at most 80 amino acids', () => {const project=p();project.fasta_selection=['a0'];const lines=FASTA.generate(project).split('\n').slice(1).filter(Boolean);assert.ok(lines.every(x=>x.length<=80));assert.equal(lines.join(''),project.clusters[0].cds[0].protein);});
test('SVG is deterministic and entirely independent top-level vectors', () => {
  const project=p(),svg=SVG.generate(project);assert.equal(svg,SVG.generate(project));
  for(const tag of ['image','mask','clipPath','g','use','foreignObject'])assert.doesNotMatch(svg,new RegExp(`<${tag}(?:[\\s>])`));
  assert.equal((svg.match(/data-kind="gene"/g)||[]).length,5);assert.equal((svg.match(/data-kind="connector"/g)||[]).length,3);
  assert.equal((svg.match(/data-kind="scale"/g)||[]).length,3);assert.doesNotMatch(svg,/transform=/);
});
test('SVG escapes all annotation/label XML', () => {const project=p();project.clusters[0].name='A<&"';project.visualization.labels=true;project.clusters[0].cds[0].gene='</text><script>bad';const svg=SVG.generate(project);assert.match(svg,/A&lt;&amp;&quot;/);assert.doesNotMatch(svg,/<script>/);});
test('SVG geometry independent of viewport state', () => {
  const project=p(),before=SVG.generate(project);project.viewport={zoom:20,x:1000,y:-300};assert.equal(SVG.generate(project),before);
});
test('SVG responds to current offsets/orientation/settings/colors', () => {
  const project=p(),before=SVG.generate(project);project.layout.rows.a.offset=500;project.layout.rows.a.orientation=-1;project.visualization.labels=true;M.setGeneColor(project,'a0','#112233');
  const after=SVG.generate(project);assert.notEqual(after,before);assert.match(after,/#112233/);assert.ok((after.match(/<text/g)||[]).length>(before.match(/<text/g)||[]).length);
});
test('SVG connector settings applied', () => {const project=p();project.visualization.connector_identity=.99;assert.doesNotMatch(SVG.generate(project),/data-kind="connector"/);});
test('compound/origin features remain one protein and independent segments', () => {
  const project=p(),gene=project.clusters[0].cds[0];gene.segments=[{start:0,end:150,strand:1},{start:9850,end:10000,strand:1}];gene.origin_spanning=true;
  const svg=SVG.generate(project);assert.equal((svg.match(/data-gene="a0"/g)||[]).length,2);project.fasta_selection=['a0'];assert.equal((FASTA.generate(project).match(/^>/gm)||[]).length,1);
});
test('exports do not mutate project data', () => {const project=p(),before=structuredClone(project);SVG.generate(project);FASTA.generate(project);assert.deepEqual(project,before);});

test('export preview removes gaps, includes only direct visible matches, and preserves colors', () => {
  const project=p(),before=structuredClone(project);
  project.color_overrides.genes.a0='#112233';
  const view=SVG.projectView(project,{hidden:['b']}),geo=M.geometry(view);
  assert.deepEqual(view.layout.order,['a','c']);
  assert.equal(geo.height,155+90);
  assert.equal(geo.rows.get('c').y-geo.rows.get('a').y,90);
  assert.equal(M.connectors(view).length,0,'Shared homolog group alone cannot create A–C connector');
  assert.deepEqual(view.homolog_groups,project.homolog_groups);
  const svg=SVG.generate(project,{hidden:['b']});
  assert.doesNotMatch(svg,/data-gene="b[01]"/);assert.match(svg,/#112233/);
  assert.deepEqual(project.clusters,before.clusters);assert.deepEqual(project.layout,before.layout);
  assert.deepEqual(project.pairwise_hits,before.pairwise_hits);
});
test('export preview reconnects newly adjacent rows only when a passing direct hit exists', () => {
  const project=p();project.pairwise_hits.push({...project.pairwise_hits[0],query:'a0',target:'c0'});
  assert.equal(M.connectors(SVG.projectView(project,{hidden:['b']})).length,1);
  assert.equal(M.connectors(SVG.projectView(project,{hidden:['b'],visualization:{connector_identity:1}})).length,0);
  project.pairwise_hits.at(-1).shorter_coverage=.1;
  assert.equal(M.connectors(SVG.projectView(project,{hidden:['b']})).length,0);
});
test('export preview reorders and scales without mutating saved figure, selection, or analysis', () => {
  const project=p(),before=structuredClone(project),options={order:['c','a','b'],hidden:['a'],
    visualization:{width_scale:2,height_scale:1.5,row_spacing_scale:2,labels:true,connectors:true,connector_identity:0}};
  const view=SVG.projectView(project,options),geo=M.geometry(view),svg=SVG.generate(project,options);
  assert.deepEqual(view.layout.order,['c','b']);assert.equal(geo.height,155*1.5+90*1.5*2);
  assert.equal(M.connectors(view)[0].query,'c0');assert.equal(M.connectors(view)[0].target,'b0');
  assert.deepEqual(view.layout.rows,project.layout.rows);assert.deepEqual(project,before);
  assert.equal(svg,SVG.generate(project,options));assert.doesNotMatch(svg,/data-gene="a[01]"/);
});
test('export preview handles one/no visible rows and rejects invalid selections/settings', () => {
  const project=p(),svg=SVG.generate(project,{hidden:['a','b']});
  assert.match(svg,/height="155"/);assert.doesNotMatch(svg,/data-kind="connector"/);
  assert.throws(()=>SVG.generate(project,{hidden:['a','b','c']}),/at least one/);
  for(const options of [{order:['a','a','c']},{order:['a','b']},{hidden:['unknown']},
    {visualization:{width_scale:Infinity}},{visualization:{labels:'false'}},{visualization:{connector_identity:-1}}]) {
    assert.throws(()=>SVG.generate(project,options));
  }
  assert.equal(SVG.generate(project,{}),SVG.generate(project));
});

test('auto-align fixes reference and uses oriented gene centers without changing biology', () => {
  const project=p();project.layout.rows.a.orientation=-1;project.layout.rows.a.offset=100;
  const before=structuredClone(project),offsets=M.alignmentOffsets(project,'a0');
  assert.deepEqual(offsets,{b:6800,c:5800});assert.deepEqual(project,before);
  for(const [id,offset]of Object.entries(offsets))project.layout.rows[id].offset=offset;
  assert.deepEqual(M.alignmentOffsets(project,'a0'),{});
  assert.deepEqual(project.clusters,before.clusters);assert.deepEqual(project.layout.order,before.layout.order);
});
test('auto-align duplicate choice uses direct bit score, independent of connector display', () => {
  const project=p(),copy=structuredClone(project.clusters[1].cds[0]);copy.id='b2';copy.segments=[{start:7000,end:7300,strand:1}];
  project.clusters[1].cds.push(copy);project.homolog_groups[0].members.push(copy.id);
  project.pairwise_hits.push({...project.pairwise_hits[0],target:'b2',bits:100,identity:.2});
  project.visualization.connectors=false;project.visualization.connector_identity=1;
  assert.equal(M.alignmentOffsets(project,'a0').b,-6000);
  // Reverse the stored edge: ranking must be symmetric.
  Object.assign(project.pairwise_hits.at(-1),{query:'b2',target:'a0'});
  assert.equal(M.alignmentOffsets(project,'a0').b,-6000);
});
test('auto-align leaves duplicate rows with no direct reference hit in place', () => {
  const project=p(),copy=structuredClone(project.clusters[2].cds[0]);copy.id='c1';
  project.clusters[2].cds.push(copy);project.homolog_groups[0].members.push(copy.id);
  project.layout.rows.c.offset=123;project.layout.rows.b.offset=456;
  assert.deepEqual(M.alignmentOffsets(project,'a0'),{b:-1000});
  assert.equal(project.layout.rows.c.offset,123);
});
test('auto-align aligns sole group members without direct reference hits and ignores missing groups', () => {
  const project=p();project.layout.rows.c.offset=123;
  assert.deepEqual(M.alignmentOffsets(project,'a0'),{c:-2000});
  assert.deepEqual(M.alignmentOffsets(project,'a1'),{});
  assert.deepEqual(M.alignmentOffsets(project,'missing'),{});
});
test('auto-align compound center weights segment lengths and preserves reverse orientation', () => {
  const project=p(),gene=project.clusters[0].cds[0];
  gene.segments=[{start:100,end:200,strand:1},{start:400,end:600,strand:1}];
  assert.equal(M.geneCenter(gene,project.clusters[0],{orientation:1}),1150/3);
  project.layout.rows.b.orientation=-1;
  assert.ok(Math.abs(M.alignmentOffsets(project,'a0').b-(1150/3-7850))<1e-6);
  gene.segments=[];assert.deepEqual(M.alignmentOffsets(project,'a0'),{});
});
test('auto-align duplicate score ties use evalue, identity, length, then stable gene id', () => {
  const project=p(),copy=structuredClone(project.clusters[1].cds[0]);copy.id='b2';copy.segments=[{start:7000,end:7300,strand:1}];
  project.clusters[1].cds.push(copy);const hit={...project.pairwise_hits[0],target:'b2'};project.pairwise_hits.push(hit);
  assert.deepEqual(M.alignmentOffsets(project,'a0'),{});
  hit.evalue=1e-20;assert.equal(M.alignmentOffsets(project,'a0').b,-6000);
  hit.evalue=1e-15;hit.identity=.7;assert.equal(M.alignmentOffsets(project,'a0').b,-6000);
  hit.identity=.6;hit.alignment_length=100;assert.equal(M.alignmentOffsets(project,'a0').b,-6000);
});
test('export offsets and auto-align operate on visible copies only', () => {
  const project=p(),before=structuredClone(project),view=SVG.projectView(project,{hidden:['b'],offsets:{c:123}});
  assert.deepEqual(M.alignmentOffsets(view,'a0'),{c:-2000});
  view.layout.rows.a.offset=999;assert.deepEqual(project,before);
  assert.throws(()=>SVG.projectView(project,{offsets:{a:Infinity}}));
  assert.throws(()=>SVG.projectView(project,{offsets:{unknown:123}}));
  assert.notEqual(SVG.generate(project,{offsets:{a:123}}),SVG.generate(project));
  assert.doesNotMatch(SVG.generate(project),/forward|reverse|CDS|10\.0 kb/);
});
