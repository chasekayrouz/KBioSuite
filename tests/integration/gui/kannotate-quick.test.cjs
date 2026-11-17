const {test}=require('node:test');
const assert=require('node:assert/strict');
const model=require('../../../src/kannotate/static/quick-model.js');
const hit=(id,source,spans,significant=true)=>({id,source,spans,significant,accession:id,name:id});
const query={id:'0',length:200,hits:[hit('a','Pfam',[[1,50],[101,150]]),hit('b','Pfam',[[40,90]]),hit('c','CDD',[[10,160]]),hit('weak','CDD',[[180,190]],false)]};
test('sources retain independent evidence and weaker matches are opt-in',()=>{
  assert.deepEqual(model.visibleHits(query).map(h=>h.id),['a','b','c']);
  assert.deepEqual(model.visibleHits(query,'CDD').map(h=>h.id),['c']);
  assert.equal(model.visibleHits(query,'All',true).length,4);
});
test('architecture uses proportional inclusive intervals, gap preservation and overlap tracks',()=>{
  const layout=model.architecture(model.visibleHits(query),200);
  const a=layout.blocks.filter(b=>b.id==='a');
  assert.deepEqual(a.map(b=>[b.x,b.width]),[[0,.25],[.5,.25]]);
  const b=layout.blocks.find(b=>b.id==='b');
  assert.notEqual(a[0].track,b.track);
  assert.equal(a[0].track,a[1].track);
  assert.notEqual(layout.blocks.find(b=>b.id==='c').track,a[0].track);
  assert.equal(model.color(a[0].hit),model.color(a[1].hit));
  assert.deepEqual(model.architecture([],100),{blocks:[],groups:[],tracks:0});
});
test('table and graphic share stable hit identity across sorting/filtering and query changes',()=>{
  const selection=new model.Selection();selection.selectQuery(query);
  const block=model.architecture(model.visibleHits(query),query.length).blocks.find(b=>b.id==='b');
  assert.equal(selection.selectHit(block.id),query.hits[1]); // graphic -> table entry
  const tableRow=model.visibleHits(query).slice().reverse()[0];
  assert.equal(selection.selectHit(tableRow.id).id,'c'); // sorted table -> graphic
  assert.equal(model.architecture(model.visibleHits(query,'CDD'),200).blocks[0].id,selection.hit);
  const second={id:'1',length:200,hits:[hit('other','Pfam',[[1,50]])]};
  selection.selectQuery(second);assert.equal(selection.hit,null);
  assert.equal(selection.selectHit('c'),null);
  selection.selectHit('other');selection.selectQuery(query);
  assert.equal(selection.hit,'c'); // each query keeps its own selected hit
});
test('significance ordering does not put nonoverlapping domains onto unnecessary tracks',()=>{
  const high={...hit('high','Pfam',[[1,30],[80,100]]),score:100};
  const middle={...hit('middle','Pfam',[[40,60]]),score:10};
  const layout=model.architecture([high,middle],100);
  assert.equal(layout.tracks,1);
  assert.equal(layout.blocks[2].width,.21);
});
test('gapped alignment chunks advance query and model coordinates independently',()=>{
  const alignment={query_start:11,query_end:17,model_start:3,model_end:9,query:'AC-DEFG-H',model:'ACD-E-GHH'};
  const text=model.alignmentText(alignment,4);
  assert.equal(text,[
    'Query     11  AC-D  13','              ||  ','Model      3  ACD-  5','',
    'Query     14  EFG-  16','              | | ','Model      6  E-GH  8','',
    'Query     17  H  17','              |','Model      9  H  9','',
  ].join('\n'));
  const gap=model.alignmentText({query_start:5,model_start:8,query:'----A',model:'ABCDE'},4);
  assert.match(gap,/Query\s+–  ----  –/);
  assert.match(gap,/Query\s+5  A  5/);
});
test('view memory keeps each query and source scroll position separate, including zero hits',()=>{
  const memory=new model.ViewMemory(),empty={id:'none',hits:[]};
  memory.save(query,'All',false,{graph:600,table:300,details:90,alignment:1});
  memory.save(query,'CDD',false,{graph:40,table:20,details:0,alignment:0});
  memory.save(empty,'All',false,{graph:0,table:0,details:0});
  assert.equal(memory.restore(query,'All',false).graph,600);
  assert.equal(memory.restore(query,'All',false).alignment,1);
  assert.equal(memory.restore(query,'CDD',false).graph,40);
  assert.equal(memory.restore(query,'All',true).graph,0);
  assert.equal(memory.restore(empty,'All',false).graph,0);
});
