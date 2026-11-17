const {test}=require('node:test');
const assert=require('node:assert/strict');
const vm=require('node:vm');
const fs=require('node:fs');
const path=require('node:path');
const QuickModel=require('../../../src/kannotate/static/quick-model.js');
const source=fs.readFileSync(path.join(__dirname,'../../../src/kannotate/static/quick.js'),'utf8');

async function page(importer) {
  const elements=new Map(),requests=[];
  const get=id=>{
    if(!elements.has(id))elements.set(id,{
      value:'',textContent:'',files:[],children:[],listeners:{},
      classList:{add(){},remove(){},toggle(){}},
      addEventListener(name,callback){this.listeners[name]=callback;},setAttribute(){},focus(){},
    });
    return elements.get(id);
  };
  const window={};
  const api=async(route,data)=>{
    requests.push({route,data});
    if(route==='quick/import')return importer(data);
    if(route==='quick/validate')return {count:1,lengths:[7]};
    return {revision:0};
  };
  vm.runInNewContext(source,{
    window,document:{getElementById:get,querySelectorAll:()=>[],addEventListener(){}},
    QuickModel,api,setTimeout(){},clearTimeout(){},ResizeObserver:class{observe(){}},
  });
  await window.quickState({job:{busy:false,quick_revision:0}},true);
  const choose=async file=>{get('quick-file').files=[file];await get('quick-file').listeners.change();};
  return {get,window,requests,choose};
}

test('file import fills editable FASTA, reports skipped CDSs and validates without searching',async()=>{
  const p=await page(async()=>({fasta:'>CDS\nMACDEFG\n',warnings:['Skipped 1 CDS.']}));
  await p.choose({name:'local.gbk',text:async()=>'LOCUS local'});
  assert.equal(p.get('quick-input').value,'>CDS\nMACDEFG\n');
  assert.match(p.get('quick-import-status').textContent,/Loaded local.gbk.*Skipped 1 CDS/);
  assert.equal(p.get('quick-search').disabled,false);
  assert.equal(p.get('quick-file').value,'');
  assert.deepEqual(p.requests.map(r=>r.route),['quick/import','quick/validate']);
  assert.equal(p.requests[0].data.text,'LOCUS local');
});

test('failed import preserves pasted queries and displays the limit explanation',async()=>{
  const p=await page(async()=>{throw Error('21 CDS proteins; up to 20. Use Annotate files.');});
  p.get('quick-input').value='>earlier\nMACDEFG';
  await p.choose({name:'large.gbk',text:async()=>''});
  assert.equal(p.get('quick-input').value,'>earlier\nMACDEFG');
  assert.match(p.get('quick-import-status').textContent,/21 CDS proteins.*Annotate files/);
  assert.equal(p.get('quick-file').disabled,false);
  assert.equal(p.get('quick-search').disabled,false);
});

test('import locks editing/search/reset, and a running job disables file selection',async()=>{
  let release;
  const p=await page(()=>new Promise(resolve=>{release=resolve;}));
  const importing=p.choose({name:'input.faa',text:async()=>'>p\nMACDEFG'});
  await new Promise(resolve=>setImmediate(resolve));
  for(const id of ['quick-search','quick-clear','quick-file','quick-input'])assert.equal(p.get(id).disabled,true);
  release({fasta:'>p\nMACDEFG\n',warnings:[]});await importing;
  assert.equal(p.get('quick-input').disabled,false);
  await p.window.quickState({job:{busy:true,quick_revision:0}},true);
  assert.equal(p.get('quick-file').disabled,true);
  const before=p.requests.length;
  await p.choose({name:'another.gbk',text:async()=>''});
  assert.equal(p.requests.length,before);
});
