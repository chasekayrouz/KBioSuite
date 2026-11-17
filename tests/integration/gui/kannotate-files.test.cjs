const {test} = require('node:test');
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../../src/kannotate/static/app.js'), 'utf8');

async function page() {
  const elements = new Map(), requests = [], revealed = [];
  function element() {
    return {value:'', checked:false, textContent:'', children:[], listeners:{},
      append(...children){this.children.push(...children);},
      replaceChildren(...children){this.children=children;},
      addEventListener(name,fn){this.listeners[name]=fn;}};
  }
  const get = id => {if(!elements.has(id))elements.set(id,element());return elements.get(id);};
  get('genbank-output').value='genbank';
  let job={busy:false,results:[],message:'Ready',error:null}, chosen=null;
  const window={setNativeState(){},kannotateBridge:{ready:true},pywebview:{api:{
    async choose_output(){return chosen;}, async reveal(p){revealed.push(p);}
  }}};
  const state=()=>({desktop:true,threads:2,job,runtime:'local',status:{executables:{},databases:{pfam:{ready:true,version:'38.2'},cdd:{ready:true,version:'3.21'}}}});
  const fetch=async (url,options)=>{
    requests.push({url,data:options.body?JSON.parse(options.body):null});
    return {ok:true,json:async()=>url.endsWith('/status')?state():url.endsWith('/job')?job:{ok:true}};
  };
  vm.runInNewContext(source,{window,document:{getElementById:get,createElement:element,querySelector:()=>({content:'token'})},fetch,setTimeout(){},clearTimeout(){}});
  await window.refresh();
  return {get,window,requests,revealed,setJob:value=>{job=value;},choose:value=>{chosen=value;}};
}

test('destination follows first input through replacement, clearing, and relative/root paths',async()=>{
  const p=await page();
  for(const [paths,wanted] of [
    [['/first folder/a.gb','/other/b.gb'],'/first folder'],
    [['/other/b.gb'],'/other'],[[],''],[['/root.gb'],'/'],[['relative.faa'],'.'],
  ]){p.window.setInputs(paths);assert.equal(p.get('output').value,wanted);}
});
test('typed and native folder overrides survive input changes for the session',async()=>{
  const p=await page();p.window.setInputs(['/in/a.gb']);
  p.get('output').value='/typed folder';p.get('output').listeners.input();
  p.window.setInputs(['/next/b.gb']);assert.equal(p.get('output').value,'/typed folder');
  p.choose('/chosen folder');await p.get('choose-output').onclick();
  p.window.setInputs(['/last/c.gb']);assert.equal(p.get('output').value,'/chosen folder');
  const reopened=await page();reopened.window.setInputs(['/fresh/d.gb']);
  assert.equal(reopened.get('output').value,'/fresh');
});
test('cancelled output picker leaves automatic destination active',async()=>{
  const p=await page();p.window.setInputs(['/in/a.gb']);
  await p.get('choose-output').onclick();p.window.setInputs(['/next/b.gb']);
  assert.equal(p.get('output').value,'/next');
});
test('run submits GenBank choice and automatic folder; busy disables format control',async()=>{
  const p=await page();p.window.setInputs(['/in/a.gb']);
  for(const mode of ['genbank','csv','both']){
    p.get('genbank-output').value=mode;await p.get('run').onclick();
    const request=p.requests.filter(r=>r.url==='/api/annotate').at(-1);
    assert.equal(request.data.genbank_output,mode);assert.equal(request.data.output,'/in');
  }
  p.setJob({busy:true,results:[],message:'Searching'});await p.window.refresh();
  assert.equal(p.get('genbank-output').disabled,true);
});
test('both output paths are displayed and each Finder button reveals its own output',async()=>{
  const p=await page();
  const outputs=['/out/a.kannotate.gbk','/out/a.kannotate.csv'];
  p.setJob({busy:false,message:'Completed',results:[{input:'/in/a.gb',output:outputs[0],outputs,proteins:1,warnings:[]}]});
  await p.window.refresh();
  const lines=p.get('results').children[0].children.filter(e=>e.className==='output-file');
  assert.equal(lines.length,2);
  for(let i=0;i<2;i++){assert.equal(lines[i].children[0].textContent,outputs[i]);await lines[i].children[1].onclick();}
  assert.deepEqual(p.revealed,outputs);
});

test('batch checkpoint option is explicit, submitted, and frozen while running',async()=>{
  const p=await page();p.window.setInputs(['/in/a.gb']);
  p.get('checkpoint').checked=true;await p.get('run').onclick();
  assert.equal(p.requests.filter(r=>r.url==='/api/annotate').at(-1).data.checkpoint,true);
  p.setJob({busy:true,results:[],message:'Searching',checkpoint_dir:'/local/checkpoints'});await p.window.refresh();
  assert.equal(p.get('checkpoint').disabled,true);
  assert.equal(p.get('checkpoint-status').textContent,'Checkpoints: /local/checkpoints');
});
