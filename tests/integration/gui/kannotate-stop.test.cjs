const {test}=require('node:test');
const assert=require('node:assert/strict');
const vm=require('node:vm');
const fs=require('node:fs');
const path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../../src/kannotate/static/stop.js'),'utf8');
function page(){
  const elements=Object.fromEntries(['job-controls','job-message','stop','stop-dialog','stop-file-warning','stop-resume-warning','stop-temporary-warning'].map(id=>[id,{hidden:true,disabled:false,textContent:''}]));
  const listeners={},requests=[],dialog=elements['stop-dialog'];
  Object.assign(dialog,{open:false,returnValue:'',addEventListener:(name,fn)=>listeners[name]=fn,
    showModal(){this.open=true;},close(value){this.returnValue=value;this.open=false;listeners.close();}});
  const window={};
  vm.runInNewContext(source,{window,document:{getElementById:id=>elements[id]},
    act:fn=>fn(),api:(route,data)=>requests.push({route,...data})});
  return {window,elements,requests,dialog};
}
const running={job_id:'run-1',kind:'files',busy:true,can_stop:true,stopping:false,stopped:false,message:'Searching file 1 of 2'};
test('opening, declining or escaping confirmation never sends Stop',()=>{
  const p=page();p.window.stopState(running);p.elements.stop.onclick();
  assert.equal(p.dialog.open,true);assert.equal(p.requests.length,0);
  assert.equal(p.elements['stop-file-warning'].hidden,false);
  p.dialog.close('continue');assert.equal(p.requests.length,0);
  p.elements.stop.onclick();p.dialog.close('');assert.equal(p.requests.length,0);
});
test('confirmed Stop names the run and stays disabled during cleanup',()=>{
  const p=page();p.window.stopState(running);p.elements.stop.onclick();p.dialog.close('stop');
  assert.deepEqual(p.requests,[{route:'stop',job_id:'run-1',confirmed:true}]);
  assert.equal(p.elements.stop.disabled,true);
  p.window.stopState({...running,can_stop:false,stopping:true,message:'Stopping…'});
  p.elements.stop.onclick();assert.equal(p.dialog.open,false);
  assert.equal(p.elements.stop.textContent,'Stopping…');
  p.window.stopState({...running,busy:false,can_stop:false,stopped:true,message:'Stopped. Temporary data removed.'});
  assert.equal(p.elements.stop.hidden,true);
  assert.match(p.elements['job-message'].textContent,/Temporary data removed/);
});
test('completion or a newer run invalidates an open confirmation',()=>{
  const p=page();p.window.stopState(running);p.elements.stop.onclick();
  p.window.stopState({...running,job_id:'run-2'});
  assert.equal(p.dialog.open,false);p.dialog.close('stop');assert.equal(p.requests.length,0);
  p.elements.stop.onclick();p.window.stopState({...running,job_id:'run-2',busy:false,can_stop:false});
  assert.equal(p.dialog.open,false);assert.equal(p.requests.length,0);
});
test('Quick Search shares Stop and setup offers no cancellation',()=>{
  const p=page();p.window.stopState({...running,kind:'quick'});p.elements.stop.onclick();
  assert.equal(p.elements['stop-file-warning'].hidden,true);p.dialog.close('continue');
  p.window.stopState({...running,kind:null,can_stop:false});p.elements.stop.onclick();
  assert.equal(p.elements['job-controls'].hidden,true);assert.equal(p.dialog.open,false);
});

test('checkpointed file Stop explains what can resume; Quick Search stays disposable',()=>{
  const p=page();p.window.stopState({...running,checkpoint_dir:'/local/checkpoints'});p.elements.stop.onclick();
  assert.match(p.elements['stop-resume-warning'].textContent,/Completed checkpoints will be kept/);
  p.dialog.close('continue');p.window.stopState({...running,kind:'quick'});p.elements.stop.onclick();
  assert.match(p.elements['stop-resume-warning'].textContent,/cannot be resumed/);
});
