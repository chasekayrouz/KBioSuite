'use strict';
const $ = id => document.getElementById(id);
const token = document.querySelector('meta[name="desktop-token"]').content;
let inputs = [], desktop = false, initialized = false;
let outputOverridden = false;
let resources = null, refreshing = null, pollTimer, recheckAfterJob = false, setupPending = false;
async function api(path, data) {
  const response = await fetch('/api/' + path, {method:data ? 'POST':'GET',headers:{'X-Desktop-Token':token,'Content-Type':'application/json'},...(data?{body:JSON.stringify(data)}:{})});
  const result=await response.json(); if(!response.ok) throw Error(result.error || 'Operation failed');return result;
}
function setInputs(paths) {
  inputs=paths;
  $('inputs').replaceChildren(...paths.map(path=>{const li=document.createElement('li');li.textContent=path;return li;}));
  if(!outputOverridden){
    const first=paths[0]||'', split=Math.max(first.lastIndexOf('/'),first.lastIndexOf('\\'));
    $('output').value=!first?'':split<0?'.':split===0?first[0]:first.slice(0,split);
  }
}
window.setInputs=setInputs;
function refresh(checkResources=false) {
  if(refreshing)return checkResources===true?refreshing.then(()=>refresh(true)):refreshing;
  refreshing=refreshState(checkResources).finally(()=>{refreshing=null;});
  return refreshing;
}
async function refreshState(checkResources) {
  let busy=false;
  try {
    let state;
    if(!resources||checkResources===true){state=await api('status');resources=state;}
    else state={...resources,job:await api('job')};
    if(recheckAfterJob&&!state.job.busy){state=await api('status');resources=state;recheckAfterJob=false;}
    // Show a setup failure beside the setup button, not only further down the page.
    if(setupPending&&!state.job.busy){setupPending=false;$('setup-error').textContent=state.job.error||'';$('setup-error').hidden=!state.job.error;}
    const job=state.job;busy=job.busy;desktop=state.desktop;
    window.stopState?.(job);
    if(!initialized){$('threads').value=state.threads;initialized=true;}
    $('browser-input').hidden=desktop;$('choose-inputs').hidden=!desktop;$('choose-output').hidden=!desktop;
    const tools=Object.entries(state.status.executables).filter(([name,value])=>['hmmscan','hmmpress','rpsblast','blastdbcmd'].includes(name)&&!value).map(([name])=>name);
    const ready=['pfam','cdd'].every(name=>state.status.databases[name]?.ready)&&tools.length===0;
    await window.quickState?.(state,ready);
    $('databases').className=ready?'ready':'';
    $('databases').textContent=ready?`Pfam ${state.status.databases.pfam.version} · CDD ${state.status.databases.cdd.version} · ready`:`Setup needed${tools.length?' · missing tools: '+tools.join(', '):''}`;
    $('runtime').textContent=state.runtime;$('message').textContent=job.message;$('error').textContent=job.error||'';
    for(const id of ['run','setup','apply-inputs'])$(id).disabled=job.busy||(id==='run'&&!ready);
    window.setNativeState(desktop, job.busy);
    $('output').disabled=job.busy;$('threads').disabled=job.busy;$('genbank-output').disabled=job.busy;$('checkpoint').disabled=job.busy;
    $('checkpoint-status').textContent=job.checkpoint_dir?'Checkpoints: '+job.checkpoint_dir:'';
    $('results').replaceChildren(...job.results.map(r=>{
      const item=document.createElement('div');item.className='result';
      const title=document.createElement('strong');title.textContent=r.error?'Failed: '+r.input:`${r.proteins} proteins · ${r.input}`;item.append(title);
      if(r.error){const p=document.createElement('p');p.className='failure';p.textContent=r.error;item.append(p);}
      for(const warning of r.warnings){const p=document.createElement('p');p.className='warning';p.textContent='Warning: '+warning;item.append(p);}
      for(const output of r.outputs||(r.output?[r.output]:[])){
        const line=document.createElement('div');line.className='output-file';
        const path=document.createElement('span');path.textContent=output;line.append(path);
        if(desktop){const button=document.createElement('button');button.textContent='Show in Finder';button.disabled=!window.kannotateBridge.ready;button.onclick=()=>act(()=>window.pywebview.api.reveal(output));line.append(button);}
        item.append(line);
      }
      return item;
    }));
  } catch(error) {$('error').textContent=String(error.message||error);}
  finally{clearTimeout(pollTimer);pollTimer=setTimeout(refresh,busy?120:1500);}
}
window.refresh=refresh;
async function act(fn){try{$('error').textContent='';await fn();await refresh();}catch(error){$('error').textContent=String(error.message||error);}}
$('choose-inputs').onclick=()=>act(async()=>{const paths=await window.pywebview.api.choose_inputs();if(paths.length)setInputs(paths);});
$('choose-output').onclick=()=>act(async()=>{const path=await window.pywebview.api.choose_output();if(path){$('output').value=path;outputOverridden=true;}});
$('output').addEventListener('input',()=>{outputOverridden=true;});
$('apply-inputs').onclick=()=>setInputs($('paths').value.split('\n').map(p=>p.trim()).filter(Boolean));
$('run').onclick=()=>act(()=>api('annotate',{inputs,output:$('output').value,threads:Number($('threads').value),genbank_output:$('genbank-output').value,checkpoint:$('checkpoint').checked}));
$('setup').onclick=()=>{$('setup-error').hidden=true;act(async()=>{await api('setup',{});recheckAfterJob=true;setupPending=true;});};
$('refresh').onclick=()=>refresh(true);
refresh();
