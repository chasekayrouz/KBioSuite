'use strict';
(() => {
  const el=id=>document.getElementById(id), dialog=el('stop-dialog');
  let current=null, confirming=null;
  window.stopState=job=>{
    current=job;
    el('job-controls').hidden=!(job.can_stop||job.stopping||job.stopped);
    el('job-message').textContent=job.message;
    el('stop').hidden=!!job.stopped&&!job.busy;
    el('stop').disabled=!job.can_stop||job.stopping;
    el('stop').textContent=job.stopping?'Stopping…':'Stop…';
    if(dialog.open&&(!job.busy||job.job_id!==confirming)){
      confirming=null;dialog.close('continue');
    }
  };
  el('stop').onclick=()=>{
    if(!current?.can_stop||current.stopping||dialog.open)return;
    confirming=current.job_id;
    el('stop-file-warning').hidden=current.kind!=='files';
    el('stop-temporary-warning').textContent=current.checkpoint_dir?'Both Pfam and CDD searches will stop. Temporary data for the unfinished group will be removed. Completed checkpoints keep local queries and evidence.':'Both Pfam and CDD searches will stop. Temporary queries, search outputs, and new cache entries from this run will be removed.';
    el('stop-resume-warning').textContent=current.checkpoint_dir?'Completed checkpoints will be kept. The unfinished group will be discarded; run the same batch with checkpoints selected to resume.':'Unfinished work cannot be resumed.';
    dialog.returnValue='continue';dialog.showModal();
  };
  dialog.addEventListener('close',()=>{
    const jobId=confirming;confirming=null;
    if(dialog.returnValue!=='stop'||!jobId)return;
    el('stop').disabled=true;
    act(()=>api('stop',{job_id:jobId,confirmed:true}));
  });
})();
