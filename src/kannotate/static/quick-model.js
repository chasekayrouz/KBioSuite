'use strict';
// Pure view model, shared by rendering and semantic interaction tests.
const QuickModel = (() => {
  function visibleHits(query, source='All', weaker=false) {
    return query.hits.filter(h => (source==='All'||h.source===source) && (weaker||h.significant));
  }
  function color(hit) {
    let hash=2166136261;
    for(const c of hit.source+':'+hit.accession) hash=Math.imul(hash^c.charCodeAt(0),16777619);
    return `hsl(${(hash>>>0)%360}, 48%, 78%)`;
  }
  function architecture(hits, length) {
    const blocks=[], groups=[];
    let offset=0;
    for(const source of ['Pfam','CDD']) {
      const intervals=hits.filter(h=>h.source===source).flatMap(hit=>hit.spans.map(span=>({hit,span})))
        .sort((a,b)=>(b.hit.score||0)-(a.hit.score||0)||(a.hit.evalue||0)-(b.hit.evalue||0)||a.span[0]-b.span[0]||a.span[1]-b.span[1]||a.hit.id.localeCompare(b.hit.id));
      const lanes=[];
      for(const {hit,span:[start,end]} of intervals) {
        let track=lanes.findIndex(intervals=>intervals.every(([a,b])=>b<start||end<a));
        if(track<0){track=lanes.length;lanes.push([]);}
        lanes[track].push([start,end]);
        blocks.push({id:hit.id,hit,start,end,track:offset+track,x:(start-1)/length,width:(end-start+1)/length,color:color(hit)});
      }
      if(lanes.length)groups.push({source,track:offset,count:lanes.length});
      offset+=lanes.length;
    }
    return {blocks,groups,tracks:offset};
  }
  class Selection {
    constructor(){this.query=null;this.hit=null;this.remembered=new Map();}
    selectQuery(query){this.query=query;this.hit=this.remembered.get(query.id)||null;}
    selectHit(id){
      if(!this.query?.hits.some(h=>h.id===id))return null;
      this.hit=id;this.remembered.set(this.query.id,id);return this.selected();
    }
    selected(){return this.query?.hits.find(h=>h.id===this.hit)||null;}
  }
  // Plain text, selected HSP only. Gaps never advance either residue axis.
  function alignmentText(alignment, columns=60) {
    const lines=[];
    let q=alignment.query_start, m=alignment.model_start;
    for(let i=0;i<alignment.query.length;i+=columns){
      const query=alignment.query.slice(i,i+columns), reference=alignment.model.slice(i,i+columns);
      const qn=query.replace(/-/g,'').length, mn=reference.replace(/-/g,'').length;
      const matches=[...query].map((c,j)=>c!=='-'&&c.toUpperCase()===reference[j].toUpperCase()?'|':' ').join('');
      const line=(name,start,sequence,end)=>`${name.padEnd(6)} ${String(start).padStart(5)}  ${sequence}  ${end}`;
      lines.push(line('Query',qn?q:'–',query,qn?q+qn-1:'–'),`              ${matches}`,line('Model',mn?m:'–',reference,mn?m+mn-1:'–'),'');
      q+=qn;m+=mn;
    }
    return lines.join('\n');
  }
  class ViewMemory {
    constructor(){this.views=new Map();}
    save(query,source,weaker,state){if(query)this.views.set(`${query.id}:${source}:${weaker}`,state);}
    restore(query,source,weaker){return this.views.get(`${query.id}:${source}:${weaker}`)||{graph:0,table:0,details:0};}
  }
  return {visibleHits,color,architecture,Selection,alignmentText,ViewMemory};
})();
if(typeof module!=='undefined')module.exports=QuickModel;
