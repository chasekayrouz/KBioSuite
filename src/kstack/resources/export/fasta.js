/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
(function(root,factory) {
  if (typeof module === 'object' && module.exports) module.exports=factory();
  else root.CSFASTA=factory();
})(typeof globalThis !== 'undefined' ? globalThis : this,function() {
  'use strict';
  const token = value => String(value).trim().replace(/[^A-Za-z0-9_.-]+/g,'_').replace(/^_+|_+$/g,'') || 'CDS';
  function generate(p, ids=p.fasta_selection) {
    const selected=new Set(ids), lines=[];
    for (const cid of p.layout.order) {
      const cluster=p.clusters.find(c=>c.id===cid);
      for (const gene of cluster.cds) {
        if (!selected.has(gene.id) || !gene.protein) continue;
        const product=gene.product.replace(/[\r\n>]/g,' ').replace(/\s+/g,' ').slice(0,100);
        lines.push(`>${token(cluster.name)}|${token(gene.locus_tag||gene.gene||'CDS')}|${gene.id}${product ? ' '+product : ''}`);
        lines.push(...gene.protein.match(/.{1,80}/g));
      }
    }
    return lines.length ? lines.join('\n')+'\n' : '';
  }
  return {generate,token};
});
