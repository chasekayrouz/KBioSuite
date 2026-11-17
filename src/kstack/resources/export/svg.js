/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
/* Model-driven export. This module never reads a browser DOM or viewport. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('../renderer/model.js'));
  else root.CSSVG = factory(root.CSModel);
})(typeof globalThis !== 'undefined' ? globalThis : this, function (M) {
  'use strict';
  const escape = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;').replace(/'/g, '&apos;');
  // A disposable export projection, never a saved/validated project. Keep full
  // group membership/colors; only rows and their direct hits are filtered.
  function projectView(p, {order = p.layout.order, hidden = [], visualization = {}, offsets = {}} = {}) {
    const ids = new Set(p.layout.order), omitted = new Set(hidden);
    if (order.length !== ids.size || new Set(order).size !== ids.size || order.some(id => !ids.has(id)) ||
        [...omitted].some(id => !ids.has(id))) throw new Error('Invalid export cluster selection.');
    const visible = order.filter(id => !omitted.has(id));
    const clusters = p.clusters.filter(c => !omitted.has(c.id));
    const genes = new Set(clusters.flatMap(c => c.cds.map(g => g.id)));
    const settings = {...p.visualization};
    for (const key of ['width_scale', 'height_scale', 'row_spacing_scale', 'connector_identity', 'labels', 'connectors']) {
      if (!Object.hasOwn(visualization, key)) continue;
      const value = visualization[key];
      const valid = ['labels', 'connectors'].includes(key) ? typeof value === 'boolean' :
        Number.isFinite(value) && (key === 'connector_identity' ? value >= 0 && value <= 1 : value >= .25 && value <= 4);
      if (!valid) throw new Error('Invalid export presentation setting.');
      settings[key] = value;
    }
    const rows = Object.fromEntries(Object.entries(p.layout.rows).map(([id, row]) => [id, {...row}]));
    for (const [id, offset] of Object.entries(offsets)) {
      if (!ids.has(id) || !Number.isFinite(offset)) throw new Error('Invalid export alignment offset.');
      rows[id].offset = offset;
    }
    return {...p, clusters, layout: {...p.layout, rows, order: visible}, visualization: settings,
      pairwise_hits: p.pairwise_hits.filter(h => genes.has(h.query) && genes.has(h.target))};
  }
  function generate(p, options) {
    if (options) p = projectView(p, options);
    if (!p.layout.order.length) throw new Error('Show at least one cluster to export.');
    const geo = M.geometry(p), {sx, sy} = geo, elements = [];
    let left = 0, right = geo.width;
    const text = (x, y, value, size = 12, anchor = 'start') => {
      const attrs = {...M.textAttributes(x, y, size, geo), 'text-anchor': anchor};
      elements.push(`<text ${Object.entries(attrs).map(([k,v]) => `${k}="${v}"`).join(' ')} font-family="Arial,Helvetica,sans-serif" fill="#253248">${escape(value)}</text>`);
      // One em per code point is a conservative Arial/sans-serif label envelope.
      // Include long annotations without clipping or a hidden measurement tree.
      const extent = Array.from(String(value)).length * size * sy, start = x - (anchor === 'middle' ? extent/2 : 0);
      left = Math.min(left, start - 4 * sy); right = Math.max(right, start + extent + 4 * sy);
    };
    for (const connector of M.connectors(p, geo)) {
      elements.push(`<polygon data-kind="connector" data-query="${escape(connector.query)}" data-target="${escape(connector.target)}" points="${connector.points.map(x => x.join(',')).join(' ')}" fill="${connector.color}" stroke="none"/>`);
    }
    for (const cid of p.layout.order) {
      const row = geo.rows.get(cid), cluster = geo.index.clusters.get(cid);
      elements.push(`<line data-kind="baseline" x1="${row.x}" x2="${row.end}" y1="${row.y}" y2="${row.y}" stroke="#bac2cc" stroke-width="${sy}"/>`);
      text(12*sy, row.y + 4*sy, cluster.name);
      for (const gene of cluster.cds) {
        const shape = geo.genes.get(gene.id);
        shape.segments.forEach((s, i) => {
          // Every exon/origin segment is an independent editable gene path.
          elements.push(`<path data-kind="gene" data-gene="${escape(gene.id)}" data-segment="${i}" d="${M.arrowPath(s)}" fill="${M.color(p,gene,geo.index)}" stroke="#344054" stroke-width="${.7*sy}" stroke-linejoin="round"/>`);
          if (i && !gene.origin_spanning) {
            elements.push(`<line data-kind="compound-link" x1="${shape.segments[i-1].x2}" x2="${s.x1}" y1="${row.y}" y2="${row.y}" stroke="${M.color(p,gene,geo.index)}" stroke-width="${1.5*sy}" stroke-dasharray="${3*sx} ${2*sx}"/>`);
          }
        });
        if (gene.origin_spanning) {
          text(row.x - 9*sy, row.y - 15*sy, '↔', 10);
          text(row.end + 3*sy, row.y - 15*sy, '↔', 10);
        }
        if (p.visualization.labels && shape.segments.length) {
          const s = shape.segments.reduce((a,b) => b.x2-b.x1 > a.x2-a.x1 ? b : a);
          text((s.x1+s.x2)/2, row.y-18*sy, M.label(gene), 9, 'middle');
        }
      }
    }
    const {bp, x, y, end} = M.scaleBar(p, geo);
    elements.push(`<line data-kind="scale" x1="${x}" x2="${end}" y1="${y}" y2="${y}" stroke="#253248" stroke-width="${2*sy}"/>`);
    for (const tick of [x,end]) elements.push(`<line data-kind="scale" x1="${tick}" x2="${tick}" y1="${y-4*sy}" y2="${y+4*sy}" stroke="#253248" stroke-width="${sy}"/>`);
    text((x+end)/2,y+18*sy,`${bp/1000} kb`,12,'middle');
    return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${right-left}" height="${geo.height}" viewBox="${left} 0 ${right-left} ${geo.height}">\n<title>${escape(p.metadata.title || 'Kstack')}</title>\n${elements.join('\n')}\n</svg>\n`;
  }
  return {generate, projectView};
});
