/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.CSModel = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const SCALE = 0.025, ROW_HEIGHT = 90, LEFT = 210, TOP = 60;
  const index = p => ({clusters: new Map(p.clusters.map(c => [c.id, c])),
    genes: new Map(p.clusters.flatMap(c => c.cds.map(g => [g.id, g]))),
    groups: new Map(p.homolog_groups.map(g => [g.id, g]))});
  function passes(p, h) {
    return h.identity + 1e-9 >= p.analysis.identity &&
      Number.isFinite(h.shorter_coverage) && h.shorter_coverage >= 0 && h.shorter_coverage <= 1 &&
      h.shorter_coverage >= p.analysis.coverage;
  }
  function color(p, gene, idx = index(p)) {
    return p.color_overrides.genes[gene.id] || p.color_overrides.groups[gene.group_id] ||
      idx.groups.get(gene.group_id)?.generated_color || '#a3aab4';
  }
  function geometry(p) {
    const idx = index(p);
    const sx = p.visualization.width_scale ?? 1, sy = p.visualization.height_scale ?? 1;
    const spacing = p.visualization.row_spacing_scale ?? 1;
    // Layout projects genomic X independently of intrinsic glyph size. Reserve
    // a natural-sized label column even when the genomic span is compressed.
    const left = LEFT * sy, rowHeight = ROW_HEIGHT * sy * spacing;
    const minimum = Math.min(...p.layout.order.map(id => p.layout.rows[id].offset));
    const maximum = Math.max(...p.layout.order.map(id => p.layout.rows[id].offset + idx.clusters.get(id).length));
    const rows = new Map(), genes = new Map();
    p.layout.order.forEach((id, i) => {
      const c = idx.clusters.get(id), state = p.layout.rows[id];
      const x = left + (state.offset - minimum) * SCALE * sx, y = TOP * sy + i * rowHeight;
      const row = {id, x, y, end: x + c.length * SCALE * sx, orientation: state.orientation};
      rows.set(id, row);
      c.cds.forEach(g => {
        const segments = g.segments.map(s => {
          const a = state.orientation === 1 ? s.start : c.length - s.end;
          const b = state.orientation === 1 ? s.end : c.length - s.start;
          return {x1: x + a * SCALE * sx, x2: x + b * SCALE * sx, y, sx, sy,
            strand: (s.strand || g.strand) * state.orientation};
        }).sort((a, b) => a.x1 - b.x1);
        genes.set(g.id, {id: g.id, segments, row});
      });
    });
    return {rows, genes, width: left + (maximum - minimum) * SCALE * sx + 65 * sy,
      height: (TOP + 95) * sy + (p.clusters.length - 1) * rowHeight,
      scale: SCALE * sx, rowHeight, left, sx, sy, spacing, minimum, index: idx};
  }
  function connectors(p, geo = geometry(p)) {
    if (!p.visualization.connectors) return [];
    const order = new Map(p.layout.order.map((id, i) => [id, i]));
    const result = [];
    const primary = gid => geo.genes.get(gid)?.segments.reduce((a, b) => !a || b.x2 - b.x1 > a.x2 - a.x1 ? b : a, null);
    for (const hit of p.pairwise_hits) {
      if (!passes(p, hit) || hit.identity + 1e-9 < p.visualization.connector_identity) continue;
      let a = geo.index.genes.get(hit.query), b = geo.index.genes.get(hit.target);
      if (Math.abs(order.get(a.cluster_id) - order.get(b.cluster_id)) !== 1) continue;
      if (order.get(a.cluster_id) > order.get(b.cluster_id)) [a, b] = [b, a];
      const top = primary(a.id), bottom = primary(b.id);
      if (!top || !bottom) continue;
      const intensity = Math.round(245 - 150 * hit.identity);
      result.push({query: a.id, target: b.id, identity: hit.identity,
        points: [[top.x1, top.y + 11 * geo.sy], [top.x2, top.y + 11 * geo.sy],
          [bottom.x2, bottom.y - 11 * geo.sy], [bottom.x1, bottom.y - 11 * geo.sy]],
        color: `rgb(${intensity},${intensity + 2},${Math.min(255, intensity + 6)})`});
    }
    return result;
  }
  function arrowPath(s) {
    const {x1, x2, y, strand, sy = 1} = s, h = 10 * sy, body = 7 * sy;
    if (!strand) return `M${x1},${y-h}L${x2},${y-h}L${x2},${y+h}L${x1},${y+h}Z`;
    const width = x2 - x1, head = Math.min(h, width);
    // Intrinsic head length is half the rendered gene height, independent of
    // genomic X projection. Short CDSs clamp to triangles within their span.
    // Reflect the same vertices so both strands have exactly mirrored geometry.
    const x = value => strand === 1 ? x1 + value : x2 - value;
    const points = width <= head ? [[0,-h],[width,0],[0,h]] :
      [[0,-body],[width-head,-body],[width-head,-h],[width,0],
        [width-head,h],[width-head,body],[0,body]];
    return points.map(([a,b],i) => `${i ? 'L' : 'M'}${x(a)},${y+b}`).join('') + 'Z';
  }
  function scaleBar(p, geo = geometry(p)) {
    const target = Math.max(...p.clusters.map(c => c.length)) / 8;
    const power = 10 ** Math.floor(Math.log10(Math.max(1, target)));
    const bp = [1, 2, 5, 10].map(v => v * power).reduce((a, b) => Math.abs(b-target)<Math.abs(a-target) ? b : a);
    return {bp, x: geo.left, y: geo.height - 28 * geo.sy, end: geo.left + bp * geo.scale};
  }
  function cameraShift(before, after, zoom) {
    return (after.minimum - before.minimum) * after.scale * zoom;
  }
  // Native font sizing is uniform within each glyph. Width scale affects the
  // supplied layout anchor, never the font's width:height proportions.
  function textAttributes(x, y, size, geo) {
    return {x, y, 'font-size': size * geo.sy};
  }
  const label = g => g.gene || g.locus_tag || g.product || g.id;
  const coordinates = g => g.segments.length ? g.segments.map(s => `${s.start + 1}–${s.end}`).join(', ') : g.location;
  function selectCluster(p, cid, selected = true) {
    const ids = p.clusters.find(c => c.id === cid).cds.filter(g => g.protein).map(g => g.id);
    const current = new Set(p.fasta_selection);
    ids.forEach(id => selected ? current.add(id) : current.delete(id));
    p.fasta_selection = [...current].sort();
  }
  function setGeneColor(p, gid, value) {
    if (!/^#[0-9a-f]{6}$/i.test(value)) throw new Error('Choose a valid RGB color.');
    p.color_overrides.genes[gid] = value;
  }
  function setGroupColor(p, gid, value) {
    if (!/^#[0-9a-f]{6}$/i.test(value)) throw new Error('Choose a valid RGB color.');
    const group = p.homolog_groups.find(g => g.id === gid);
    if (!group) throw new Error('Unknown homolog group.');
    p.color_overrides.groups[gid] = value;
    group.members.forEach(id => delete p.color_overrides.genes[id]);
  }
  function resetColor(p, gene, whole = false) {
    delete p.color_overrides.genes[gene.id];
    if (whole) {
      delete p.color_overrides.groups[gene.group_id];
      p.homolog_groups.find(g => g.id === gene.group_id).members.forEach(id => delete p.color_overrides.genes[id]);
    }
  }
  function moveRow(p, cid, delta) {
    const i = p.layout.order.indexOf(cid), j = i + delta;
    if (i < 0 || j < 0 || j >= p.layout.order.length) return false;
    [p.layout.order[i], p.layout.order[j]] = [p.layout.order[j], p.layout.order[i]];
    return true;
  }
  // Segment-length-weighted genomic center, matching the backend layout rule.
  // Presentation scales never participate in choosing or aligning genes.
  function geneCenter(gene, cluster, row) {
    const length = gene.segments.reduce((n, s) => n + s.end - s.start, 0);
    if (!length) return null;
    const center = gene.segments.reduce((n, s) => n + (s.start + s.end) / 2 * (s.end - s.start), 0) / length;
    return row.orientation === -1 ? cluster.length - center : center;
  }
  function alignmentOffsets(p, referenceId) {
    const idx = index(p), reference = idx.genes.get(referenceId), offsets = {};
    if (!reference || !p.layout.order.includes(reference.cluster_id)) return offsets;
    const referenceRow = p.layout.rows[reference.cluster_id];
    const center = geneCenter(reference, idx.clusters.get(reference.cluster_id), referenceRow);
    if (center === null) return offsets;
    const anchor = referenceRow.offset + center, hits = new Map();
    for (const hit of p.pairwise_hits) {
      if (hit.query === referenceId) hits.set(hit.target, hit);
      else if (hit.target === referenceId) hits.set(hit.query, hit);
    }
    for (const cid of p.layout.order) {
      if (cid === reference.cluster_id) continue;
      const cluster = idx.clusters.get(cid), row = p.layout.rows[cid];
      const candidates = cluster.cds.filter(g => g.group_id === reference.group_id);
      let chosen = candidates[0];
      if (candidates.length > 1) {
        chosen = candidates.filter(g => hits.has(g.id)).sort((a, b) => {
          const x = hits.get(a.id), y = hits.get(b.id);
          return y.bits - x.bits || x.evalue - y.evalue || y.identity - x.identity ||
            y.alignment_length - x.alignment_length || (a.id < b.id ? -1 : a.id > b.id ? 1 : 0);
        })[0];
      }
      if (!chosen) continue;
      const target = geneCenter(chosen, cluster, row);
      if (target === null) continue;
      const offset = Math.round((anchor - target) * 1e6) / 1e6;
      if (Math.abs(offset - row.offset) > 1e-6) offsets[cid] = offset;
    }
    return offsets;
  }
  return {SCALE, ROW_HEIGHT, LEFT, TOP, index, passes, color, geometry, connectors, arrowPath, scaleBar, cameraShift, textAttributes,
    label, coordinates, selectCluster, setGeneColor, setGroupColor, resetColor, moveRow, geneCenter, alignmentOffsets};
});
