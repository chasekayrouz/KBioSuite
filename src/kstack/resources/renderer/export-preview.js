/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
/* Disposable presentation state shared by the desktop and portable HTML. */
(function () {
  'use strict';
  const M = window.CSModel, SVG = window.CSSVG, $ = id => document.getElementById(id);
  function create({getProject, onExport}) {
    const dialog = $('export-preview'), list = $('export-clusters'), canvas = $('export-preview-canvas');
    const numeric = {width: ['width_scale', 25, 400, 100], height: ['height_scale', 25, 400, 100],
      spacing: ['row_spacing_scale', 25, 400, 100], identity: ['connector_identity', 0, 100, 100]};
    let base = null, options, nodes = new Map(), frame = 0, content = '', geometry, drag = null, opener;
    const announce = text => { $('export-preview-status').textContent = text; };
    const visibleOrder = () => options.order.filter(id => !options.hidden.has(id));
    function render() {
      if (frame) cancelAnimationFrame(frame);
      frame = 0;
      if (!base) return;
      const count = visibleOrder().length;
      $('export-preview-count').textContent = `${count} of ${base.clusters.length} clusters`;
      $('export-save').disabled = !count;
      canvas.classList.toggle('fit-width', $('export-fit').checked);
      if (!count) {
        content = ''; geometry = null;
        const empty = document.createElement('p'); empty.className = 'export-empty';
        empty.textContent = 'Show at least one cluster to preview and export.';
        canvas.replaceChildren(empty); return;
      }
      const view = SVG.projectView(base, options);
      geometry = M.geometry(view);
      content = SVG.generate(view);
      // Parse only our XML-escaped, model-generated SVG, never source HTML.
      const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
      svg.setAttribute('role', 'group'); svg.setAttribute('aria-label', 'SVG export figure');
      for (const gene of svg.querySelectorAll('path[data-gene]')) {
        gene.setAttribute('tabindex', '0'); gene.setAttribute('role', 'button');
        gene.setAttribute('aria-haspopup', 'menu');
        gene.setAttribute('aria-label', `${M.label(geometry.index.genes.get(gene.dataset.gene))}; gene actions`);
      }
      if (drag?.axis === 'x') {
        const box = [...drag.viewBox];
        box[0] += (drag.minimum - geometry.minimum) * geometry.scale;
        svg.setAttribute('viewBox', box.join(' '));
        svg.style.width = `${drag.width}px`; svg.style.height = `${drag.height}px`;
      }
      canvas.replaceChildren(document.importNode(svg, true));
    }
    function invalidate() { if (!frame) frame = requestAnimationFrame(render); }
    function syncList() {
      options.order.forEach((id, i) => {
        const node = nodes.get(id);
        node.querySelector('input').checked = !options.hidden.has(id);
        node.classList.toggle('omitted', options.hidden.has(id));
        node.classList.toggle('reordering', drag?.id === id);
        node.querySelector('[data-direction="-1"]').disabled = i === 0;
        node.querySelector('[data-direction="1"]').disabled = i === options.order.length - 1;
        if (list.children[i] !== node) list.insertBefore(node, list.children[i] || null);
      });
    }
    function move(id, next) {
      const current = options.order.indexOf(id);
      next = Math.max(0, Math.min(options.order.length - 1, next));
      if (current === next) return;
      options.order.splice(current, 1); options.order.splice(next, 0, id);
      syncList(); invalidate();
      announce(`${nodes.get(id).dataset.name}, position ${next + 1} of ${options.order.length}.`);
    }
    function button(label, text, handler) {
      const node = document.createElement('button'); node.textContent = text;
      node.setAttribute('aria-label', label); node.title = label; node.onclick = handler; return node;
    }
    function startDrag(event, id, kind) {
      if (event.button !== 0) return;
      event.preventDefault();
      const host = kind === 'list' ? list : canvas;
      host.setPointerCapture(event.pointerId);
      const svg = canvas.querySelector('svg'), rect = svg?.getBoundingClientRect();
      drag = {id, kind, pointer: event.pointerId, host, order: [...options.order], offsets: {...options.offsets},
        x: event.clientX, y: event.clientY, axis: kind === 'list' ? 'y' : null,
        offset: options.offsets[id] ?? base.layout.rows[id].offset, minimum: geometry?.minimum,
        units: svg?.getScreenCTM().a * geometry?.scale, width: rect?.width, height: rect?.height,
        viewBox: svg?.getAttribute('viewBox').split(' ').map(Number)};
      canvas.classList.add('reordering'); syncList();
    }
    function finishDrag(cancel = false) {
      if (!drag) return;
      const previous = drag; drag = null;
      if (cancel) { options.order = previous.order; options.offsets = previous.offsets; announce('Drag canceled.'); }
      invalidate();
      if (previous.host.hasPointerCapture(previous.pointer)) previous.host.releasePointerCapture(previous.pointer);
      canvas.classList.remove('reordering'); syncList();
    }
    function buildList() {
      nodes = new Map(); list.replaceChildren();
      for (const cluster of base.clusters) {
        const item = document.createElement('li'); item.dataset.cluster = cluster.id; item.dataset.name = cluster.name;
        const handle = button(`Reorder ${cluster.name}`, '⋮⋮', () => {}); handle.className = 'export-handle';
        handle.title = `Drag ${cluster.name} to reorder, or use arrow keys`;
        handle.onpointerdown = event => startDrag(event, cluster.id, 'list');
        handle.onkeydown = event => {
          if (!['ArrowUp', 'ArrowDown'].includes(event.key)) return;
          event.preventDefault(); move(cluster.id, options.order.indexOf(cluster.id) + (event.key === 'ArrowUp' ? -1 : 1));
          handle.focus();
        };
        const label = document.createElement('label'), checkbox = document.createElement('input'); checkbox.type = 'checkbox';
        checkbox.setAttribute('aria-label', `Include ${cluster.name}`);
        checkbox.onchange = () => {
          checkbox.checked ? options.hidden.delete(cluster.id) : options.hidden.add(cluster.id);
          syncList(); invalidate();
        };
        const name = document.createElement('span'); name.textContent = cluster.name; name.title = cluster.name;
        label.append(checkbox, name); item.append(handle, label);
        for (const direction of [-1, 1]) {
          const control = button(`Move ${cluster.name} ${direction < 0 ? 'up' : 'down'}`, direction < 0 ? '↑' : '↓', () => {
            move(cluster.id, options.order.indexOf(cluster.id) + direction);
            // The focused button may become disabled at a boundary.
            handle.focus();
          });
          control.dataset.direction = direction; item.append(control);
        }
        nodes.set(cluster.id, item);
      }
      syncList();
    }
    function reset() {
      finishDrag(); geneMenu.close();
      options = {order: [...base.layout.order], hidden: new Set(), offsets: {}, visualization: {...base.visualization}};
      for (const [id, [key, , , factor]] of Object.entries(numeric)) {
        $(`export-${id}`).value = (options.visualization[key] ?? (id === 'identity' ? 0 : 1)) * factor;
      }
      for (const key of ['labels', 'connectors']) $(`export-${key}`).checked = options.visualization[key];
      buildList(); render();
    }
    function open() {
      if (dialog.open) return;
      base = JSON.parse(JSON.stringify(getProject())); opener = document.activeElement;
      $('export-fit').checked = true; announce(''); reset();
      dialog.showModal(); $('export-preview-close').focus();
      canvas.scrollTop = 0; canvas.scrollLeft = 0;
    }
    function close() {
      geneMenu.close(); finishDrag();
      if (frame) cancelAnimationFrame(frame);
      frame = 0; base = null; content = ''; geometry = null;
      canvas.replaceChildren(); list.replaceChildren(); nodes.clear();
      if (dialog.open) { dialog.close(); if (opener?.isConnected) opener.focus(); }
    }
    for (const [id, [key, min, max, factor]] of Object.entries(numeric)) {
      const input = $(`export-${id}`);
      const update = () => { options.visualization[key] = Number(input.value) / factor; invalidate(); };
      input.oninput = () => { if (base && input.value !== '' && input.validity.valid) update(); };
      input.onchange = () => {
        if (!base) return;
        if (input.value !== '' && Number.isFinite(input.valueAsNumber)) {
          input.value = Math.max(min, Math.min(max, input.valueAsNumber)); update();
        } else input.value = (options.visualization[key] ?? 1) * factor;
      };
    }
    for (const key of ['labels', 'connectors']) $(`export-${key}`).onchange = () => {
      options.visualization[key] = $(`export-${key}`).checked; invalidate();
    };
    $('export-fit').onchange = () => { canvas.classList.toggle('fit-width', $('export-fit').checked); };
    $('export-show-all').onclick = () => { options.hidden.clear(); syncList(); invalidate(); };
    $('export-hide-all').onclick = () => { options.hidden = new Set(options.order); syncList(); invalidate(); };
    $('export-reset').onclick = () => { reset(); announce('Restored the stack’s current presentation.'); };
    for (const id of ['export-preview-close', 'export-cancel']) $(id).onclick = close;
    dialog.addEventListener('cancel', event => { event.preventDefault(); if (drag) finishDrag(true); else close(); });
    $('export-save').onclick = () => { render(); if (content) onExport(content); };
    canvas.onpointerdown = event => {
      const svg = canvas.querySelector('svg');
      if (!geometry || !svg?.contains(event.target)) return;
      const matrix = svg.getScreenCTM(); if (!matrix) return;
      const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse());
      const index = Math.round((point.y - M.TOP * geometry.sy) / geometry.rowHeight);
      const order = visibleOrder();
      if (index >= 0 && index < order.length) startDrag(event, order[index], 'figure');
    };
    window.addEventListener('pointermove', event => {
      if (!drag || drag.pointer !== event.pointerId) return;
      const dx = event.clientX - drag.x, dy = event.clientY - drag.y;
      if (!drag.axis) {
        if (Math.max(Math.abs(dx), Math.abs(dy)) < 4) return;
        drag.axis = Math.abs(dx) > Math.abs(dy) ? 'x' : 'y';
      }
      if (drag.axis === 'x') {
        options.offsets[drag.id] = Math.round((drag.offset + dx / drag.units) * 1e6) / 1e6;
        invalidate(); return;
      }
      if (Math.abs(dy) < 4) return;
      if (drag.kind === 'list') {
        const box = list.getBoundingClientRect();
        if (event.clientY < box.top + 24) list.scrollTop -= 18;
        if (event.clientY > box.bottom - 24) list.scrollTop += 18;
        const closest = [...list.children].reduce((best, node) => {
          const rect = node.getBoundingClientRect(), distance = Math.abs(event.clientY - (rect.top + rect.height / 2));
          return !best || distance < best.distance ? {node, distance} : best;
        }, null);
        if (closest) move(drag.id, options.order.indexOf(closest.node.dataset.cluster));
      } else {
        // Refresh geometry before hit testing when a preceding move awaits a frame.
        if (frame) render();
        const svg = canvas.querySelector('svg'), matrix = svg?.getScreenCTM(); if (!matrix) return;
        const box = canvas.getBoundingClientRect();
        if (event.clientY < box.top + 24) canvas.scrollTop -= 18;
        if (event.clientY > box.bottom - 24) canvas.scrollTop += 18;
        const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(matrix.inverse());
        const order = visibleOrder(), index = Math.max(0, Math.min(order.length - 1,
          Math.round((point.y - M.TOP * geometry.sy) / geometry.rowHeight)));
        move(drag.id, options.order.indexOf(order[index]));
      }
    });
    window.addEventListener('pointerup', event => { if (drag?.pointer === event.pointerId) finishDrag(); });
    window.addEventListener('pointercancel', event => { if (drag?.pointer === event.pointerId) finishDrag(true); });
    for (const host of [list, canvas]) host.addEventListener('lostpointercapture', () => finishDrag());
    const geneMenu = window.CSGeneMenu.attach({surface: canvas, host: dialog, onAlign: id => {
      const offsets = M.alignmentOffsets(SVG.projectView(base, options), id);
      if (Object.keys(offsets).length) { Object.assign(options.offsets, offsets); invalidate(); }
    }});
    return {open, close};
  }
  window.CSExportPreview = {create};
})();
