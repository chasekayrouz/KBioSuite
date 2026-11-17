/* Copyright (c) 2026 Chase Kayrouz. MIT License. */
/* Gene actions shared by the main editor and the modal export preview. */
(function () {
  'use strict';
  function attach({surface, host, onAlign, enabled = () => true}) {
    const menu = document.createElement('div'), action = document.createElement('button');
    menu.className = 'gene-menu'; menu.hidden = true;
    menu.setAttribute('role', 'menu'); menu.setAttribute('aria-label', 'Gene actions');
    action.textContent = 'Auto-align'; action.setAttribute('role', 'menuitem');
    menu.append(action); host.append(menu);
    let reference = null, trigger = null;
    function close(restore = false) {
      if (menu.hidden) return;
      menu.hidden = true;
      if (restore && trigger?.isConnected) trigger.focus();
      reference = null;
    }
    function open(event, keyboard = false) {
      const gene = event.target.closest('path[data-gene]');
      if (!enabled() || !gene || !surface.contains(gene)) return;
      event.preventDefault(); event.stopPropagation();
      reference = gene.dataset.gene; trigger = gene;
      const rect = gene.getBoundingClientRect();
      const x = keyboard ? rect.left + rect.width / 2 : event.clientX;
      const y = keyboard ? rect.top + rect.height / 2 : event.clientY;
      menu.hidden = false;
      menu.style.left = `${Math.max(8, Math.min(x, window.innerWidth - menu.offsetWidth - 8))}px`;
      menu.style.top = `${Math.max(8, Math.min(y, window.innerHeight - menu.offsetHeight - 8))}px`;
      action.focus();
    }
    surface.addEventListener('contextmenu', event => open(event));
    surface.addEventListener('keydown', event => {
      if (event.key === 'ContextMenu' || (event.shiftKey && event.key === 'F10')) open(event, true);
    });
    action.onclick = () => {
      const id = reference; close();
      if (id !== null) {
        onAlign(id);
        requestAnimationFrame(() => {
          [...surface.querySelectorAll('path[data-gene]')].find(node => node.dataset.gene === id)?.focus({preventScroll:true});
        });
      }
    };
    menu.addEventListener('keydown', event => {
      if (event.key === 'Escape' || event.key === 'Tab') {
        event.preventDefault(); event.stopPropagation(); close(true);
      } else if (['ArrowUp', 'ArrowDown', 'Home', 'End'].includes(event.key)) {
        event.preventDefault(); action.focus();
      }
    });
    window.addEventListener('pointerdown', event => { if (!menu.contains(event.target)) close(); }, true);
    window.addEventListener('scroll', () => close(), true);
    window.addEventListener('resize', () => close());
    window.addEventListener('blur', () => close());
    surface.addEventListener('wheel', () => close(), {passive:true});
    return {close};
  }
  window.CSGeneMenu = {attach};
})();
