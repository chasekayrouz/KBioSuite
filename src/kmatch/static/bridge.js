'use strict';
// pywebviewready, not DOM load, establishes the native API contract.
// The page learns desktop vs browser from /api/status and reports it with setNativeState.
(() => {
  const required = ['choose_msa', 'export_structure', 'reveal_run', 'export_alignment', 'open_pdb'];
  const state = window.kmatchBridge = {
    event_seen: false, ready: false, desktop: false, methods: {}, error: null
  };
  const failed = 'Native controls could not initialize. Close and reopen Kmatch; if this persists, rebuild/reinstall the suite.';
  const missing = 'Native controls did not connect. Close and reopen Kmatch; if this persists, rebuild/reinstall the suite.';
  let expired = false;
  function show() {
    if (expired && state.desktop && !state.ready) state.error ||= missing;
    const node = document.getElementById('bridge-error');
    node.textContent = state.desktop ? state.error || '' : '';
    node.hidden = !node.textContent;
  }
  window.setNativeState = (desktop) => {
    state.desktop = desktop === true; show();
  };
  // Native controls call through here; browser pages never need the bridge.
  window.kmatchNative = () => {
    if (!state.ready) throw new Error(state.error || 'Native controls are still connecting. Try again shortly.');
    return window.pywebview.api;
  };
  window.addEventListener('pywebviewready', () => {
    state.event_seen = true;
    state.methods = Object.fromEntries(required.map(name => [name,
      typeof window.pywebview?.api?.[name] === 'function']));
    state.ready = required.every(name => state.methods[name]);
    state.error = state.ready ? null : failed;
    show();
  });
  setTimeout(() => {expired = true; show();}, 15000);
  show();
})();
