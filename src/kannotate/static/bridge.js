'use strict';
// pywebviewready, not DOM load, establishes the native API contract.
(() => {
  const required = ['choose_inputs', 'choose_output', 'reveal'];
  const state = window.kannotateBridge = {
    event_seen: false, ready: false, methods: {}, error: null, desktop: false, busy: false
  };
  function controls() {
    for (const id of ['choose-inputs', 'choose-output']) {
      document.getElementById(id).disabled = !state.ready || state.busy;
    }
    document.getElementById('bridge-error').textContent = state.desktop ? state.error || '' : '';
  }
  window.setNativeState = (desktop, busy) => {
    state.desktop = desktop; state.busy = busy; controls();
  };
  window.addEventListener('pywebviewready', () => {
    state.event_seen = true;
    state.methods = Object.fromEntries(required.map(name => [name,
      typeof window.pywebview?.api?.[name] === 'function']));
    state.ready = required.every(name => state.methods[name]);
    state.error = state.ready ? null :
      'Native file controls could not initialize. Close and reopen Kannotate; if this persists, rebuild/reinstall the suite.';
    controls();
  });
  setTimeout(() => {
    if (state.desktop && !state.ready) {
      state.error ||= 'Native file controls did not connect. Close and reopen Kannotate; if this persists, rebuild/reinstall the suite.';
      controls();
    }
  }, 15000);
  controls();
})();
