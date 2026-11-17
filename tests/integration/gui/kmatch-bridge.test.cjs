const {test} = require('node:test');
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../../src/kmatch/static/bridge.js'), 'utf8');
const METHODS = ['choose_msa', 'export_structure', 'reveal_run', 'export_alignment', 'open_pdb'];
function api(names = METHODS) {return Object.fromEntries(names.map(name => [name, () => name]));}
function page(methods = {}) {
  const elements = {'bridge-error': {textContent: '', hidden: true}};
  const listeners = {};
  let timeout;
  const window = {pywebview: {api: methods}, addEventListener: (name, callback) => {listeners[name] = callback;}};
  vm.runInNewContext(source, {window, document: {getElementById: id => elements[id]}, setTimeout: fn => {timeout = fn;}});
  return {window, banner: elements['bridge-error'], ready: () => listeners.pywebviewready(), timeout: () => timeout()};
}
test('pywebviewready, not DOM load, establishes the native contract', () => {
  const p = page(api());
  p.window.setNativeState(true);
  assert.equal(p.window.kmatchBridge.event_seen, false);
  assert.equal(p.window.kmatchBridge.ready, false);
  assert.throws(() => p.window.kmatchNative(), /still connecting/);
  p.ready();
  const bridge = p.window.kmatchBridge;
  assert.equal(bridge.event_seen, true);
  assert.equal(bridge.ready, true);
  assert.equal(bridge.desktop, true);
  assert.equal(bridge.error, null);
  assert.deepEqual(Object.keys(bridge.methods).sort(), [...METHODS].sort());
  assert.ok(Object.values(bridge.methods).every(Boolean));
  assert.equal(p.window.kmatchNative().open_pdb(), 'open_pdb');
  assert.equal(p.banner.hidden, true);
  // The installed bridge check serializes this exact object.
  assert.deepEqual(Object.keys(JSON.parse(JSON.stringify(bridge))).sort(),
    ['desktop', 'error', 'event_seen', 'methods', 'ready']);
});
test('missing native methods leave the bridge unready with an actionable error', () => {
  const p = page(api(METHODS.filter(name => name !== 'export_alignment')));
  p.window.setNativeState(true); p.ready();
  const bridge = p.window.kmatchBridge;
  assert.equal(bridge.ready, false);
  assert.equal(bridge.methods.export_alignment, false);
  assert.equal(bridge.methods.choose_msa, true);
  assert.match(bridge.error, /Close and reopen Kmatch/);
  assert.match(p.banner.textContent, /reinstall/);
  assert.equal(p.banner.hidden, false);
  assert.throws(() => p.window.kmatchNative(), /could not initialize/);
  p.timeout();
  assert.match(p.banner.textContent, /could not initialize/);
});
test('timeout is actionable and a later ready event recovers', () => {
  const p = page(); p.window.setNativeState(true); p.timeout();
  assert.match(p.window.kmatchBridge.error, /did not connect/);
  assert.match(p.banner.textContent, /did not connect/);
  assert.throws(() => p.window.kmatchNative(), /did not connect/);
  p.window.pywebview.api = api(); p.ready();
  assert.equal(p.window.kmatchBridge.ready, true);
  assert.equal(p.window.kmatchBridge.error, null);
  assert.equal(p.banner.textContent, '');
  assert.equal(p.banner.hidden, true);
});
test('timeout before the page learns it is native still reports once desktop is known', () => {
  const p = page(); p.timeout();
  assert.equal(p.banner.textContent, '');
  p.window.setNativeState(true);
  assert.match(p.banner.textContent, /did not connect/);
});
test('browser mode requires no native bridge and shows no timeout error', () => {
  const p = page(); p.window.setNativeState(false); p.timeout();
  assert.equal(p.banner.textContent, '');
  assert.equal(p.banner.hidden, true);
  assert.equal(p.window.kmatchBridge.ready, false);
  assert.equal(p.window.kmatchBridge.desktop, false);
});
