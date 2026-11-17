const {test} = require('node:test');
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../../src/kannotate/static/bridge.js'), 'utf8');
function page(methods = {}) {
  const elements = Object.fromEntries(['choose-inputs', 'choose-output', 'bridge-error'].map(id => [id, {disabled:false,textContent:''}]));
  const listeners = {};
  let timeout;
  const window = {pywebview:{api:methods},addEventListener:(name, callback)=>{listeners[name]=callback;}};
  vm.runInNewContext(source, {window, document:{getElementById:id=>elements[id]}, setTimeout:fn=>{timeout=fn;}});
  return {window,elements,ready:()=>listeners.pywebviewready(),timeout:()=>timeout()};
}
test('DOM load cannot enable native controls before pywebviewready', () => {
  const p = page({choose_inputs(){}, choose_output(){}, reveal(){}});
  p.window.setNativeState(true, false);
  assert.equal(p.elements['choose-inputs'].disabled, true);
  assert.equal(p.window.kannotateBridge.event_seen, false);
  p.ready();
  assert.equal(p.window.kannotateBridge.ready, true);
  assert.equal(p.elements['choose-output'].disabled, false);
  p.window.setNativeState(true, true);
  assert.equal(p.elements['choose-inputs'].disabled, true);
});
test('incomplete bridge remains disabled with persistent actionable error', () => {
  const p = page({choose_output(){}, reveal(){}});
  p.window.setNativeState(true,false);p.ready();
  assert.equal(p.window.kannotateBridge.ready,false);
  assert.equal(p.elements['choose-inputs'].disabled,true);
  assert.match(p.elements['bridge-error'].textContent,/Close and reopen/);
  p.window.setNativeState(true,false);
  assert.match(p.elements['bridge-error'].textContent,/reinstall/);
});
test('timeout is actionable and a later ready event recovers', () => {
  const p = page();p.window.setNativeState(true,false);p.timeout();
  assert.match(p.elements['bridge-error'].textContent,/did not connect/);
  p.window.pywebview.api = {choose_inputs(){}, choose_output(){}, reveal(){}};
  p.ready();assert.equal(p.elements['choose-inputs'].disabled,false);
  assert.equal(p.elements['bridge-error'].textContent,'');
});
test('browser/debug mode requires no native bridge and displays no timeout error', () => {
  const p = page();p.window.setNativeState(false,false);p.timeout();
  assert.equal(p.elements['bridge-error'].textContent,'');
  assert.equal(p.window.kannotateBridge.ready,false);
});
