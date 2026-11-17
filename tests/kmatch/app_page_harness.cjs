'use strict';
// Loads src/kmatch/static/app.js against a minimal fake DOM in which only index.html's ids
// exist, answers /api requests from a JSON scenario (optionally overriding the job state that
// /api/status reports), presses the first dependency's setup button and prints what the page
// showed, asked and sent. Driven by tests/kmatch/test_smoke_page.py.
const fs = require('node:fs');
const vm = require('node:vm');
const [appPath, scenarioPath] = process.argv.slice(2);
const {ids, dependencies, job: reported} = JSON.parse(fs.readFileSync(scenarioPath, 'utf8'));
const known = new Set(ids);
function node(tag) {
  const listeners = {};
  return {
    tag, textContent: '', hidden: false, disabled: false, open: false, value: '', checked: false,
    className: '', title: '', dataset: {}, children: [],
    classList: {toggle() {}, add() {}, remove() {}, contains() {return false;}},
    addEventListener(type, callback) {(listeners[type] ||= []).push(callback);},
    async dispatch(type) {
      for (const callback of listeners[type] || []) await callback({stopPropagation() {}, preventDefault() {}});
    },
    append(...nodes) {this.children.push(...nodes);},
    replaceChildren(...nodes) {this.children = nodes;},
    querySelector() {return node('option');},
  };
}
const elements = new Map();
const document = {
  getElementById(id) {
    if (!known.has(id)) throw new Error('app.js referenced an element missing from index.html: #' + id);
    if (!elements.has(id)) elements.set(id, node(id));
    return elements.get(id);
  },
  querySelector(selector) {
    if (selector === 'meta[name=kmatch-token]') return {content: 'token'};
    throw new Error('Unexpected selector: ' + selector);
  },
  querySelectorAll() {return [];},
  createElement: node,
  createTextNode: text => ({textContent: text}),
};
const job = {busy: false, phase: 'ready', message: 'Ready', error: null, prediction: null, hits: [],
  hit_rankings: {}, directory: null, run_id: null, alignment_export: null, operation: null,
  pdb_details: {}, pdb_details_status: '', ...reported};
const status = {job, desktop: false, home: '/home', runs: '/home/runs',
  defaults: {predictor: 'boltz2', samples: 3},
  backends: [{id: 'boltz2', name: 'Boltz-2', label: 'Boltz-2 · Apple GPU', note: ''}],
  alignment_support: {pse: false, message: ''}};
const requests = [], confirms = [];
async function fetch(url, options) {
  requests.push({url, body: options.body});
  const body = url === '/api/status' ? status : url === '/api/dependencies' ? dependencies : {started: true};
  return {ok: true, json: async () => body};
}
const context = {
  document, fetch, URL, console,
  ResizeObserver: class {observe() {}},
  setInterval() {}, setTimeout() {}, clearTimeout() {},
  confirm(message) {confirms.push(message); return false;},
  setNativeState() {}, kmatchBridge: {ready: false},
};
context.window = context;
vm.createContext(context);
vm.runInContext(fs.readFileSync(appPath, 'utf8'), context, {filename: appPath});
const settle = async () => {for (let i = 0; i < 25; i++) await new Promise(resolve => setImmediate(resolve));};
(async () => {
  await settle();
  const shown = {summary: elements.get('dependency-summary').textContent,
    panel_open: elements.get('setup-panel').open};
  const card = elements.get('dependencies').children[0];
  const button = card.children.find(child => child.dataset.component);
  const before = requests.length;
  await button.dispatch('click');
  await settle();
  process.stdout.write(JSON.stringify({
    ...shown, component: button.dataset.component, confirms,
    setup: requests.slice(before).filter(item => item.url === '/api/setup').map(item => JSON.parse(item.body)),
    error: elements.get('error')?.textContent || '',
    setup_error: elements.get('setup-error')?.hidden ? '' : elements.get('setup-error')?.textContent || '',
  }));
})().catch(error => {console.error(error.stack); process.exit(1);});
