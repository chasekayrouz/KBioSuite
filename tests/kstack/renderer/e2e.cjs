// Cross-language end-to-end smoke: actual GenBank → MMseqs2 → HTML → model → SVG/FASTA.
const {spawnSync} = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const SVG = require('../../../src/kstack/resources/export/svg.js'), FASTA = require('../../../src/kstack/resources/export/fasta.js'), M = require('../../../src/kstack/resources/renderer/model.js');
const root=path.resolve(__dirname,'../../..');
const python=process.env.KSTACK_PYTHON || path.join(root,'.venv/bin/python');
const inputs=['A. alliaceus CBS 536.65.gbk','A. mulundensis DSM 5745.gbk','A. versicolor CBS 583.65.gbk'].map(name=>({path:path.join(root,'tests/kstack/fixtures/realistic',name)}));
const run=request=>{
  const result=spawnSync(python,['-m','kstack.cli'],{cwd:root,env:{...process.env,PYTHONPATH:path.join(root,'src')},input:JSON.stringify(request),encoding:'utf8',maxBuffer:64*1024*1024});
  assert.equal(result.status,0,result.stderr||result.stdout);const response=JSON.parse(result.stdout);assert.ok(response.ok);return response.result;
};
const result=run({action:'generate',inputs});
const output=path.join(root,'.artifacts/e2e');fs.mkdirSync(output,{recursive:true});const file=path.join(output,'Stack.html');fs.writeFileSync(file,result.html);
const project=run({action:'open',path:file}).project;assert.deepEqual(project,result.project);
M.selectCluster(project,project.clusters[0].id);project.layout.rows[project.clusters[0].id].offset+=1250;
project.layout.rows[project.clusters[0].id].orientation*=-1;project.color_overrides.genes[project.clusters[0].cds[0].id]='#123456';project.visualization.labels=true;
Object.assign(project.visualization,{width_scale:.65,height_scale:1.75,row_spacing_scale:1.5});
const saved=run({action:'serialize',project});fs.writeFileSync(file,saved.html);assert.deepEqual(run({action:'open',path:file}).project,project);
const recomputed=run({action:'recompute',project,settings:{identity:.4,coverage:.8}}).project;
assert.equal(recomputed.analysis.identity,.4);assert.deepEqual(recomputed.layout,project.layout);
assert.deepEqual(recomputed.visualization,project.visualization);
const svg=SVG.generate(project),fasta=FASTA.generate(project);assert.ok(svg.includes('data-kind="gene"'));assert.ok(fasta.startsWith('>'));
fs.writeFileSync(path.join(output,'Figure.svg'),svg);fs.writeFileSync(path.join(output,'Selected.faa'),fasta);
const validation=spawnSync(python,['-c','import sys,xml.etree.ElementTree as E; r=E.parse(sys.argv[1]).getroot(); assert all(e.tag.split("}")[-1] not in ("image","mask","clipPath","g","use") for e in r.iter()); print("XML vector structure valid")',path.join(output,'Figure.svg')],{encoding:'utf8'});
assert.equal(validation.status,0,validation.stderr);
console.log(JSON.stringify({ok:true,clusters:project.clusters.length,genes:project.clusters.reduce((n,c)=>n+c.cds.length,0),pairs:project.pairwise_hits.length,selected:project.fasta_selection.length,elapsed:result.elapsed_seconds,svg:validation.stdout.trim()},null,2));
