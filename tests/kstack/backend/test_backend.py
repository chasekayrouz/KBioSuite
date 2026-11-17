import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from kstack import __version__
from kstack.parser import parse_genbank, inspect_inputs, ImportError
from kstack.homology import passing, connected_components, search, parse_alignment, shorter_protein_coverage
from kstack.layout import automatic_layout
from kstack.colors import generate_colors
from kstack.model import create_project, recompute, validate
from kstack.persistence import to_html, from_html, load
from kstack.cli import handle

ROOT=Path(__file__).resolve().parents[1]
SETTINGS={'identity':.3,'coverage':.7}

def cluster(cid, values):
    return {'id':cid,'name':cid,'source_name':cid+'.gbk','accession':cid,'description':'Synthetic test cluster',
            'topology':'linear','annotations':{},'warnings':[],'length':10000,'cds':[
        {'id':f'{cid}{i}','cluster_id':cid,'group_id':group,'segments':[{'start':start,'end':start+300,'strand':strand}],
         'cds_length':300,'strand':strand,'protein':'M'+'A'*99,'protein_length':100,'gene':str(i),
         'product':'enzyme','locus_tag':'','warnings':[],'partial':False,'qualifiers':{},'location':str(start),'translation_source':'translation','origin_spanning':False}
        for i,(group,start,strand) in enumerate(values)]}

def hit(q,t,identity=.6,length=90,qlen=100,tlen=100,qcov=.9,tcov=.9):
    if q>t:q,t,qlen,tlen,qcov,tcov=t,q,tlen,qlen,tcov,qcov
    return {'query':q,'target':t,'identity':identity,'alignment_length':length,'query_length':qlen,'target_length':tlen,
            'query_coverage':qcov,'target_coverage':tcov,
            'shorter_coverage':qcov if qlen<tlen else tcov if tlen<qlen else min(qcov,tcov),'bits':80.,'evalue':1e-15}

def simple_project():
    cs=[cluster('a',[('h1',1000,1),('h2',5000,-1)]),cluster('b',[('h1',2000,-1),('h2',6000,1)]),cluster('c',[('h1',3000,1)])]
    hits=[hit('a0','b0'),hit('b0','c0'),hit('a1','b1')]
    groups=connected_components(cs,hits,SETTINGS);generate_colors(groups)
    return {'format':'Kstack','schema_version':1,'app_version':__version__,'metadata':{'title':'Test','created':'2026-01-01','modified':'2026-01-01'},
            'clusters':cs,'analysis':{**SETTINGS,'preserve_input':False},'visualization':{'connectors':True,'connector_identity':.3,'labels':False},
            'homolog_groups':groups,'pairwise_hits':hits,'layout':automatic_layout(cs),
            'color_overrides':{'genes':{},'groups':{}},'fasta_selection':[]}

class ParsingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.c=parse_genbank(ROOT/'fixtures/synthetic/edge-cases.gbk')
    def test_normal_translation_and_coordinates(self):
        g=self.c['cds'][0];self.assertEqual(g['protein'],'M'+'A'*99);self.assertEqual(g['segments'],[{'start':90,'end':393,'strand':1}])
    def test_minus_strand_missing_translation(self):
        g=self.c['cds'][1];self.assertEqual(g['strand'],-1);self.assertEqual(g['protein'],'M'+'C'*99);self.assertEqual(g['translation_source'],'calculated')
    def test_hypothetical_retained(self): self.assertEqual(self.c['cds'][0]['product'],'hypothetical protein')
    def test_partial_retained(self):
        g=self.c['cds'][3];self.assertTrue(g['partial']);self.assertEqual(g['protein'],'M'+'P'*99);self.assertTrue(g['warnings'])
    def test_compound_protein(self):
        g=self.c['cds'][2];self.assertEqual(len(g['segments']),2);self.assertEqual(g['cds_length'],303);self.assertEqual(g['protein'],'M'+'G'*99)
    def test_malformed_retained_with_warning(self):
        g=self.c['cds'][4];self.assertTrue(g['warnings']);self.assertEqual(g['segments'],[]);self.assertEqual(len(self.c['cds']),5);self.assertEqual(g['translation_source'],'unavailable')
    def test_non_cds_not_visualized(self): self.assertFalse(any(g['product']=='not shown' for g in self.c['cds']))
    def test_multirecord_rejection(self):
        with self.assertRaisesRegex(ImportError,'found 3'): parse_genbank(ROOT/'fixtures/upstream/A. burnettii MST-FP2249.gbk')
    def test_over_100kb_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            r=SeqRecord(Seq('A'*100001),id='too_long');r.annotations['molecule_type']='DNA';p=Path(d)/'long.gbk';SeqIO.write(r,p,'genbank')
            with self.assertRaisesRegex(ImportError,'100,001'): parse_genbank(p)
    def test_empty_file_rejected(self):
        with tempfile.NamedTemporaryFile() as f:
            with self.assertRaisesRegex(ImportError,'found 0'): parse_genbank(f.name)
    def test_circular_origin_plus(self):
        c=parse_genbank(ROOT/'fixtures/synthetic/circular-origin.gbk');g=c['cds'][0]
        self.assertEqual(c['topology'],'circular');self.assertTrue(g['origin_spanning']);self.assertEqual(g['protein'],'M'+'V'*199)
        self.assertEqual(g['cds_length'],603);self.assertEqual(g['segments'][0]['start'],2700)
    def test_circular_origin_minus(self):
        c=parse_genbank(ROOT/'fixtures/synthetic/circular-origin-minus.gbk');g=c['cds'][0]
        self.assertEqual(g['protein'],'M'+'V'*199);self.assertEqual(g['strand'],-1)
    def test_nucleotide_absent(self):
        self.assertNotIn('sequence',self.c);self.assertNotIn('nucleotide',json.dumps(self.c))
    def test_ids_independent_of_display_name(self):
        p=ROOT/'fixtures/synthetic/edge-cases.gbk';self.assertEqual(parse_genbank(p,'x')['id'],parse_genbank(p,'y')['id'])
        self.assertNotEqual(parse_genbank(p,slot=0)['id'],parse_genbank(p,slot=1)['id'])
    def test_inspect_reports_error_without_discarding_other_files(self):
        r=inspect_inputs([str(ROOT/'fixtures/synthetic/edge-cases.gbk'),'/tmp/no-kstack-file']);self.assertIsNone(r[0]['error']);self.assertTrue(r[1]['error'])
    def test_input_count_limit(self):
        with self.assertRaises(ImportError):inspect_inputs([])
        with self.assertRaises(ImportError):inspect_inputs(['x']*51)
    def test_realistic_reference_biological_results(self):
        for p in sorted((ROOT/'fixtures/realistic').glob('*.gbk')):
            record=SeqIO.read(p,'genbank');parsed=parse_genbank(p);features=[f for f in record.features if f.type=='CDS']
            self.assertEqual(len(features),len(parsed['cds']))
            for feature,gene in zip(features,parsed['cds']):
                self.assertEqual(gene['strand'],feature.location.strand)
                self.assertEqual(gene['cds_length'],len(feature.location))
                self.assertEqual(gene['segments'][0]['start'],int(feature.location.parts[0].start))
                protein=feature.qualifiers.get('translation',[str(feature.extract(record.seq).translate())])[0].rstrip('*')
                self.assertEqual(gene['protein'],protein)

class HomologyTests(unittest.TestCase):
    def test_passing(self):self.assertTrue(passing(hit('a','b'),SETTINGS))
    def test_identity_failure(self):self.assertFalse(passing(hit('a','b',identity=.2999),SETTINGS))
    def test_coverage_failure(self):self.assertFalse(passing(hit('a','b',length=69,qcov=.69,tcov=.69),SETTINGS))
    def test_inclusive_threshold_and_shorter_coverage(self):self.assertTrue(passing(hit('a','b',identity=.3,length=70,tlen=500,qcov=.7,tcov=.14),SETTINGS))
    def test_gapped_alignment_columns_do_not_inflate_coverage(self):
        h=parse_alignment('a\tb\t0.6\t150\t100\t200\t0.69\t0.46\t80\t1e-15')
        self.assertGreater(h['alignment_length']/min(h['query_length'],h['target_length']),1)
        self.assertEqual(h['shorter_coverage'],.69);self.assertFalse(passing(h,SETTINGS))
    def test_query_shorter_uses_query_sequence_coverage(self):
        h=parse_alignment('a\tb\t0.6\t90\t100\t200\t0.8\t0.4\t80\t1e-15')
        self.assertEqual(h['shorter_coverage'],.8)
    def test_target_shorter_uses_target_sequence_coverage(self):
        h=parse_alignment('a\tb\t0.6\t90\t200\t100\t0.4\t0.8\t80\t1e-15')
        self.assertEqual(h['shorter_coverage'],.8)
    def test_equal_lengths_use_conservative_symmetric_coverage(self):
        h=parse_alignment('a\tb\t0.6\t120\t100\t100\t0.9\t0.69\t80\t1e-15')
        self.assertEqual(h['shorter_coverage'],.69);self.assertFalse(passing(h,SETTINGS))
    def test_reversing_query_target_preserves_coverage_and_decision(self):
        for qlen,tlen,qcov,tcov in ((100,200,.7,.35),(200,100,.35,.7),(100,100,.9,.69),(100,100,.9,.7)):
            with self.subTest(qlen=qlen,tlen=tlen,qcov=qcov,tcov=tcov):
                a=parse_alignment(f'a\tb\t0.6\t120\t{qlen}\t{tlen}\t{qcov}\t{tcov}\t80\t1e-15')
                b=parse_alignment(f'b\ta\t0.6\t120\t{tlen}\t{qlen}\t{tcov}\t{qcov}\t80\t1e-15')
                self.assertEqual(a,b);self.assertEqual(passing(a,SETTINGS),passing(b,SETTINGS))
                self.assertEqual(shorter_protein_coverage(qlen,tlen,qcov,tcov),shorter_protein_coverage(tlen,qlen,tcov,qcov))
    def test_coverage_threshold_boundaries(self):
        for coverage,expected in ((.699999,False),(math.nextafter(.7,0.),False),(.7,True),(.700001,True),(1.,True)):
            with self.subTest(coverage=coverage):
                self.assertEqual(passing(hit('a','b',length=150,qcov=coverage,tcov=coverage),SETTINGS),expected)
    def test_sequence_coverage_bounds_enforced(self):
        for value in (-.01,1.01,float('nan'),float('inf'),True):
            for qcov,tcov in ((value,.8),(.8,value)):
                with self.subTest(qcov=qcov,tcov=tcov):
                    with self.assertRaisesRegex(ValueError,'between 0 and 1'):shorter_protein_coverage(100,200,qcov,tcov)
    def test_actual_search_persists_bounded_sequence_coverages(self):
        cs=[{'cds':[{'id':'a','protein':'MACDEFGHIKLMNPQ'},{'id':'b','protein':'MACDEFGHIKLMNPQ'}]}]
        hits,_=search(cs,threads=2)
        self.assertTrue(hits)
        for h in hits:
            self.assertEqual(h['shorter_coverage'],min(h['query_coverage'],h['target_coverage']))
            self.assertTrue(0<=h['shorter_coverage']<=1)
    def test_abc_transitivity(self):
        p=simple_project();self.assertEqual(p['clusters'][0]['cds'][0]['group_id'],p['clusters'][2]['cds'][0]['group_id']);self.assertFalse(any(h['query']=='a0' and h['target']=='c0' for h in p['pairwise_hits']))
    def test_duplications(self):
        cs=[cluster('a',[('x',100,1),('x',500,1)])];groups=connected_components(cs,[hit('a0','a1')],SETTINGS);self.assertEqual(len(groups),1)
    def test_unrelated(self):self.assertEqual(len(connected_components([cluster('a',[('x',100,1),('y',500,1)])],[],SETTINGS)),2)
    def test_real_mmseqs_identical_and_unrelated(self):
        import random
        rng=random.Random(900);seq='M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=250))
        cs=[cluster('a',[('x',100,1),('y',500,1),('x',900,1)])]
        cs[0]['cds'][0]['protein']=seq;cs[0]['cds'][2]['protein']=seq
        cs[0]['cds'][1]['protein']='M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=250))
        for g in cs[0]['cds']:g['protein_length']=len(g['protein'])
        hits,meta=search(cs,threads=2)
        self.assertTrue(any(h['query']=='a0' and h['target']=='a2' and h['identity']==1 and passing(h,SETTINGS) for h in hits))
        self.assertFalse(any('a1' in (h['query'],h['target']) and passing(h,SETTINGS) for h in hits));self.assertEqual(meta['max_seqs'],3)
    def test_short_partial_proteins_seedless_search(self):
        cs=[{'cds':[{'id':'a','protein':'MACDEFGHIKLMNPQ'},{'id':'b','protein':'MACDEFGHIKLMNPQ'}]}]
        hits,meta=search(cs,threads=2);self.assertEqual(meta['prefilter_mode'],2);self.assertTrue(any(passing(h,SETTINGS)for h in hits))
    def test_mixed_short_protein_supplement(self):
        cs=[{'cds':[{'id':'a','protein':'MACDEFGHIKLMNPQ'},{'id':'b','protein':'MACDEFGHIKLMNPQ'},{'id':'c','protein':'MACDEFGHIKLMNPQ'+'ACDEFGHIKLMNPQ'*20}]}]
        hits,meta=search(cs,threads=2);self.assertTrue(meta['short_protein_supplement']);self.assertTrue(any(h['query']=='a'and h['target']=='b'for h in hits))

class LayoutTests(unittest.TestCase):
    def test_deterministic_order_orientation_offset(self):
        cs=simple_project()['clusters'];self.assertEqual(automatic_layout(cs),automatic_layout(list(reversed(cs))))
    def test_reversal_and_alignment(self):
        cs=[cluster('a',[('x',1000,1),('y',5000,1)]),cluster('b',[('y',4700,-1),('x',8700,-1)])]
        layout=automatic_layout(cs);self.assertEqual(layout['rows']['a']['orientation'],1);self.assertEqual(layout['rows']['b']['orientation'],-1)
        self.assertEqual(layout['rows']['b']['offset'],0)
    def test_offsets(self):
        cs=[cluster('a',[('x',1000,1),('y',5000,1)]),cluster('b',[('x',3000,1),('y',7000,1)])]
        self.assertEqual(automatic_layout(cs)['rows']['b']['offset'],-2000)
    def test_coordinate_preservation(self):
        cs=simple_project()['clusters'];before=copy.deepcopy(cs);automatic_layout(cs);self.assertEqual(cs,before)
    def test_no_homolog_align_left(self):
        cs=[cluster('a',[('x',1000,1)]),cluster('b',[('y',3000,-1)])];l=automatic_layout(cs)
        self.assertEqual(l['rows']['a']['offset'],l['rows']['b']['offset']);self.assertEqual(l['rows']['b']['orientation'],1)
    def test_preserve(self):
        cs=list(reversed(simple_project()['clusters']));l=automatic_layout(cs,True);self.assertEqual(l['order'],[c['id']for c in cs]);self.assertTrue(all(r=={'orientation':1,'offset':0.0}for r in l['rows'].values()))

class PersistenceTests(unittest.TestCase):
    def test_presentation_scales_roundtrip_without_biological_changes(self):
        p=simple_project();before=copy.deepcopy(p)
        p['visualization'].update(width_scale=.65,height_scale=1.75,row_spacing_scale=1.5)
        self.assertEqual(from_html(to_html(p)),p)
        for key in ('clusters','layout','homolog_groups','pairwise_hits'):
            self.assertEqual(p[key],before[key])
    def test_optional_scale_defaults_do_not_mutate_project(self):
        p=simple_project();self.assertEqual(validate(p),p)
        for key in ('width_scale','height_scale','row_spacing_scale'):
            self.assertNotIn(key,p['visualization'])
    def test_presentation_scales_are_bounded_finite_numbers(self):
        for key in ('width_scale','height_scale','row_spacing_scale'):
            for value in (0,.249,4.001,float('nan'),float('inf'),True,'1',None):
                with self.subTest(key=key,value=value):
                    p=simple_project();p['visualization'][key]=value
                    with self.assertRaisesRegex(ValueError,'between 0.25 and 4'):validate(p)
            for value in (.25,1,4):
                p=simple_project();p['visualization'][key]=value;self.assertEqual(validate(p),p)
    def test_model_rejects_unbounded_or_nonfinite_shorter_coverage(self):
        for coverage in (-.01,1.01,float('nan'),float('inf'),True):
            with self.subTest(coverage=coverage):
                p=simple_project();p['pairwise_hits'][0]['shorter_coverage']=coverage
                with self.assertRaisesRegex(ValueError,'between 0 and 1'):to_html(p)
    def test_model_validates_sequence_coverage_derivation(self):
        p=simple_project();p['pairwise_hits'][0]['shorter_coverage']=.8
        with self.assertRaisesRegex(ValueError,'Inconsistent shorter-protein coverage'):validate(p)
    def test_gapped_alignment_roundtrip_uses_sequence_coverage(self):
        p=simple_project();p['pairwise_hits'][0]['alignment_length']=150
        self.assertEqual(from_html(to_html(p)),p)
        self.assertEqual(p['pairwise_hits'][0]['shorter_coverage'],.9)
    def test_complete_html_roundtrip(self):
        p=simple_project();p['layout']['order'].reverse();p['layout']['rows']['a']={'orientation':-1,'offset':-1234.5}
        p['color_overrides']['genes']['a0']='#123456';p['color_overrides']['groups'][p['homolog_groups'][0]['id']]='#abcdef'
        p['visualization']={'connectors':False,'connector_identity':.75,'labels':True};p['fasta_selection']=['a0']
        self.assertEqual(p,from_html(to_html(p)))
    def test_safe_annotation_serialization(self):
        p=simple_project();p['clusters'][0]['cds'][0]['product']='</script><img src=x onerror=alert(1)>& α'
        html=to_html(p);self.assertNotIn('</script><img',html);self.assertEqual(from_html(html),p)
    def test_schema_version_rejected(self):
        p=simple_project();p['schema_version']=2
        with self.assertRaisesRegex(ValueError,'Unsupported'):to_html(p)
    def test_marker_required(self):
        with self.assertRaisesRegex(ValueError,'marker'):from_html('<html>hello</html>')
    def test_duplicate_marker_rejected(self):
        p=simple_project();html=to_html(p);marker='<script id="kstack-project" type="application/json">{}</script>'
        with self.assertRaisesRegex(ValueError,'one'):from_html(html+marker)
    def test_validation_rejects_bad_coordinates(self):
        p=simple_project();p['clusters'][0]['cds'][0]['segments'][0]['end']=20000
        with self.assertRaisesRegex(ValueError,'coordinates'):validate(p)
    def test_validation_rejects_invalid_order(self):
        p=simple_project();p['layout']['order']=['a','a','b']
        with self.assertRaisesRegex(ValueError,'permutation'):validate(p)
    def test_validation_rejects_nucleotide(self):
        p=simple_project();p['clusters'][0]['sequence']='ACGT'
        with self.assertRaisesRegex(ValueError,'nucleotide'):validate(p)
    def test_real_reopen_recompute_from_proteins(self):
        c=parse_genbank(ROOT/'fixtures/realistic/A. alliaceus CBS 536.65.gbk');p=create_project([c]);restored=from_html(to_html(p))
        r=recompute(restored,{'identity':.4,'coverage':.8});self.assertEqual(r['analysis']['identity'],.4);self.assertEqual(r['layout'],p['layout']);self.assertNotIn('sequence',r['clusters'][0])
    def test_group_override_survival(self):
        from unittest.mock import patch
        p=simple_project();p['color_overrides']['groups'][p['homolog_groups'][0]['id']]='#123456'
        with patch('kstack.model.search',return_value=(p['pairwise_hits'],{})):
            self.assertEqual(recompute(p,SETTINGS)['color_overrides'],p['color_overrides'])
    def test_cli_json_interface_and_circular_choice(self):
        payload={'action':'inspect','paths':[str(ROOT/'fixtures/synthetic/edge-cases.gbk')]}
        r=subprocess.run([sys.executable,'-m','kstack.cli'],input=json.dumps(payload),capture_output=True,text=True)
        self.assertEqual(r.returncode,0);self.assertTrue(json.loads(r.stdout)['ok'])
        with self.assertRaisesRegex(ValueError,'Circular'):handle({'action':'generate','inputs':[{'path':str(ROOT/'fixtures/synthetic/circular-origin.gbk')}]})
    def test_circular_linearize_and_omit_generation(self):
        circular=str(ROOT/'fixtures/synthetic/circular-origin.gbk');linear=str(ROOT/'fixtures/realistic/A. alliaceus CBS 536.65.gbk')
        result=handle({'action':'generate','inputs':[{'path':circular,'circular':'linearize'}]})
        c=result['project']['clusters'][0];self.assertTrue(c['linearized_at_origin']);self.assertEqual(c['cds'][0]['protein'],'M'+'V'*199)
        self.assertEqual(c['cds'][0]['segments'][0]['start'],2700)
        omitted=handle({'action':'generate','inputs':[{'path':circular,'circular':'omit'},{'path':linear}]})
        self.assertEqual(len(omitted['project']['clusters']),1);self.assertEqual(omitted['project']['clusters'][0]['topology'],'linear')
    def test_inspected_circular_requires_explicit_decision(self):
        with self.assertRaisesRegex(ValueError,'Circular'):handle({'action':'generate','clusters':[parse_genbank(ROOT/'fixtures/synthetic/circular-origin.gbk')]})

class ColorTests(unittest.TestCase):
    def test_deterministic_palette_and_singletons(self):
        a=[{'id':f'h{i:03}','members':['a','b'] if i<100 else ['c']}for i in range(101)];b=copy.deepcopy(a);generate_colors(a);generate_colors(b)
        self.assertEqual(a,b);self.assertEqual(len({g['generated_color']for g in a[:100]}),100);self.assertEqual(a[-1]['generated_color'],'#a3aab4')
    def test_palette_does_not_cycle_after_288_groups(self):
        groups=[{'id':f'h{i:04}','members':['a','b']}for i in range(350)];generate_colors(groups)
        self.assertEqual(len({g['generated_color']for g in groups}),350)

if __name__=='__main__':unittest.main()
