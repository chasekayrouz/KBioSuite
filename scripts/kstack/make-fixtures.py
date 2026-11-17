"""Deterministic synthetic GenBanks; no private inputs. Also used by tests/benchmarks."""
from pathlib import Path
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, BeforePosition, AfterPosition
from Bio import SeqIO

CODONS = dict(zip('ACDEFGHIKLMNPQRSTVWY', ['GCT','TGT','GAT','GAA','TTT','GGT','CAT','ATT','AAA','CTG','ATG','AAT','CCT','CAA','CGT','TCT','ACT','GTT','TGG','TAT']))

def edge_record(circular=False):
    dna=list('A'*3600)
    def put(start,nt): dna[start:start+len(nt)]=nt
    translated='ATG'+'GCT'*99+'TAA'
    put(90,translated)
    minus='ATG'+'TGT'*99+'TAA';put(600,str(Seq(minus).reverse_complement()))
    compound='ATG'+'GGT'*99+'TAA';put(1200,compound[:150]);put(1500,compound[150:])
    partial='ATG'+'CCT'*99+'GC';put(1800,partial)
    rec=SeqRecord(Seq(''.join(dna)),id='SYNTHETIC',name='SYNTHETIC',description='Kstack exact edge cases')
    rec.annotations={'molecule_type':'DNA','topology':'circular' if circular else 'linear'}
    rec.features=[
        SeqFeature(SimpleLocation(90,393,1),type='CDS',qualifiers={'gene':['normal'],'locus_tag':['SYN_001'],'product':['hypothetical protein'],'translation':['M'+'A'*99]}),
        SeqFeature(SimpleLocation(600,903,-1),type='CDS',qualifiers={'gene':['minus'],'locus_tag':['SYN_002'],'transl_table':['11']}),
        SeqFeature(CompoundLocation([SimpleLocation(1200,1350,1),SimpleLocation(1500,1653,1)]),type='CDS',qualifiers={'gene':['compound']}),
        SeqFeature(SimpleLocation(BeforePosition(1800),AfterPosition(2102),1),type='CDS',qualifiers={'gene':['partial']}),
        SeqFeature(SimpleLocation(2400,2500,1),type='CDS',qualifiers={'gene':['malformed']}),
        SeqFeature(SimpleLocation(2700,2750,1),type='tRNA',qualifiers={'product':['not shown']}),
    ]
    return rec

def origin_record(strand=1):
    protein='M'+'V'*199;nt='ATG'+'GTT'*199+'TAA'
    dna=list('A'*3000)
    if strand==1:
        dna[2700:3000]=nt[:300];dna[:303]=nt[300:]
        parts=[SimpleLocation(2700,3000,1),SimpleLocation(0,303,1)]
    else:
        rev=str(Seq(nt).reverse_complement());dna[:303]=rev[300:];dna[2700:3000]=rev[:300]
        # Biological extraction order on minus strand is reversed.
        parts=[SimpleLocation(0,303,-1),SimpleLocation(2700,3000,-1)]
    r=SeqRecord(Seq(''.join(dna)),id='ORIGIN',name='ORIGIN',description='Origin-spanning CDS')
    r.annotations={'molecule_type':'DNA','topology':'circular'}
    r.features=[SeqFeature(CompoundLocation(parts),type='CDS',qualifiers={'gene':['origin'],'transl_table':['11']})]
    return r

def bgc_records(count=20, genes=30, seed=717):
    rng=random.Random(seed);alphabet=list(CODONS)
    families=['M'+''.join(rng.choices(alphabet,k=length-1)) for length in [240+(i%5)*45 for i in range(genes)]]
    records=[]
    for c in range(count):
        dna=list('A'*44000);features=[]
        order=list(range(genes))
        if c%4==1:order.reverse()
        for i,family in enumerate(order):
            sequence=list(families[family])
            # Independent ~20% amino-acid substitutions per cluster; shared families
            # include both long proteins and duplicates, with real intergenic spacing.
            for pos in range(1,len(sequence)):
                if rng.random()<.20:sequence[pos]=rng.choice(alphabet)
            protein=''.join(sequence);nt=''.join(CODONS[a] for a in protein)+'TAA'
            start=200+i*1400;strand=1 if (i+c)%3 else -1
            dna[start:start+len(nt)]=nt if strand==1 else str(Seq(nt).reverse_complement())
            features.append(SeqFeature(SimpleLocation(start,start+len(nt),strand),type='CDS',qualifiers={
                'gene':[f'bgc{family:02d}'],'locus_tag':[f'BGC{c+1:02d}_{i+1:03d}'],
                'product':['biosynthetic enzyme' if i%4 else 'hypothetical protein'],'translation':[protein]}))
        record=SeqRecord(Seq(''.join(dna)),id=f'BGC{c+1:02d}',name=f'BGC{c+1:02d}',description='Synthetic BGC-like benchmark with homologous protein families')
        record.annotations={'molecule_type':'DNA','topology':'linear'};record.features=features;records.append(record)
    return records

if __name__=='__main__':
    root=Path(__file__).resolve().parents[2]/'tests/kstack/fixtures/synthetic';root.mkdir(parents=True,exist_ok=True)
    SeqIO.write(edge_record(),root/'edge-cases.gbk','genbank')
    text=(root/'edge-cases.gbk').read_text().replace('2401..2500','nonsense')
    (root/'edge-cases.gbk').write_text(text)
    SeqIO.write(origin_record(),root/'circular-origin.gbk','genbank')
    SeqIO.write(origin_record(-1),root/'circular-origin-minus.gbk','genbank')
