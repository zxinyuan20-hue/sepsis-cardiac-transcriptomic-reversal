"""Dependencies: cutadapt 5.2, Subread2.1.1, pandas. Config centralized.
Inputs: verified14FASTQ+GENCODEv47; outputs: trimmedFASTQ,BAM,counts,QC/logs.
Reproducibility: frozen parameters, deterministic sample order, stage hashes.
No disease effects computed; all strand modes retained for technical audit.
"""
import gzip,json,logging,subprocess,sys,collections
from pathlib import Path
import pandas as pd
from pipeline_utils import ROOT,sha256,write_json,run_standard_module,now
QPATH=ROOT/'config/raw_heart_quantification_v1.json'
Q=json.loads(QPATH.read_text(encoding='utf8'));OUT=ROOT/Q['out']

def raw_quality_audit(fq,run,expected_reads,expected_bases):
    dst=OUT/'raw_qc'/f'{run}.json';dst.parent.mkdir(exist_ok=True)
    if dst.exists():
        old=json.loads(dst.read_text());assert old['input_sha256']==sha256(fq);return old
    n=0;bases=0;qual=collections.Counter();lengths=collections.Counter()
    with gzip.open(fq,'rb') as f:
        while True:
            head=f.readline()
            if not head:break
            seq=f.readline().rstrip(b'\r\n');plus=f.readline();q=f.readline().rstrip(b'\r\n')
            assert head.startswith(b'@') and plus.startswith(b'+') and len(seq)==len(q),(run,n)
            n+=1;bases+=len(seq);lengths[len(seq)]+=1
            if q==b'?'*len(q):qual[63]+=len(q)
            else:qual.update(q)
    assert n==expected_reads and bases==expected_bases,(run,n,bases,expected_reads,expected_bases)
    result={'run':run,'reads':n,'bases':bases,'lengths':dict(lengths),'phred33_histogram':{str(k-33):v for k,v in qual.items()},'constant_quality':len(qual)==1,'quality_interpretation':'Constant encoded scores cannot establish original instrument base quality' if len(qual)==1 else 'Encoded archive qualities; instrument provenance not independently verified','input_sha256':sha256(fq),'gzip_crc_and_records_pass':True}
    write_json(dst,result);return result

def stage(name,args,outputs,inputs):
    meta=OUT/'stages'/f'{name}.json';meta.parent.mkdir(exist_ok=True)
    inputs={str(p.relative_to(ROOT)):sha256(p) for p in inputs}
    signature={'command':[str(a) for a in args],'inputs':inputs,'config_sha256':sha256(QPATH)}
    if meta.exists():
        old=json.loads(meta.read_text(encoding='utf8'))
        assert old['signature']==signature, f'Stage signature changed: {name}'
        assert all((ROOT/p).exists() and sha256(ROOT/p)==h for p,h in old['outputs'].items()),name
        logging.info('REUSE_STAGE %s',name);return
    for p in outputs:
        if p.exists():raise RuntimeError(f'Unverified existing derived output preserved: {p}')
    log=OUT/'logs'/f'{name}.log';log.parent.mkdir(exist_ok=True)
    logging.info('RUN %s',name)
    with log.open('wb') as f:
        r=subprocess.run([str(a) for a in args],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    assert r.returncode==0, f'{name} exit={r.returncode}; inspect {log}'
    assert all(p.exists() and p.stat().st_size>0 for p in outputs),name
    if name=='build_index':outputs=list((OUT/'index').glob('GRCh38_v47.*'))
    write_json(meta,{'signature':signature,'outputs':{p.relative_to(ROOT).as_posix():sha256(p) for p in outputs},'finished_utc':now()})

def main():
    runtime=json.loads((OUT/'runtime_reference_audit.json').read_text(encoding='utf8'))
    raw=json.loads((OUT/'raw_acquisition_audit.json').read_text(encoding='utf8'))
    assert raw['verified']==14 and not raw['failures']
    import cutadapt
    assert cutadapt.__version__==Q['cutadapt_version']
    freeze=OUT/'quantification_plan_freeze.json'
    sig={'config_sha256':sha256(QPATH),'runtime_audit_sha256':sha256(OUT/'runtime_reference_audit.json'),'fastq_manifest_sha256':sha256(OUT/'frozen_fastq_manifest.tsv')}
    inherited=[ROOT/'outputs/analysis/human/primary_sepsis_vs_nonfailing.tsv',ROOT/'outputs/analysis/external_human/external_human_summary.json',ROOT/'outputs/analysis/external_human/GSE141864_fixed_pathway_tests.tsv',ROOT/'outputs/analysis/mechanism/fixed_pathways.tsv']
    inherited += [ROOT/f'outputs/analysis/human/frozen_signatures/landmarks_{size}_{direction}.tsv' for size in [50,100] for direction in ['up','down']]
    sig['inherited_analysis_inputs']={p.relative_to(ROOT).as_posix():sha256(p) for p in inherited}
    if freeze.exists():assert json.loads(freeze.read_text())['signature']==sig
    else:write_json(freeze,{'signature':sig,'parameters':Q,'frozen_utc':now(),'disease_effects_computed':False})
    for d in ['index','trimmed','bam','counts','logs','stages']: (OUT/d).mkdir(exist_ok=True)
    ref=OUT/'reference';fasta=ref/'GRCh38.primary_assembly.genome.fa';gtf=ref/'gencode.v47.primary_assembly.annotation.gtf'
    for name,h in runtime['reference_sha256'].items():assert sha256(ref/name)==h,name
    tools={k:ROOT/v['path'] for k,v in runtime['tools'].items()}
    for k,p in tools.items():assert sha256(p)==runtime['tools'][k]['sha256']
    chroms=set()
    with fasta.open() as f:
        for line in f:
            if line.startswith('>'):chroms.add(line[1:].split()[0])
    anno=set()
    with gtf.open() as f:
        for line in f:
            if not line.startswith('#'):anno.add(line.split('\t')[0])
    assert anno<=chroms,anno-chroms
    idx=OUT/'index'/'GRCh38_v47'
    # Sparse index can have multiple blocks under the configured memory limit.
    idxfiles=[Path(str(idx)+s) for s in ['.00.b.array','.00.b.tab','.files','.reads']]
    stage('build_index',[tools['subread-buildindex.exe'],'-M',str(Q['index_memory_mb']),'-o',idx,fasta],idxfiles,[fasta,tools['subread-buildindex.exe']])
    idxfiles=sorted((OUT/'index').glob('GRCh38_v47.*'))
    m=pd.read_csv(OUT/'frozen_fastq_manifest.tsv',sep='\t');rows=[]
    ena=pd.read_csv(OUT/'sources/ena_submission_format.tsv',sep='\t').set_index('run_accession')
    for row in m.itertuples():
        run=row.run_accession;fq=OUT/'sources/fastq'/f'{run}.fastq.gz'
        check=json.loads((OUT/'checksums'/f'{run}.json').read_text())
        assert sha256(fq)==check['sha256']
        rawqc=raw_quality_audit(fq,run,int(ena.loc[run,'read_count']),int(ena.loc[run,'base_count']))
        trimmed=OUT/'trimmed'/f'{run}.fastq.gz';tj=OUT/'trimmed'/f'{run}.cutadapt.json'
        stage(run+'_trim',[sys.executable,'-m','cutadapt',*Q['trim_arguments'],'--json',tj,'-o',trimmed,fq],[trimmed,tj],[fq])
        bam=OUT/'bam'/f'{run}.bam'
        stage(run+'_align',[tools['subjunc.exe'],*Q['subjunc_arguments'],'-a',gtf,'-i',idx,'-r',trimmed,'-o',bam],[bam],[trimmed,gtf,tools['subjunc.exe'],*idxfiles])
        assigned={};totals={}
        for strand in Q['audit_strands']:
            dst=OUT/'counts'/f'{run}.s{strand}.tsv';summary=Path(str(dst)+'.summary')
            stage(run+f'_count_s{strand}',[tools['featureCounts.exe'],*Q['featurecounts_arguments'],'-s',str(strand),'-a',gtf,'-o',dst,bam],[dst,summary],[bam,gtf,tools['featureCounts.exe']])
            s=pd.read_csv(summary,sep='\t',index_col=0).iloc[:,0];assigned[strand]=int(s['Assigned']);totals[strand]=int(s.sum())
        t=json.loads(tj.read_text());n=t['read_counts']['input'];kept=t['read_counts']['output'];assert n==rawqc['reads']
        assert all(totals[x]==kept for x in totals),(run,totals,kept)
        metrics={'title':row.title,'run':run,'raw_reads':n,'retained_reads':kept,'retained_read_fraction':kept/n,'assigned_s0':assigned[0],'assigned_s1':assigned[1],'assigned_s2':assigned[2],'assigned_forward_fraction':assigned[1]/kept,'assigned_forward_reads':assigned[1],'forward_share_of_forward_plus_reverse':assigned[1]/max(1,assigned[1]+assigned[2])}
        holds=[k for k,v in Q['qc_hold_thresholds'].items() if metrics[k.removesuffix('_min')]<v]
        metrics['qc_hold']=';'.join(holds);metrics['constant_archive_quality']=rawqc['constant_quality'];rows.append(metrics)
        pd.DataFrame(rows).to_csv(OUT/'technical_qc.tsv',sep='\t',index=False)
        logging.info('SAMPLE_COUNTED %s assigned_forward=%s hold=%s',run,assigned[1],holds)
    matrix=None
    for row in m.itertuples():
        t=pd.read_csv(OUT/'counts'/f'{row.run_accession}.s1.tsv',sep='\t',comment='#')
        z=t.set_index('Geneid').iloc[:,-1].rename(row.title)
        matrix=z.to_frame() if matrix is None else matrix.join(z,how='outer',validate='one_to_one')
    assert matrix.shape[1]==14 and not matrix.isna().any().any()
    matrix.to_csv(OUT/'gene_counts_14.tsv.gz',sep='\t',compression='gzip')
    holds=[r['title'] for r in rows if r['qc_hold']]
    report={'samples':14,'genes':len(matrix),'counts_sha256':sha256(OUT/'gene_counts_14.tsv.gz'),'qc_hold_samples':holds,'disease_inference_ready':not holds,'all_samples_retained':True}
    write_json(OUT/'quantification_audit.json',report)
    return report
if __name__=='__main__':run_standard_module('48_quantify_raw_heart',main)
