"""Packages: pandas/numpy; local source path in runtime config, membership in local_rat.json.
Outputs: outputs/preprocessing/rat/source_files and audit/provenance, never modifies originals.
Seed: shared configuration. Import -> sample contracts -> count QC; no group-effect tests.
"""
import hashlib,json,tarfile
from pathlib import Path
import numpy as np
import pandas as pd
from pipeline_utils import ROOT,RUNTIME,path,sha256,write_json,now,run_standard_module

def main():
    cfg=json.loads((ROOT/'config/local_rat.json').read_text(encoding='utf-8'))
    source=Path(RUNTIME[cfg['source_location_key']])
    out=ROOT/'outputs/preprocessing/rat';raw=out/'source_files';raw.mkdir(parents=True,exist_ok=True)
    entries=[]
    for archive,desired in [('Summary.tar.gz',cfg['summary_members']),('Report.tar.gz',cfg['report_members'])]:
        archive_path=source/archive;archive_hash=sha256(archive_path);seen=[]
        with tarfile.open(archive_path,'r|gz') as t:
            for member in t:
                if member.name not in desired:continue
                if not member.isfile():raise ValueError('Expected regular member')
                payload=t.extractfile(member).read();digest=hashlib.sha256(payload).hexdigest()
                dest=raw/Path(member.name).name
                if dest.exists():
                    if sha256(dest)!=digest:raise ValueError('Existing imported source differs')
                else:dest.write_bytes(payload)
                entries.append({'source_archive':str(archive_path),'archive_sha256':archive_hash,'member':member.name,'bytes':len(payload),'sha256':digest,'extracted_to':dest.relative_to(ROOT).as_posix(),'retrieved_utc':now(),'access':'private; local only'})
                seen.append(member.name)
        if set(seen)!=set(desired):raise ValueError('Missing requested archive members')
    data=pd.read_csv(raw/'gene_count_matrix.txt',sep='\t',index_col=0)
    expected=[s for group in cfg['groups'].values() for s in group]
    if set(data.columns)!=set(expected) or data.index.duplicated().any():raise ValueError('Sample/gene contract mismatch')
    data=data[expected];a=data.to_numpy(dtype=float)
    if not np.isfinite(a).all() or (a<0).any():raise ValueError('Invalid counts')
    sample_table=pd.read_csv(raw/'sample_info.txt',sep='\t').rename(columns={'#SampleID':'sample_id','COND1':'group'})
    known={s:g for g,ss in cfg['groups'].items() for s in ss}
    if any(known.get(row.sample_id)!=row.group for row in sample_table.itertuples()):raise ValueError('User/source group mismatch')
    sample_table['independent_biological_sample']=True
    sample_table['material']=cfg['material'];sample_table['batch']='not_reported'
    sample_table.to_csv(out/'sample_manifest.tsv',sep='\t',index=False)
    duplicates=[]
    for i,left in enumerate(expected):
        for right in expected[i+1:]:
            if np.array_equal(data[left].values,data[right].values):duplicates.append([left,right])
    if duplicates:raise ValueError(f'Exact duplicate samples: {duplicates}')
    data.to_csv(out/'rat_counts.tsv.gz',sep='\t',compression='gzip')
    qc={'shape':list(data.shape),'groups':{g:len(ss) for g,ss in cfg['groups'].items()},'nonfinite':int((~np.isfinite(a)).sum()),'negative':int((a<0).sum()),'fractional_elements':int((a!=np.floor(a)).sum()),'exact_duplicate_columns':duplicates,'library_sums':data.sum().to_dict(),'zero_genes':int((data.sum(axis=1)==0).sum()),'limitations':['StringTie count derivation not yet documented','LPS dose/duration and exact cell origin/batches pending'],'scope':'input QC only; held-out effects remain unread'}
    write_json(out/'source_manifest.json',entries);write_json(out/'input_qc.json',qc)
    return qc
if __name__=='__main__':run_standard_module('08_import_local_rat',main)
