"""Packages: shared pipeline utilities and standard library.
Inputs/outputs: config/robustness_amendment_v1.json; ROOT-relative inputs.
Reproducibility: record seed, file hashes and post-result timing before new scores.
Pipeline: check historical manifests, freeze inputs, register dated amendment.
Saved files: plan_freeze.json and amendment.md; shared module execution log.
"""
import json
from pipeline_utils import ROOT,run_standard_module,write_json,sha256,now

def main():
 cp=ROOT/'config/robustness_amendment_v1.json';cfg=json.loads(cp.read_text(encoding='utf8'));out=ROOT/cfg['output'];out.mkdir(parents=True,exist_ok=True)
 checks=[]
 for name in cfg['historical_manifests']:
  data=json.loads((ROOT/name).read_text(encoding='utf8'))
  rows=data if isinstance(data,list) else data.get('files',[])
  assert rows,name
  drift=[]
  for row in rows:
   actual=sha256(ROOT/row['path'])
   if actual!=row['sha256']:
    assert row['path'] in ['scripts/19_release_screening.py','docs/03_后续分析路线与任务表.md','docs/04_协议修订记录.md'],row['path']
    drift.append({'path':row['path'],'historical_sha256':row['sha256'],'current_sha256':actual,'status':'pre-existing non-input code/document drift; not repaired or silently re-frozen as historical'})
  checks.append({'manifest':name,'files':len(rows),'preexisting_drift':drift,'all_scientific_data_and_results_match':True})
 files=[cp,ROOT/'CURRENT_RELEASE.json',ROOT/cfg['source_manuscript']]
 for folder in ['outputs/analysis/screening','outputs/analysis/expansion','outputs/analysis/external_compound_v1']:
  files.extend(p for p in (ROOT/folder).glob('*') if p.is_file())
 files.extend([ROOT/'outputs/preparation/lincs_compound_dictionary.tsv',ROOT/'data/reference/mechanism/PMC11390749_fulltext.xml'])
 freeze={'recorded_utc':now(),'version':cfg['version'],'prior_knowledge':'482-entry zero hits; cardiac 21-drug zero hits; two duplicated source-key groups known before this extension','historical_checks':checks,'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha256(p)} for p in files]}
 target=out/'plan_freeze.json'
 if target.exists():
  prev=json.loads(target.read_text(encoding='utf8'));assert prev['files']==freeze['files'];return {'status':'existing freeze verified'}
 write_json(target,freeze)
 (out/'amendment.md').write_text('# 有限稳健性扩展方案\n\n记录时间：'+freeze['recorded_utc']+'\n\n本次为已知历史结果后的探索性扩展，不是前瞻性注册。已知482条目和21药物未获得正式支持，且来源字典存在两组完整结构键重复。\n\n'+ '\n\n'.join('## '+k+'\n\n'+str(v) for k,v in cfg.items() if k not in ['historical_manifests'])+'\n',encoding='utf8')
 return {'frozen_files':len(files),'historical_manifests':checks}

if __name__=='__main__':run_standard_module('77_freeze_robustness_amendment',main)
