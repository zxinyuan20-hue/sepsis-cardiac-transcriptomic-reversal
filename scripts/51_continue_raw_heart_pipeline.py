"""Continue this local execution after active module48; never a recurring job.
Inputs: completed48status; outputs: execution_state.json and existing49/50 outputs.
Dependencies: standard library plus shared pipeline. Fixed scientific configs.
Stops on failure or QC hold; no unapproved threshold changes or sample deletion.
"""
import json,logging,subprocess,sys,time
from pipeline_utils import ROOT,write_json,run_standard_module,now
OUT=ROOT/'outputs/analysis/raw_heart_v1'
def main():
    state=OUT/'execution_state.json';deadline=time.monotonic()+18*3600
    status=ROOT/'outputs/audit/48_quantify_raw_heart_status.json'
    write_json(state,{'phase':'WAITING_FOR_ACTIVE_QUANTIFICATION','updated_utc':now(),'disease_inference':'NOT_RUN','release':'PREVIOUS_RESULTS_REMAIN_CURRENT'})
    while not status.exists():
        if time.monotonic()>deadline:raise TimeoutError('Module48 did not finish within18hours')
        time.sleep(30)
    result=json.loads(status.read_text(encoding='utf8'))
    if result['status']!='PASS':
        write_json(state,{'phase':'STOPPED_QUANTIFICATION_FAILURE','module48':result,'updated_utc':now()})
        raise RuntimeError('Module48 failed; preserved data and logs')
    qa=json.loads((OUT/'quantification_audit.json').read_text())
    tasks=['49_analyze_reconstructed_heart.py','50_audit_reconstructed_heart.py'] if qa['disease_inference_ready'] else ['50_audit_reconstructed_heart.py']
    for task in tasks:
        write_json(state,{'phase':'RUNNING_'+task,'updated_utc':now(),'qc_hold_samples':qa['qc_hold_samples']})
        r=subprocess.run([sys.executable,str(ROOT/'scripts'/task)],cwd=ROOT)
        if r.returncode:
            write_json(state,{'phase':'STOPPED_FAILURE','failed_task':task,'returncode':r.returncode,'updated_utc':now()})
            raise RuntimeError(f'{task} failed; no further task executed')
    write_json(state,{'phase':'FINISHED_WITH_QC_HOLD' if qa['qc_hold_samples'] else 'FINISHED_AUDITED','updated_utc':now(),'qc_hold_samples':qa['qc_hold_samples'],'report':'outputs/audit/完整人类心肌重建与疾病复现审计.md'})
    return json.loads(state.read_text(encoding='utf8'))
if __name__=='__main__':run_standard_module('51_continue_raw_heart_pipeline',main)
