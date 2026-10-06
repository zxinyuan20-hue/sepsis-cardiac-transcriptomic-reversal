"""Create ignored local settings. Standard library; no network/data operations."""
import json,os,shutil,sys
from pathlib import Path
root=Path(__file__).resolve().parent
p=root/'config/runtime.local.json'
if p.exists():
    print('Existing local runtime preserved.')
else:
    d=json.loads((root/'config/runtime.example.json').read_text())
    d['python']=sys.executable
    d['rscript']=os.environ.get('RSCRIPT',shutil.which('Rscript') or 'Rscript')
    d['local_rat_source']=os.environ.get('H9C2_SOURCE_DIR',str(root/'data/private/H9c2'))
    p.write_text(json.dumps(d,indent=2),encoding='utf8')
    print('Created ignored config/runtime.local.json. Check Rscript before running R modules.')
