"""Offline source/integrity checks. No study data required or downloaded."""
import ast,hashlib,json,sys
from pathlib import Path
root=Path(__file__).resolve().parent
errors=[]
for p in (root/'scripts').glob('*.py'):
    try: ast.parse(p.read_text(encoding='utf8'))
    except Exception as e: errors.append(str(p.relative_to(root))+': '+str(e))
for p in (root/'config').glob('*.json'):
    try: json.loads(p.read_text(encoding='utf8'))
    except Exception as e: errors.append(str(p.relative_to(root))+': '+str(e))
manifest=root/'MANIFEST_SHA256.json'
if manifest.exists():
    for row in json.loads(manifest.read_text())['files']:
        p=root/row['path']
        if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']: errors.append(row['path'])
print(json.dumps({'status':'FAIL' if errors else 'PASS','errors':errors,'scope':'source syntax, JSON and file integrity only; no biological rerun'},indent=2))
sys.exit(bool(errors))
