"""Package final outputs and source evidence; validate every archived file."""
from pathlib import Path
import zipfile,json,hashlib
P=Path(__file__).resolve().parents[1];OUT=P/'Research_package_20261117.zip'
def include(f):
 r=f.relative_to(P);parts=r.parts
 if f.match('Research_package_*.zip') or f.name in ['MANIFEST.json','PACKAGE_SHA256.txt'] or f.name.startswith('~$'):return False
 if any(v in parts for v in ['__pycache__','.git','vendor']):return False
 if f.suffix.lower() in ['.xpt','.tmp','.log','.pyc']:return False
 if parts[0]=='qa' and len(parts)>1 and parts[1] in ['Manuscript','Supplement','canonical_render']:return False
 if f.name=='access_probe.py':return False
 return True
def main():
 files=sorted(f for f in P.rglob('*') if f.is_file() and include(f));records=[]
 for f in files:
  with f.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
  records.append({'path':f.relative_to(P).as_posix(),'bytes':f.stat().st_size,'sha256':h})
 manifest={'project':P.name,'date':'2026-11-17','file_count':len(files),'files':records,'omitted':'Earlier package archives, uncompressed XPT files, internal rendered QA page images, caches and temporary files. Original public ZIP source archives and full analysis extracts are included.'}
 (P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
 with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files+[P/'MANIFEST.json']:z.write(f,arcname=P.name+'/'+f.relative_to(P).as_posix())
 with zipfile.ZipFile(OUT) as z:
  assert z.testzip() is None
  for r in records:assert hashlib.sha256(z.read(P.name+'/'+r['path'])).hexdigest()==r['sha256'],r['path']
 with OUT.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
 (P/'PACKAGE_SHA256.txt').write_text(h+'  '+OUT.name+'\n',encoding='ascii');print(json.dumps({'archive':str(OUT),'bytes':OUT.stat().st_size,'files':len(files)+1,'sha256':h,'integrity':'PASS'},indent=2))
if __name__=='__main__':main()
