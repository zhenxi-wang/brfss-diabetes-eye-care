"""Extract a fixed panel from the unmodified official XPT archive."""
from pathlib import Path
import zipfile,shutil,json,hashlib
import pandas as pd
P=Path(__file__).resolve().parents[1]
COLS=['_STATE','SEQNO','IYEAR','_PSU','_STSTR','_LLCPWT','DIABETE4','DIABTYPE','INSULIN1','INSULIN','CHKHEMO3','EYEEXAM1','DIABEYE1','DIABAGE4','MEDCOST1','CHECKUP1','PERSDOC3','PRIMINSR','_HLTHPLN','_AGE80','_AGEG5YR','_AGE65YR','SEXVAR','_SEX','_RACE','_RACEGR3','_RACEGR4','_IMPRACE','_RACEPRV','_HISPANC','EDUCA','INCOME3','BLIND','_METSTAT','_URBSTAT','_BMI5','_SMOKER3','GENHLTH','DROCDY4_','DIABEDU1']
def main():
 for y in [2025,2022]:
  path=P/f'data/raw/LLCP{y}.xpt';archive=P/f'data/raw/LLCP{y}XPT.zip'
  if not path.exists():
   with zipfile.ZipFile(archive) as z:
    n=next(n for n in z.namelist() if n.strip().lower().endswith('.xpt'))
    with z.open(n) as a,path.open('wb') as b:shutil.copyfileobj(a,b,2**22)
  chunks=[];cols=None;total=0
  for d in pd.read_sas(path,format='xport',chunksize=50000,encoding='utf8'):
   if cols is None:
    cols=[c for c in [*COLS,'DISPCODE','QSTVER','FEETSORE'] if c in d.columns];(P/f'sources/{y}_all_variables.json').write_text(json.dumps(d.columns.tolist(),indent=2));print(y,'selected',cols,flush=True)
   chunks.append(d[cols]);total+=len(d)
  d=pd.concat(chunks,ignore_index=True);d.to_parquet(P/f'data/derived/brfss{y}_panel.parquet',index=False)
  with path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
  (P/f'sources/{y}_extraction.json').write_text(json.dumps({'records':total,'selected_columns':cols,'uncompressed_bytes':path.stat().st_size,'xpt_sha256':sha},indent=2));print(y,total,'rows',flush=True)
if __name__=='__main__':main()
