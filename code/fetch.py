"""Retrieve official public BRFSS archives and documentation with hashes."""
from pathlib import Path
from urllib.request import urlopen,Request
from urllib.parse import urljoin
from concurrent.futures import ThreadPoolExecutor
from lxml import html
import json,hashlib,zipfile,io,time
from urllib.error import HTTPError
import pandas as pd
import pymupdf
P=Path(__file__).resolve().parents[1]
PINNED={'LLCP2025XPT.zip':'d26a8b391aa2966b58dc2a3ffddf8da050e9f38eb143d0f5ca7b3d21ece08f98','LLCP2022XPT.zip':'b6319c949b7289f8020fc903c3c2d38fe7e0b28ad703dcd81547f5cd0c0dcfeb'}

def get(url,path):
 path=P/path;path.parent.mkdir(parents=True,exist_ok=True)
 retrieval_url=url.replace('www.cdc.gov','portugues.cdc.gov') if path.exists() else url
 if not path.exists():
  for attempt in range(4):
   try:
    try:r=urlopen(Request(url,headers={'User-Agent':'Mozilla/5.0','Accept':'*/*'}),timeout=180)
    except HTTPError as e:
     if e.code!=403 or 'www.cdc.gov' not in url:raise
     retrieval_url=url.replace('www.cdc.gov','portugues.cdc.gov');r=urlopen(Request(retrieval_url,headers={'User-Agent':'Mozilla/5.0'}),timeout=180)
    with r, path.with_suffix(path.suffix+'.tmp').open('wb') as out:
     while b:=r.read(2**20):out.write(b)
    path.with_suffix(path.suffix+'.tmp').replace(path);break
   except Exception:
    if attempt==3:raise
    time.sleep(2)
 with path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
 if path.name in PINNED and sha!=PINNED[path.name]:raise ValueError('Source archive differs from the pinned research version: '+path.name)
 return {'url':url,'retrieval_url':retrieval_url,'path':path.relative_to(P).as_posix(),'bytes':path.stat().st_size,'sha256':sha,'retrieved':'2026-11-12'}

def main():
 records=[];jobs=[]
 for year in [2025,2022]:
  url=f'https://www.cdc.gov/brfss/annual_data/annual_{year}.html';r=get(url,f'sources/annual_{year}.html');records.append(r);root=html.fromstring((P/r['path']).read_bytes())
  jobs.append((f'https://www.cdc.gov/brfss/annual_data/{year}/files/LLCP{year}XPT.zip',f'data/raw/LLCP{year}XPT.zip'))
  for a in root.xpath('//a[@href]'):
   label=' '.join(a.itertext()).strip();href=urljoin(url,a.get('href'))
   if 'Codebook' in label:jobs.append((href,f'sources/{year}_codebook'+('.zip' if '.zip' in href.lower() else '.pdf')))
   if 'Complex Sampling Weights' in label:jobs.append((href,f'sources/{year}_module_weights.pdf'))
   if 'Overview' in label:jobs.append((href,f'sources/{year}_overview.pdf'))
  jobs.extend([(f'https://www.cdc.gov/brfss/questionnaires/modules/state{year}.htm',f'sources/module_states_{year}.html'),(f'https://www.cdc.gov/brfss/questionnaires/pdf-ques/{year}-BRFSS-Questionnaire-508.pdf',f'sources/questionnaire_{year}.pdf')])
 with ThreadPoolExecutor(max_workers=5) as pool:
  for r in pool.map(lambda z:get(*z),jobs):records.append(r);print(r['path'],r['bytes'],flush=True)
 for year in [2025,2022]:
  f=P/f'sources/{year}_codebook.zip'
  if f.exists():
   with zipfile.ZipFile(f) as z:
    for n in z.namelist():
     if n.lower().endswith('.pdf'):(P/f'sources/{year}_codebook.pdf').write_bytes(z.read(n))
     if n.lower().endswith(('.html','.htm')):
      b=z.read(n);(P/f'sources/{year}_codebook.html').write_bytes(b);(P/f'sources/{year}_codebook.txt').write_text(html.fromstring(b).text_content(),encoding='utf8')
  for kind in ['codebook','module_weights','overview']:
   f=P/f'sources/{year}_{kind}.pdf'
   if f.exists():
    with pymupdf.open(f) as d:(P/f'sources/{year}_{kind}.txt').write_text('\n'.join(p.get_text() for p in d),encoding='utf8')
  with pymupdf.open(P/f'sources/questionnaire_{year}.pdf') as d:(P/f'sources/questionnaire_{year}.txt').write_text('\n'.join(p.get_text() for p in d),encoding='utf8')
  tabs=pd.read_html(io.StringIO((P/f'sources/module_states_{year}.html').read_text(encoding='utf8')));table=next(t for t in tabs if len(t.columns)>=5);table.columns=['state','description','dataset','weight','modules'];table['state']=table.state.ffill();table.to_csv(P/f'sources/module_states_{year}.csv',index=False)
  elig=table[table.modules.fillna('').str.split(',').map(lambda x:'Diabetes' in [z.strip() for z in x]) & (table.dataset==f'LLCP{year}') & ~table.state.isin(['Guam','Puerto Rico','Virgin Islands','Northern Mariana Isl'])]
  elig.to_csv(P/f'sources/eligible_states_{year}.csv',index=False);print(year,'full-module states:',elig.state.to_list(),flush=True)
 (P/'sources/download_manifest.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
