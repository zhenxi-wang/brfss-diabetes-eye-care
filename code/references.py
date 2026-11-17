"""Verify bibliographic records and preserve a bounded, reproducible novelty search."""
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import quote
from concurrent.futures import ThreadPoolExecutor
import json
P=Path(__file__).resolve().parents[1]
REFS={
'ada2026':'DOI:10.2337/dc26-S012',
'telemedicine':'DOI:10.1001/jamaophthalmol.2015.1',
'implementation':'DOI:10.1001/jamainternmed.2017.0204',
'pivotalAI':'DOI:10.1038/s41746-018-0040-6',
'accessRCT':'DOI:10.1038/s41467-023-44676-z',
'equity2024':'DOI:10.1038/s41746-024-01197-3',
'followup2026':'DOI:10.1038/s41746-026-02460-5',
'Wagner2025':'DOI:10.1136/bmjph-2024-002202',
'rural2025':'EXT_ID:39693586 AND SRC:MED',
'validation':'TITLE:"Medical Record Validation of Self-Reported Eye Diseases and Eye Care Utilization among Older Adults"',
'agreement':'TITLE:"Agreement between self-reports and medical records was only fair in a cross-sectional study of performance of annual eye examinations among adults with diabetes in managed care"',
'STROBE':'DOI:10.1371/journal.pmed.0040296',
'burden':'TITLE:"Prevalence of Diabetic Retinopathy in the US in 2021"',
}
SEARCHES={
'paired_modality':'(BRFSS OR "Behavioral Risk Factor Surveillance System") AND ("retinal photography" OR "retinal imaging") AND (dilat* OR ascertain* OR misclassif* OR "either")',
'photo_coverage':'(BRFSS OR "Behavioral Risk Factor Surveillance System") AND ("retinal photography" OR "eye photo")',
'dilation_metric':'diabet* AND ("dilation only" OR "dilation-only" OR "dilated examination") AND (surveillance OR undercount* OR misclassif*)',
'cgm_vision':'(BRFSS OR NHIS) AND ("continuous glucose monitoring" OR CGM) AND (vision OR visual OR blind*)',
}
def query(q,size=20):
 u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=core&pageSize='+str(size)+'&query='+quote(q)
 with urlopen(u,timeout=90) as r:x=json.load(r)
 rows=x.get('resultList',{}).get('result',[])
 return {'query':q,'url':u,'date':'2026-11-12','hitCount':x.get('hitCount'),'results':[{k:v.get(k) for k in ['title','doi','pubYear','authorList','journalInfo','pageInfo','id','source']} for v in rows]}
def main():
 out={}
 with ThreadPoolExecutor(max_workers=4) as pool:
  keys=list(REFS)
  for key,v in zip(keys,pool.map(lambda k:query(REFS[k],3),keys)):
   if not v['results']:print('NO MATCH',key,flush=True)
   else:out[key]=v['results'][0];print(key,out[key]['title'],out[key]['doi'],flush=True)
 (P/'sources/reference_metadata.json').write_text(json.dumps(out,indent=2),encoding='utf8')
 with ThreadPoolExecutor(max_workers=4) as pool:audits=dict(zip(SEARCHES,pool.map(lambda k:query(SEARCHES[k],100),SEARCHES)))
 (P/'sources/novelty_search.json').write_text(json.dumps(audits,indent=2),encoding='utf8');print({k:v['hitCount'] for k,v in audits.items()},flush=True)
if __name__=='__main__':main()
