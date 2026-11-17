from pathlib import Path
from lxml import html
import re
P=Path(__file__).resolve().parents[1]
out=[]
for y in [2025,2022]:
 s=' '.join(html.fromstring((P/f'sources/{y}_codebook.html').read_bytes()).itertext());s=re.sub(r'\s+',' ',s)
 for c in ['DIABEYE1','EYEEXAM1','_RACEGR3','_RACEGR4','_IMPRACE','INCOME3','_METSTAT']:
  m=re.search(r'SAS Variable Name:\s*'+re.escape(c)+r'\b',s)
  part=s[m.start():m.start()+1800].split('Label:')[0] if m else 'NOT FOUND'
  out.append(f'{y}: {c}\n{part}')
(P/'sources/coding_excerpts.txt').write_text('\n\n'.join(out),encoding='utf8');print('\n\n'.join(out))
