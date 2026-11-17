from pathlib import Path
import json
import pymupdf as fitz
from PIL import Image,ImageDraw
P=Path(__file__).resolve().parents[1]
def main():
 report={}
 for name in ['Manuscript','Supplement']:
  doc=fitz.open(P/'deliverables'/(name+'.pdf'));out=P/'qa'/name;out.mkdir(exist_ok=True)
  images=[];pages=[]
  for i,page in enumerate(doc):
   path=out/f'page-{i+1:02d}.png';page.get_pixmap(dpi=150).save(path);images.append(path)
   text=page.get_text();bad=[]
   for block in page.get_text('dict')['blocks']:
    for line in block.get('lines',[]):
     for span in line['spans']:
      x0,y0,x1,y1=span['bbox']
      if x0<-1 or y0<-1 or x1>page.rect.width+1 or y1>page.rect.height+1:bad.append(span['text'])
   pages.append({'page':i+1,'characters':len(text),'replacement_characters':text.count('\ufffd'),'out_of_page_text':bad})
  for i in range(0,len(images),2):
   ims=[Image.open(f).convert('RGB') for f in images[i:i+2]];w=sum(im.width for im in ims)+30;h=max(im.height for im in ims)+35;sheet=Image.new('RGB',(w,h),'#D8D8D8');dr=ImageDraw.Draw(sheet);x=10
   for k,im in enumerate(ims):sheet.paste(im,(x,25));dr.text((x+10,6),f'{name} | page {i+k+1}',fill='black');x+=im.width+10
   sheet.save(out/f'pair-{i//2+1:02d}.png')
  report[name]={'pages':len(doc),'page_details':pages}
 (P/'qa/page_layout_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:dict(pages=v['pages'],out_of_page=sum(len(z['out_of_page_text']) for z in v['page_details']),replacement_chars=sum(z['replacement_characters'] for z in v['page_details'])) for k,v in report.items()},indent=2))
if __name__=='__main__':main()
