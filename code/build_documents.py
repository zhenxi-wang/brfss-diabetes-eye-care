"""Create editable journal manuscript and supplement with the bundled runtime."""
from pathlib import Path
import json,re
from PIL import Image
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION_START,WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_ALIGN_VERTICAL
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE as RT
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'deliverables'

def make(title,main=False):
    d=Document();sec=d.sections[0]
    sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=Inches(.85);sec.left_margin=sec.right_margin=Inches(1)
    sec.footer_distance=Inches(.35)
    for name in ['Normal','Title','Subtitle','Heading 1','Heading 2','Heading 3','Caption']:
        s=d.styles[name];s.font.name='Times New Roman';s.font.color.rgb=RGBColor(0,0,0)
        rp=s.element.get_or_add_rPr();fonts=rp.find(qn('w:rFonts'))
        if fonts is not None:
            for key in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme']:fonts.attrib.pop(qn('w:'+key),None)
        pp=s.element.find(qn('w:pPr'))
        if pp is not None:
            for b in list(pp.findall(qn('w:pBdr'))):pp.remove(b)
    normal=d.styles['Normal'];normal.font.size=Pt(12 if main else 11.5)
    normal.paragraph_format.line_spacing=2 if main else 1.12
    normal.paragraph_format.space_after=Pt(0 if main else 6)
    normal.paragraph_format.widow_control=True
    d.styles['Title'].font.size=Pt(18);d.styles['Title'].font.bold=True;d.styles['Title'].paragraph_format.line_spacing=1.1;d.styles['Title'].paragraph_format.space_after=Pt(16)
    for name,size in [('Heading 1',13),('Heading 2',12),('Heading 3',11.5)]:
        s=d.styles[name];s.font.size=Pt(size);s.font.bold=True;s.paragraph_format.space_before=Pt(12);s.paragraph_format.space_after=Pt(6);s.paragraph_format.keep_with_next=True;s.paragraph_format.line_spacing=1.1
    d.styles['Caption'].font.size=Pt(11);d.styles['Caption'].font.italic=False;d.styles['Caption'].font.bold=False;d.styles['Caption'].paragraph_format.line_spacing=1.15
    foot=sec.footer.paragraphs[0];foot.alignment=WD_ALIGN_PARAGRAPH.CENTER
    field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');foot._p.append(field)
    d.core_properties.title=title;d.core_properties.author='';d.core_properties.subject='Paired survey ascertainment of diabetes eye care';d.core_properties.keywords='diabetes; retinal photography; dilation; BRFSS; surveillance; equity'
    d.add_paragraph(title,'Title');return d

def para(d,text,lead=None,style=None):
    text=re.sub(r'(?<![\w/])-(?=\d)','−',text)
    p=d.add_paragraph(style=style)
    if lead:p.add_run(lead+'. ').bold=True
    parts=re.split(r'(https://\S+|\[\d+(?:[-,–]\d+)*\])',text)
    for part in parts:
        if part.startswith('https://'):
            link=OxmlElement('w:hyperlink');link.set(qn('r:id'),p.part.relate_to(part,RT.HYPERLINK,is_external=True))
            rr=OxmlElement('w:r');pr=OxmlElement('w:rPr');cc=OxmlElement('w:color');cc.set(qn('w:val'),'000000');pr.append(cc);rr.append(pr);tt=OxmlElement('w:t');tt.text=part;rr.append(tt);link.append(rr);p._p.append(link)
        elif re.fullmatch(r'\[\d+(?:[-,–]\d+)*\]',part):
            run=p.add_run();run.font.superscript=True
            for character in part[1:-1]:
                if character in '-–':run._r.append(OxmlElement('w:noBreakHyphen'))
                else:run.add_text(character)
        else:p.add_run(part)
    return p

def table(d,t,newpage=True):
    h=d.add_heading(t['title'],2);h.paragraph_format.page_break_before=newpage
    tb=d.add_table(rows=1,cols=len(t['headers']));tb.alignment=WD_TABLE_ALIGNMENT.CENTER;tb.autofit=False
    for c,width in zip(tb.columns,t['widths']):c.width=Inches(width)
    for c,s in zip(tb.rows[0].cells,t['headers']):c.text=s
    for row in t['rows']:
        cells=tb.add_row().cells
        for c,s in zip(cells,row):c.text=str(s)
    for i,row in enumerate(tb.rows):
        trpr=row._tr.get_or_add_trPr();trpr.append(OxmlElement('w:cantSplit'))
        if i==0:trpr.append(OxmlElement('w:tblHeader'))
        for j,c in enumerate(row.cells):
            c.width=Inches(t['widths'][j]);c.vertical_alignment=WD_ALIGN_VERTICAL.CENTER
            cp=c._tc.get_or_add_tcPr();sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'E7EDF3' if i==0 else 'FFFFFF');cp.append(sh)
            mar=OxmlElement('w:tcMar')
            is_definitions=t['title'].startswith('eTable 1 Operational')
            pad='35'
            for side,val in [('top',pad),('bottom',pad),('left','85'),('right','85')]:
                mm=OxmlElement('w:'+side);mm.set(qn('w:w'),val);mm.set(qn('w:type'),'dxa');mar.append(mm)
            cp.append(mar);bo=OxmlElement('w:tcBorders')
            for side in ['top','bottom','left','right']:
                bb=OxmlElement('w:'+side);bb.set(qn('w:val'),'single');bb.set(qn('w:sz'),'3');bb.set(qn('w:color'),'D9D9D9');bo.append(bb)
            cp.append(bo)
            for p in c.paragraphs:
                p.paragraph_format.space_after=Pt(0);p.paragraph_format.space_before=Pt(0);p.paragraph_format.line_spacing=1.04 if is_definitions else 1.07
                prose_table=t['title'].startswith(('eTable 1 Operational','eTable 10 STROBE','Supplementary Table S9','Supplementary Table S10'))
                p.alignment=WD_ALIGN_PARAGRAPH.LEFT if j==0 or prose_table else WD_ALIGN_PARAGRAPH.CENTER
                for run in p.runs:run.font.size=Pt(10 if is_definitions else 10.5);run.bold=i==0
                if t['title'].startswith('eTable 2 Complete') and i==23 and j==0:
                    p.paragraph_format.page_break_before=True
    p=para(d,t['note']);p.paragraph_format.line_spacing=1.12;p.paragraph_format.space_before=Pt(8);p.paragraph_format.space_after=Pt(8)
    for run in p.runs:run.font.size=Pt(10.5)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    x=json.loads((ROOT/'manuscript/manuscript.json').read_text(encoding='utf8'));d=make(x['title'],True)
    p=para(d,x['subtitle']);p.paragraph_format.line_spacing=1.15
    p=para(d,f"Main text {x['main_word_count']} words | Abstract {x['abstract_word_count']} words");p.paragraph_format.line_spacing=1.15
    h=d.add_heading('Abstract',1);h.paragraph_format.page_break_before=True
    for h,t in x['abstract']:para(d,t,lead=h)
    for i,(h,ps) in enumerate(x['sections']):
        hh=d.add_heading(h,1);hh.paragraph_format.page_break_before=i==0
        for text in ps:para(d,text)
    for title,text in x['declarations']:
        d.add_heading(title,1);para(d,text)
    hh=d.add_heading('References',1);hh.paragraph_format.page_break_before=True
    for i,text in enumerate(x['references'],1):
        p=para(d,f'{i}. '+text);p.paragraph_format.keep_together=True
    for t in x['tables']:table(d,t)
    hh=d.add_heading('Figure Legends',1);hh.paragraph_format.page_break_before=True
    for f in x['figures']:
        d.add_heading(f['title'],2);p=para(d,f['caption'])
        for run in p.runs:run.font.size=Pt(11)
    # Journal submission file keeps figures separate.
    d.save(OUT/'Manuscript.docx')
    y=json.loads((ROOT/'manuscript/supplement.json').read_text(encoding='utf8'));d=make(y['title'])
    para(d,y['references_note'])
    for h,ps in y['sections']:
        d.add_heading(h,1)
        for text in ps:para(d,text)
    for idx,t in enumerate(y['tables']):table(d,t,newpage=not t['title'].startswith(('Supplementary Table S2.','Supplementary Table S10.')))
    for i,f in enumerate(y['figures']):
        sec=d.add_section(WD_SECTION_START.NEW_PAGE)
        portrait=f['file'].startswith('FigureS2')
        sec.orientation=WD_ORIENT.PORTRAIT if portrait else WD_ORIENT.LANDSCAPE
        sec.page_width=Inches(8.5 if portrait else 11);sec.page_height=Inches(11 if portrait else 8.5)
        sec.left_margin=sec.right_margin=Inches(1 if portrait else .7);sec.top_margin=sec.bottom_margin=Inches(.65)
        hh=d.add_heading(f['title'],2)
        p=d.add_paragraph();p.paragraph_format.keep_with_next=True
        path=ROOT/'figures'/f['file'];im=Image.open(path);width=min(6.5 if portrait else 9.5,(7.4 if portrait else 5.55)*im.width/im.height)
        p.alignment=WD_ALIGN_PARAGRAPH.CENTER
        run=p.add_run();run.add_picture(str(path),width=Inches(width))
        d.inline_shapes[-1]._inline.docPr.set('descr',f['caption'])
        para(d,f['caption'],style='Caption')
    d.save(OUT/'Supplement.docx')
    print('Created two editable documents')

if __name__=='__main__':main()
