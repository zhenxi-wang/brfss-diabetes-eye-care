"""Export paginated PDF reading copies from the same source as the editable documents."""
from argparse import ArgumentParser
from html import escape
from io import BytesIO
from pathlib import Path
import json
import re
import pymupdf as pdf

ROOT = Path(__file__).resolve().parents[1]
STYLE = """
body { font-family: serif; margin: 0; color: #15191d; font-size: 11pt; line-height: 1.25; }
h1 { font-size: 17pt; line-height: 1.15; margin: 0 0 15pt; }
h2 { font-size: 13pt; line-height: 1.15; margin: 0 0 9pt; }
p { margin: 0 0 8pt; }
table { width: 100%; border-collapse: collapse; font-size: 9.1pt; line-height: 1.15; }
th, td { border: 1px solid #cbd2d8; padding: 4pt 4pt; vertical-align: middle; }
th { background: #e7edf3; font-weight: bold; }
td { background: white; }
.note { margin-top: 8pt; font-size: 9.5pt; }
"""


def document_story(writer, markup, pagesize=(612, 792)):
    story = pdf.Story(markup, user_css=STYLE)
    page = pdf.Rect(0, 0, *pagesize)
    body = page + (54, 49, -54, -47)
    more = True
    while more:
        device = writer.begin_page(page)
        more, filled = story.place(body)
        if pdf.Rect(filled).height == 0:
            raise RuntimeError('Supplement contains content that does not fit a page')
        story.draw(device)
        writer.end_page()


def table_markup(table):
    heading = '<h2>' + escape(table['title']) + '</h2>'
    header = ''.join(f'<th width="{w/6.5:.2%}">{escape(h)}</th>'
                     for w, h in zip(table['widths'], table['headers']))
    rows = ['<tr>' + ''.join(f'<td>{escape(v)}</td>' for v in row) + '</tr>' for row in table['rows']]
    # Balance long tables across pages and repeat their column headings.
    limit = 24
    if table['title'].startswith('Supplementary Table S2'):
        limit = 18
    if table['title'].startswith('Supplementary Table S11'):
        limit = 12
    if table['title'].startswith('Supplementary Table S20'):
        limit = 12
    segments = (len(rows) + limit - 1) // limit
    chunk = (len(rows) + segments - 1) // segments
    for start in range(0, len(rows), chunk):
        title = heading if start == 0 else heading.replace('</h2>', ' continued</h2>')
        note = f'<p class="note">{escape(table["note"])}</p>' if start+chunk >= len(rows) else ''
        yield title + '<table><tr>' + header + '</tr>' + ''.join(rows[start:start+chunk]) + '</table>' + note


def save_reading_copy(document, path, title, subject):
    for i, page in enumerate(document, 1):
        page.insert_textbox(pdf.Rect(0, page.rect.height-28, page.rect.width, page.rect.height-7),
                            str(i), fontname='tiro', fontsize=10, align=1)
    document.set_metadata({'title': title, 'subject': subject})
    document.save(path, garbage=4, deflate=True)


def export_supplement():
    source = json.loads((ROOT / 'manuscript/supplement.json').read_text(encoding='utf8'))
    path = ROOT / 'deliverables/Supplement.pdf'
    buffer = BytesIO()
    writer = pdf.DocumentWriter(buffer)
    text = f'<h1>{escape(source["title"])}</h1><p>{escape(source["references_note"])}</p>'
    for title, paragraphs in source['sections']:
        text += '<h2>' + escape(title) + '</h2>'
        text += ''.join(f'<p>{escape(paragraph)}</p>' for paragraph in paragraphs)
    document_story(writer, text)
    for table in source['tables']:
        for fragment in table_markup(table):
            document_story(writer, fragment)
    writer.close()
    with pdf.open(stream=buffer.getvalue(), filetype='pdf') as document:
        for figure in source['figures']:
            landscape = not figure['file'].startswith('FigureS2')
            page = document.new_page(width=792 if landscape else 612, height=612 if landscape else 792)
            width, height = page.rect.width, page.rect.height
            page.insert_text((54, 59), figure['title'], fontname='tibo', fontsize=13)
            bounds = pdf.Rect(50, 77, width-50, height-110)
            image = ROOT / 'figures' / Path(figure['file']).with_suffix('.pdf')
            with pdf.open(image) as vector:
                page.show_pdf_page(bounds, vector, 0)
            page.insert_htmlbox(pdf.Rect(54, height-96, width-54, height-36),
                                f'<p>{escape(figure["caption"])}</p>', css=STYLE)
        save_reading_copy(document, path, source['title'], 'Supplementary analyses and figures')


def prose(value):
    text = re.sub(r'(?<![\w/])-(?=\d)', '−', value)
    text = re.sub(r'\[(\d+(?:[-–,]\d+)*)\]', r'<sup>\1</sup>', escape(text))
    return re.sub(r'(https://[^\s<]+)', r'<a href="\1">\1</a>', text)


def export_manuscript():
    source = json.loads((ROOT / 'manuscript/manuscript.json').read_text(encoding='utf8'))
    buffer = BytesIO()
    writer = pdf.DocumentWriter(buffer)
    cover = f'<h1>{escape(source["title"])}</h1><p>{escape(source["subtitle"])}</p>'
    cover += f'<p>Main text {source["main_word_count"]:,} words | Abstract {source["abstract_word_count"]:,} words</p>'
    abstract = '<h2>Abstract</h2>' + ''.join(f'<p><b>{escape(lead)}.</b> {prose(text)}</p>'
                                              for lead, text in source['abstract'])
    document_story(writer, cover + abstract)
    main = ''
    for title, paragraphs in [*source['sections'], *[(h, [p]) for h, p in source['declarations']]]:
        main += '<h2>' + escape(title) + '</h2>'
        main += ''.join(f'<p>{prose(paragraph)}</p>' for paragraph in paragraphs)
    references = '<h2>References</h2>' + ''.join(f'<p>{i}. {prose(ref)}</p>'
                                                 for i, ref in enumerate(source['references'], 1))
    document_story(writer, main + references)
    for table in source['tables']:
        for fragment in table_markup(table):
            document_story(writer, fragment)
    legends = '<h2>Figure legends</h2>' + ''.join(
        f'<h2>{escape(item["title"])}</h2><p>{prose(item["caption"])}</p>' for item in source['figures'])
    document_story(writer, legends)
    writer.close()
    with pdf.open(stream=buffer.getvalue(), filetype='pdf') as document:
        save_reading_copy(document, ROOT / 'deliverables/Manuscript.pdf', source['title'],
                          '2025 BRFSS analysis with separate 2022 comparison')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--manuscript', action='store_true', help='Also regenerate the manuscript reading copy.')
    args = parser.parse_args()
    export_supplement()
    if args.manuscript:
        export_manuscript()

