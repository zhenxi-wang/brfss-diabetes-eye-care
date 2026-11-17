"""Join the revised manuscript and the three vector figure PDFs."""
from pathlib import Path
import pymupdf as pdf

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ['Figure1_paired_eye_care_classification', 'Figure2_photography_only_by_group',
           'Figure3_adjusted_contact_disparities']

with pdf.open(ROOT / 'deliverables/Manuscript.pdf') as document:
    toc = [[1, 'Manuscript tables and figure legends', 1]]
    for number, name in enumerate(FIGURES, 1):
        toc.append([1, f'Figure {number}', len(document) + 1])
        with pdf.open(ROOT / 'figures' / (name + '.pdf')) as figure:
            document.insert_pdf(figure)
    document.set_toc(toc)
    document.set_metadata({
        'title': 'Retinal photography and the ascertainment of diabetes eye care',
        'subject': '2025 BRFSS paired eye-care classification with separate 2022 analysis',
    })
    document.save(ROOT / 'deliverables/Illustrated_review.pdf', garbage=4, deflate=True)
    print(f'Illustrated reading copy: {len(document)} pages')

