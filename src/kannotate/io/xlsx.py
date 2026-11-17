from pathlib import Path

from openpyxl import Workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from .. import __version__
from ..cancellation import check_cancelled
from ..scoring import SCORE_VERSION
from .tabular import GENBANK_HEADERS, HEADERS, annotation_rows


def write_xlsx(annotations, path: Path, *, database_versions=None, genbank=None, input_name=None):
    annotations = list(annotations)
    headers = GENBANK_HEADERS if genbank is not None else HEADERS
    workbook = Workbook()
    workbook.properties.creator = "Kannotate"
    workbook.properties.description = (
        f"Kannotate {__version__}; score {SCORE_VERSION}; {database_versions or {}}"
    )
    sheet = workbook.active
    sheet.title = "Annotations"
    sheet.append(headers)
    if genbank is not None:
        sheet["D1"].comment = Comment(
            "From/To: 1-based inclusive lower/upper reported bounds. Joined CDSs use "
            "outer bounds; partial/complex locations have an original-location cell note. "
            "Unknown or remote bounds are blank. Strand is +, -, or blank when unknown/mixed.",
            "Kannotate",
        )
    for row_index, (values, location_note) in enumerate(
        annotation_rows(annotations, genbank=genbank, input_name=input_name), 2
    ):
        sheet.append(values)
        if location_note:
            sheet.cell(row_index, 4).comment = Comment(location_note, "Kannotate")
        # Prevent arbitrary headers/descriptions from becoming spreadsheet formulas.
        for cell in sheet.iter_rows(min_row=row_index, max_row=row_index, max_col=len(headers)):
            for item in cell:
                if isinstance(item.value, str):
                    item.data_type = "s"
        # Each row index is known; avoid scanning the entire growing worksheet.
        sheet.row_dimensions[row_index].height = 45
    for cell in sheet[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="34495E")
        cell.alignment = Alignment(wrap_text=True, vertical="center")
    sheet.row_dimensions[1].height = 34
    widths = (
        [24, 24, 24, 12, 12, 10, 52, 44, 52, 18] if genbank is not None else [38, 44, 52, 18]
    ) + [28, 58, 16, 18, 22] * 4
    for column, width in enumerate(widths, 1):
        sheet.column_dimensions[get_column_letter(column)].width = width
    body_font = Font(name="Arial", size=11)
    plain = Alignment(vertical="top")
    wrapped = Alignment(vertical="top", wrap_text=True)
    for row in sheet.iter_rows(min_row=2):
        check_cancelled()
        for cell in row:
            cell.font = body_font
            cell.alignment = plain
        for header, cell in zip(headers, row):
            if header != "Normalized/search sequence":
                cell.alignment = wrapped
            if header.endswith("score"):
                cell.number_format = "0.00"
            elif header.endswith("E-value"):
                cell.number_format = "0.00E+00"
    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions
    check_cancelled()
    workbook.save(path)
    check_cancelled()
