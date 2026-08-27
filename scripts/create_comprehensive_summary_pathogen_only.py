import pandas as pd

pathogen_xlsx = snakemake.input["pathogen"]
out_xlsx = snakemake.output["excel"]

# In screening-only mode we may not have host/mtDNA summaries, but we still want a
# single "comprehensive" workbook with at least the pathogen sheet.
try:
    xl = pd.ExcelFile(pathogen_xlsx)
    sheet_names = xl.sheet_names
except Exception:
    sheet_names = []

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    if sheet_names:
        for sheet in sheet_names:
            df = pd.read_excel(pathogen_xlsx, sheet_name=sheet)
            df.to_excel(writer, index=False, sheet_name=sheet[:31])
    else:
        df = pd.read_excel(pathogen_xlsx)
        df.to_excel(writer, index=False, sheet_name="pathogen_summary")

