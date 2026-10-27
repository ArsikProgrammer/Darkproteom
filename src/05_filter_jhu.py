"""Filter JHU mutations to coding SNPs."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "processed" / "sclc_jhu_mutations.csv"
OUT = ROOT / "data" / "processed" / "sclc_jhu_coding_snps.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(RAW, low_memory=False)
coding = {
    "Missense_Mutation", "Silent", "Nonsense_Mutation",
    "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "Nonstop_Mutation", "NON_SYNONYMOUS_CODING",
}
out = df[(df.Variant_Type == "SNP") & (df.Variant_Classification.isin(coding))].copy()
keep = [
    "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
    "End_Position", "Variant_Classification", "Variant_Type",
    "Reference_Allele", "Tumor_Seq_Allele2", "Transcript_ID", "HGVSp_Short",
    "Protein_position",
]
out[keep].to_csv(OUT, index=False)
print("JHU coding SNPs:", len(out), "samples:", out.Tumor_Sample_Barcode.nunique())
print("Saved:", OUT)
