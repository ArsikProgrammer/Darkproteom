"""Build darkness annotation from downloaded pLDDT cache."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations", nargs="+", required=True)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    cache = pd.read_csv("data/processed/_plddt_needed.csv")
    mapping = pd.read_csv("data/processed/_uniprot_gene_accession.csv")

    # unique (gene, position) needed
    frames = [pd.read_csv(m, low_memory=False)[["Hugo_Symbol", "Protein_position"]] for m in args.mutations]
    mut = pd.concat(frames, ignore_index=True)
    mut["Protein_position"] = pd.to_numeric(mut["Protein_position"], errors="coerce")
    mut = mut.dropna(subset=["Protein_position"])
    mut["Protein_position"] = mut["Protein_position"].astype(int)
    uniq = mut.drop_duplicates()

    rows = []
    for gene, grp in uniq.groupby("Hugo_Symbol"):
        positions = set(grp["Protein_position"].tolist())
        accs = mapping[mapping["gene"] == gene]["accession"].tolist()
        for pos in positions:
            val = None
            for acc in accs:
                hit = cache[(cache["accession"] == acc) & (cache["position"] == pos)]
                if len(hit):
                    val = float(hit.iloc[0]["plddt"])
                    break
            if val is None:
                continue
            rows.append({
                "Hugo_Symbol": gene,
                "Protein_position": pos,
                "is_dark": int(val < args.threshold),
                "plddt": val,
                "source": "alphafold_plddt",
            })

    out = pd.DataFrame(rows).drop_duplicates()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print("Rows:", len(out), "dark:", int(out.is_dark.sum()) if len(out) else 0)
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
