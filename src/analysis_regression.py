"""Regress dark-mutation fraction on smoking signature and covariates."""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="sample features CSV")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, low_memory=False)
    df = df.dropna(subset=["dark_fraction", "SBS4", "total_mutations"])

    # dataset dummies (george as reference)
    if "dataset" in df.columns:
        dummies = pd.get_dummies(df["dataset"], prefix="dataset", drop_first=True, dtype=float)
        df = pd.concat([df, dummies], axis=1)

    predictors = ["SBS4", "total_mutations"] + [c for c in df.columns if c.startswith("dataset_")]
    X = df[predictors].astype(float)
    X = sm.add_constant(X)
    y = df["dark_fraction"].astype(float)

    model = sm.OLS(y, X).fit()
    result = pd.DataFrame({
        "term": model.params.index,
        "coef": model.params.values,
        "pvalue": model.pvalues.values,
        "ci_low": model.conf_int()[0].values,
        "ci_high": model.conf_int()[1].values,
    })
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(model.summary().as_text()[:1500])
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
