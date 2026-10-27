"""Run the implemented part of the Dark Proteome pipeline."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
PYTHON = sys.executable

STEPS = [
    "01_inspect_george.py",
    "02_convert_jhu.py",
    "04_filter_clcgp.py",
    "05_filter_jhu.py",
    "03_decompress_hg19.py",
    "03_extract_primary.py",
    "09_parse_gencode.py",
]


def main():
    for step in STEPS:
        path = SRC / step
        print(f"\n=== Running {step} ===")
        result = subprocess.run([PYTHON, str(path)], cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"Step failed: {step}")

    # Context extraction is run per dataset via the unified CLI.
    context_jobs = [
        ("data/processed/george_coding_snps.csv", "data/processed/george_coding_snps_context.csv"),
        ("data/processed/sclc_clcgp_coding_snps.csv", "data/processed/sclc_clcgp_coding_snps_context.csv"),
        ("data/processed/sclc_jhu_coding_snps.csv", "data/processed/sclc_jhu_coding_snps_context.csv"),
    ]
    for inp, out in context_jobs:
        print(f"\n=== Context: {inp} -> {out} ===")
        result = subprocess.run(
            [PYTHON, str(SRC / "06_extract_context.py"), "--input", inp, "--output", out],
            cwd=ROOT,
        )
        if result.returncode != 0:
            raise SystemExit(f"Context step failed: {inp}")

    for step in ["07_build_context_counts.py", "08_estimate_signatures.py"]:
        path = SRC / step
        print(f"\n=== Running {step} ===")
        result = subprocess.run([PYTHON, str(path)], cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"Step failed: {step}")

    print("\nPipeline finished. Note: dark-proteome core analysis is not implemented yet.")


if __name__ == "__main__":
    main()
