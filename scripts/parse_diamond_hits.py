#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def score_hit(row):
    score = 0
    ev = float(row["evalue"])
    pident = float(row["pident"])
    aln_len = int(row["length"])
    qspan = int(row["qspan_nt"])

    if ev < 1e-10:
        score += 5
    elif ev < 1e-7:
        score += 3
    elif ev <= 1e-5:
        score += 1

    if pident >= 90:
        score += 3
    elif pident >= 80:
        score += 2
    elif pident >= 70:
        score += 1

    if 50 <= qspan <= 150:
        score += 3
    elif 30 <= qspan <= 200:
        score += 1

    if aln_len >= 60:
        score += 2
    elif aln_len >= 30:
        score += 1

    return score


def classify(score):
    if score >= 12:
        return "high_confidence"
    if score >= 8:
        return "medium_confidence"
    if score >= 5:
        return "low_confidence"
    return "rejected"


def main():
    p = argparse.ArgumentParser(description="Filter and score DIAMOND hits for aDNA screening.")
    p.add_argument("--sample", required=True)
    p.add_argument("--in-tsv", required=True)
    p.add_argument("--out-all", required=True)
    p.add_argument("--out-filtered", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--min-aa-len", type=int, default=30)
    p.add_argument("--min-pident", type=float, default=70.0)
    p.add_argument("--min-qspan-nt", type=int, default=90)
    args = p.parse_args()

    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "staxids",
        "stitle",
    ]

    df = pd.read_csv(args.in_tsv, sep="\t", header=None, names=cols)
    if df.empty:
        os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)
        pd.DataFrame(
            [{
                "sample": args.sample,
                "total_hits": 0,
                "filtered_hits": 0,
                "high_confidence": 0,
                "medium_confidence": 0,
                "low_confidence": 0,
            }]
        ).to_csv(args.out_summary, index=False)
        df.to_csv(args.out_all, sep="\t", index=False)
        df.to_csv(args.out_filtered, sep="\t", index=False)
        return

    df["qspan_nt"] = (df["qend"] - df["qstart"]).abs() + 1
    df["sample"] = args.sample
    df["authenticity_score"] = df.apply(score_hit, axis=1)
    df["confidence"] = df["authenticity_score"].map(classify)

    filtered = df[
        (df["evalue"] <= args.evalue)
        & (df["length"] >= args.min_aa_len)
        & (df["pident"] >= args.min_pident)
        & (df["qspan_nt"] >= args.min_qspan_nt)
    ].copy()

    filtered["authenticity_score"] = filtered.apply(score_hit, axis=1)
    filtered["confidence"] = filtered["authenticity_score"].map(classify)

    summary = {
        "sample": args.sample,
        "total_hits": int(len(df)),
        "filtered_hits": int(len(filtered)),
        "high_confidence": int((filtered["confidence"] == "high_confidence").sum()),
        "medium_confidence": int((filtered["confidence"] == "medium_confidence").sum()),
        "low_confidence": int((filtered["confidence"] == "low_confidence").sum()),
    }

    os.makedirs(os.path.dirname(args.out_all), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_filtered), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    df.to_csv(args.out_all, sep="\t", index=False)
    filtered.to_csv(args.out_filtered, sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(args.out_summary, index=False)


if __name__ == "__main__":
    main()
