#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import pandas as pd
from typeguard import typechecked

CODON_LENGTH = 3
STOP_CODON_NT = 3


@typechecked
def main(
    best_orf_path: str,
    sqanti_protein_path: str,
    output_prefix: str,
) -> int:
    orfs = pd.read_csv(
        # f"{output_prefix}_best_orf.tsv", sep="\t")
        f"{best_orf_path}",
        sep="\t",
    ).sort_values("transcript_id")
    sqanti_protein = pd.read_csv(
        # f"{output_prefix}.sqanti_protein_classification.tsv", sep="\t"
        f"{sqanti_protein_path}",
        sep="\t",
    ).sort_values("pb")
    pd.DataFrame(
        {
            "transcript_id": orfs["transcript_id"].values,
            "transcript_length_nt": orfs["len"].values,
            "orf_length_nt": orfs["orf_len"].values,
            "protein_length_cd": (orfs["orf_len"].values - STOP_CODON_NT)
            / CODON_LENGTH,
            "protein_is_nmd": sqanti_protein.is_nmd.values,
            "protein_splice_category": sqanti_protein.pr_splice_cat.values,
            "protein_splice_subcategory": sqanti_protein.pr_splice_subcat.values,
            "number_of_junctions_after_stop_codon": sqanti_protein.num_junc_after_stop_codon.values,
            "number_of_nt_after_stop_codon": sqanti_protein.num_nt_after_stop_codon.values,
        }
    ).to_csv(f"{output_prefix}_protein_annotation.tsv", index=False, sep="\t")
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--best_orf_path")
    parser.add_argument("--sqanti_protein_path")
    parser.add_argument("--output_prefix")

    args = parser.parse_args()
    main(
        best_orf_path=args.best_orf_path,
        sqanti_protein_path=args.sqanti_protein_path,
        output_prefix=args.output_prefix,
    )
