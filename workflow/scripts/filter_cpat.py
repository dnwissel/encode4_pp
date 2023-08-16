#!/usr/bin/env python
# -*- coding: utf-8 -*-


from argparse import ArgumentParser

import pandas as pd
from typeguard import typechecked


@typechecked
def main(
    input_file_path: str,
    output_path: str,
    first_cutoff: float,
    second_cutoff: float,
) -> int:
    df = pd.read_csv(f"{input_file_path}", sep="\t")
    df["transcript_id"] = df.ID.str.rsplit("_").str[0]
    df_first_cutoff = df.loc[df["Coding_prob"] > first_cutoff]
    df_second_cutoff = df.loc[df["Coding_prob"] > second_cutoff]
    pd.concat(
        [
            df_first_cutoff.sort_values(
                ["ORF", "Coding_prob"], ascending=False
            )
            .groupby("transcript_id")
            .head(1),
            df_second_cutoff.sort_values(
                ["ORF", "Coding_prob"], ascending=False
            )
            .groupby("transcript_id")
            .head(1),
            df.sort_values(["ORF", "Coding_prob"], ascending=False)
            .groupby("transcript_id")
            .head(1),
        ],
        axis=0,
        ignore_index=True,
    ).sort_values(["Coding_prob"], ascending=False).sort_values(
        "Coding_prob", ascending=False
    ).groupby(
        "transcript_id"
    ).head(
        1
    ).to_csv(
        f"{output_path}",
        sep="\t",
        index=False,
        header=True,
    )
    return 0


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--input_file_path")
    parser.add_argument("--output_path")
    parser.add_argument("--first_cutoff")
    parser.add_argument("--second_cutoff")

    args = parser.parse_args()
    main(
        input_file_path=args.input_file_path,
        output_path=args.output_path,
        first_cutoff=float(args.first_cutoff),
        second_cutoff=float(args.second_cutoff),
    )
