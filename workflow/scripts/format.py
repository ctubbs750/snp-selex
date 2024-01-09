#!/usr/bin/python

import pandas as pd


# Snakemake parameters
INPUT = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

###
# Functions
###


def main():
    """D"""
    # Read in params and ensure dtypes
    fields = [
        "tf",
        "vid",
        "obs",
        "obs_p",
        "pbs",
        "pbs_p",
        "batch",
        "profile",
        "profile_bp",
        "oligo_chrm",
        "oligo_pos0",
        "oligo_pos1",
        "obs_bound",
        "pbSNP",
        "oligo_seq",
    ]
    dtypes = [
        str,
        str,
        float,
        float,
        float,
        float,
        str,
        str,
        int,
        str,
        int,
        int,
        int,
        int,
        str,
    ]

    # Read input
    selex = pd.read_csv(
        INPUT,
        header=None,
        sep="\t",
        engine="c",
        dtype=dict(zip(fields, dtypes)),
        names=fields,
    )

    # Get variant info from vid
    selex[["vid_chrm", "vid_pos1", "ref", "alt"]] = selex["vid"].str.split(
        "_", expand=True
    )
    selex["vid_pos0"] = selex["vid_pos1"].astype(int) - 1

    # Format into BED and discard uneeded fields
    order = [
        "vid_chrm",
        "vid_pos0",
        "vid_pos1",
        "profile",
        "oligo_seq",
        "tf",
        "profile_bp",
        "batch",
        "obs",
        "obs_p",
        "pbs",
        "pbs_p",
        "obs_bound",
        "pbSNP",
    ]
    selex = selex[order]

    # Write out
    selex.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
