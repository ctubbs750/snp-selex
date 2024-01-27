"""d"""

import pandas as pd


# Snakemake parameters
INPUT = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_selex_data(input_file, fields, dtypes):
    """Reads the selex data from the input file.

    Args:
        input_file (str): The path to the input file.
        fields (list): The list of field names.
        dtypes (list): The list of data types for the fields.

    Returns:
        DataFrame: The selex data.
    """
    return pd.read_csv(
        input_file,
        header=None,
        sep="\t",
        engine="c",
        dtype=dict(zip(fields, dtypes)),
        names=fields,
    )


def format_selex_data(selex):
    """Formats the selex data.

    Args:
        selex (DataFrame): The selex data.

    Returns:
        DataFrame: The formatted selex data.
    """
    # Get variant info from vid
    selex[["vid_chrm", "vid_pos1", "ref", "alt"]] = selex["vid"].str.split(
        "_", expand=True
    )

    # Add pos0
    selex["vid_pos0"] = selex["vid_pos1"].astype(int) - 1

    # change seq to upper just cause
    selex["oligo_seq"] = selex["oligo_seq"].str.upper()

    # Format into BED and discard uneeded fields
    order = [
        "vid_chrm",
        "vid_pos0",
        "vid_pos1",
        "ref",
        "alt",
        "oligo_seq",
        "tf",
        "batch",
        "obs",
        "obs_p",
        "pbs",
        "pbs_p",
        "obs_bound",
        "pbSNP",
    ]

    # Return bedplus
    return selex[order]


def write_selex_data(selex, output_file):
    """Writes the selex data to the output file.

    Args:
        selex (DataFrame): The selex data.
        output_file (str): The path to the output file.
    """
    selex.to_csv(output_file, sep="\t", index=False)


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
        "oligo_chrm",
        "oligo_pos0",
        "oligo_pos1",
        "obs_bound",
        "pbSNP",
        "oligo_seq",
    ]
    # Dtypes
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
        int,
        int,
        int,
        str,
    ]

    # Read input
    selex = read_selex_data(INPUT, fields, dtypes)

    # Format data
    selex = format_selex_data(selex)

    # Write out
    write_selex_data(selex, OUTPUT)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
