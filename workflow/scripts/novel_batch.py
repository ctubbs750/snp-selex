"""Formats original batch selex data from Yan et al. 2020"""

import pandas as pd
from pyliftover import LiftOver


# Snakemake parameters
IP_SNVS = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_snvs(filepath: str):
    """Reads selex data"""
    return pd.read_csv(filepath, sep="\t", engine="c")


def main():
    """Main"""
    # Read input
    snvs = read_snvs(IP_SNVS)

    # Extract info from snp field
    snvs[
        [
            "chrm",
            "pos1",
            "ref",
            "alt",
        ]
    ] = snvs[
        "snp"
    ].str.split("_", expand=True)

    # Ensure dtype
    snvs["pos1"] = snvs["pos1"].astype(int)

    # Add null expn id
    snvs["expn_id"] = "."

    # Cache liftchain
    lo = LiftOver("hg19", "hg38")

    # Map variants
    chrms = []
    pos1s = []
    for chrm, pos1 in zip(snvs["chrm"], snvs["pos1"]):
        coord = lo.convert_coordinate(chrm, pos1)
        if len(coord) == 0:
            chrms.append("NA")
            pos1s.append("NA")
        else:
            chrms.append(coord[0][0])
            pos1s.append(coord[0][1])

    # Update df
    snvs["chrm_hg38"] = chrms
    snvs["pos1_hg38"] = pos1s

    # Make hg38 vid
    snvs["vid_hg38"] = (
        snvs["chrm_hg38"]
        + ":"
        + snvs["pos1_hg38"].astype(str)
        + "-"
        + snvs["ref"]
        + "-"
        + snvs["alt"]
        + "-"
        + snvs["tf"]
    )

    # Reorder fields
    fields = ["vid_hg38", "expn_id", "oligo_auc", "oligo_pval", "pbs", "pval"]
    snvs = snvs[fields].rename(columns={"pval": "pbs_pval"})

    # Add rsid and batch tags
    snvs["rsid"] = "."
    snvs["batch"] = "novel"

    # Write out
    snvs.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
