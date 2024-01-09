#!/usr/bin/python


# Snakemake parameters
OUTPUT= snakemake.output[0]  # type: ignore
TARGETS = snakemake.params[0] # type: ignore

###
# Functions
###


def main():
    """D"""
    with open(OUTPUT, "w") as outp:
            for target in TARGETS:
                info = target.split("|")
                line = f"{info[0]}\t{info[1]}\t{info[2]}\n"
                outp.write(line)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()