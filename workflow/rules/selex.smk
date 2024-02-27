from os import listdir, path
from snakemake.utils import min_version


# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #


# Parameters
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]
BATCH_1_URL = config["file_urls"]["batch_1"]
BATCH_2_URL = config["file_urls"]["batch_2"]
LIFTCHAIN_URL = config["file_urls"]["liftchain"]
DELTASVM__URL = config["file_urls"]["deltasvm_supplement"]

# ------------- #
# I/O           #
# ------------- #

# Reference genome
GENOME = path.join("resources/data/genome/hg19", "hg19.fa.gz")

# Liftchain
LIFTCHAIN = path.join(INSTALL_DIR, "hg19ToHg38.over.chain.gz")

# DeltaSVM supplement
DELTASVM_INSTALL = path.join(INSTALL_DIR, "41586_2021_3211_MOESM14_ESM.csv")
DELTASVM_PROCESS = path.join(PROCESS_DIR, "deltasvm.supplement.snvs.hg19.bed")

# Raw PBS data
ORIGINAL_BATCH_INSTALL = path.join(INSTALL_DIR, "GSE118725_pbs.obs_pval05.tsv")
NOVEL_BATCH_INSTALL = path.join(INSTALL_DIR, "GSE118725_pbs.novel_batch.tsv")

# RSID hg38 coords
ORIGINAL_BATCH_RSID_HG38 = path.join(PROCESS_DIR, "original_batch_rsid_hg38.bed")
# ORIGINAL_BATCH_RSID_HG38 = path.join("workflow", "misc", "original_batch_rsid_hg38.bed")

# Original batch format
ORIGINAL_BATCH_FORMAT = path.join(PROCESS_DIR, "format_batch_original.tsv")
NOVEL_BATCH_FORMAT = path.join(PROCESS_DIR, "format_batch_novel.tsv")

# Combined batches
COMBINED_BATCHES = path.join(PROCESS_DIR, "combine_batches.hg38.tsv")
# # Formatted PBS data
# BATCH_FORMAT = path.join(PROCESS_DIR, "format_batch_{batch}.tsv")

# # With oligo bounds
# BATCH_BOUNDS = path.join(PROCESS_DIR, "add_oligo_bounds_{batch}.tsv")

# # Combined batches
# COMBINED_BATCHES = path.join(PROCESS_DIR, "combine_batches.tsv")

# # Flags
# FLAG_BOUND = path.join(PROCESS_DIR, "flag_bound.tsv")
# FLAG_PBVAR = path.join(PROCESS_DIR, "flag_pbvar.tsv")

# # Fasta seqs
# OLIGO_FASTA = path.join(PROCESS_DIR, "oligo_seqs.fa")
# SEQUENCES = path.join(PROCESS_DIR, "sequences.tsv")

# # Final - hg19
# FINAL_OUTPUT = path.join(PROCESS_DIR, "snp-selex.combined.final.hg19.tsv")
# FINAL_FORMAT = path.join(PROCESS_DIR, "snp-selex.combined.final.hg19.bed")

# # Final for real - hg38
# LIFTUP_OUTPUT = path.join(PROCESS_DIR, "snp-selex.combined.final.hg38.tsv")

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        DELTASVM_PROCESS,
        COMBINED_BATCHES,


# rule download_liftchain:
#     message:
#         """
#         Downloads liftchain from UCSC
#         """
#     output:
#         LIFTCHAIN,
#     params:
#         url=LIFTCHAIN_URL,
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/download_deltasvm_supplement.stdout",
#         stderr="workflow/logs/download_deltasvm_supplement.stderr",
#     shell:
#         "wget {params.url} -O {output}"


rule download_deltasvm_supplement:
    message:
        """
        Downloads deltasvm supplement from Nature Genetics paper.
        """
    output:
        DELTASVM_INSTALL,
    params:
        url=DELTASVM__URL,
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/download_deltasvm_supplement.stdout",
        stderr="workflow/logs/download_deltasvm_supplement.stderr",
    shell:
        "wget {params.url} -O {output}"


rule extract_deltasvm_supplement:
    message:
        """
        Formats raw data download to bed format.
        """
    input:
        rules.download_deltasvm_supplement.output,
    output:
        DELTASVM_PROCESS,
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/extract_deltasvm_supplement.stdout",
        stderr="workflow/logs/extract_deltasvm_supplement.stderr",
    shell:
        """
        awk -F',' '(NR>1) {{split($1, oligo, "_"); print oligo[1], oligo[2], oligo[3], oligo[4], $2, $3, $4, $5}}' OFS="\t" {input} |
        vawk '{{ print $1, $2-1, $2, $3, $4, $5, $6, $7, $8 }}' > {output}
        """


rule install_datasets:
    message:
        """
        Download and upack selex data across both batches. 
        Batch_1 is original and batch_2 is novel.
        """
    output:
        batch_1=ORIGINAL_BATCH_INSTALL,
        batch_2=NOVEL_BATCH_INSTALL,
    params:
        batch_1_url=BATCH_1_URL,
        batch_2_url=BATCH_2_URL,
    log:
        stdout="workflow/logs/install_datasets.stdout",
        stderr="workflow/logs/install_datasets.stderr",
    conda:
        "../envs/snp-selex.yaml"
    shell:
        """
        curl {params.batch_1_url} -o {output.batch_1}.gz && gunzip {output.batch_1}.gz
        curl {params.batch_2_url} -o {output.batch_2}.gz && gunzip {output.batch_2}.gz
        """


# rule original_rsid_hg38:
#     message:
#         """
#         Fetches coords in hg38 for rsid list
#         """
#     input:
#         ORIGINAL_BATCH_INSTALL,
#     output:
#         ORIGINAL_BATCH_RSID_HG38,
#     params:
#         script="workflow/scripts/rsid_coord.sh",
#     log:
#         stdout="workflow/logs/original_rsid_hg38.stdout",
#         stderr="workflow/logs/original_rsid_hg38.stderr",
#     conda:
#         "../envs/snp-selex.yaml"
#     shell:
#         "sh {params.script} {input} {output}"


rule original_format:
    message:
        """
        formatting original batch
        """
    input:
        rules.install_datasets.output.batch_1,
    output:
        temp(ORIGINAL_BATCH_FORMAT),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/original_format.stdout",
        stderr="workflow/logs/original_format.stderr",
    conda:
        "../envs/snp-selex.yaml"
    script:
        "../scripts/original_batch.py"


rule novel_format:
    message:
        """
        formatting original batch
        """
    input:
        rules.install_datasets.output.batch_2,
    output:
        temp(NOVEL_BATCH_FORMAT),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/novel_format.stdout",
        stderr="workflow/logs/novel_format.stderr",
    conda:
        "../envs/snp-selex.yaml"
    script:
        "../scripts/novel_batch.py"


rule combine_batches:
    message:
        """
        aggregates both batches into summary file
        """
    input:
        batch_1=rules.original_format.output,
        batch_2=rules.novel_format.output,
    output:
        COMBINED_BATCHES,
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/combine_batches.stdout",
        stderr="workflow/logs/combine_batches.stderr",
    conda:
        "../envs/snp-selex.yaml"
    shell:
        "head -n 2 {input.batch_1} > {output}; tail -n +3 -q {input.batch_2} >> {output}"


# rule original_format:
#     message:
#         """
#         Formatting for original batch
#         """
#     input:
#         selex=rules.install_datasets.output.batch_1,
#         rsids=rules.original_rsid_hg38.output,
#     output:
#        ORIGINAL_BATCH_FORMAT,
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/original_format.stdout",
#         stderr="workflow/logs/original_format.stderr",
#     conda:
#         "../envs/snp-selex.yaml"
#     shell:
#         """
#         tail -n +2 {input} | tr ":" "_" |
#         vawk '{{split($1, id, "."); print id[1], $2, $3, $4, $9, $10, $5, $6}}' |
#         vawk '{{gsub(/\-/, "_", $2); print $0}}' |
#         vawk '{{split($2, oligo, "_"); print $1, oligo[1]"_"oligo[2]+20"_"$7"_"$8, $3, $4, $5, $6, "original"}}' > {output}
#         """

# rule format_batch_1:
#     message:
#         """
#         Columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pbs-pval, batch
#         Notes:
#         - tail removes header
#         - tr replaces :
#         - vawk first split grabs TF from oligo ID
#         - vawk gsub replace - in column 2
#         - last vawk split formats oligo coords into vid to match novel batch format
#         """
#     input:
#         rules.install_datasets.output.batch_1,
#     output:
#         temp(expand(BATCH_FORMAT, batch=1)),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/format_batch_1.stdout",
#         stderr="workflow/logs/format_batch_1.stderr",
#     conda:
#         "../envs/snp-selex.yaml"
#     shell:
#         """
#         tail -n +2 {input} | tr ":" "_" |
#         vawk '{{split($1, id, "."); print id[1], $2, $3, $4, $9, $10, $5, $6}}' |
#         vawk '{{gsub(/\-/, "_", $2); print $0}}' |
#         vawk '{{split($2, oligo, "_"); print $1, oligo[1]"_"oligo[2]+20"_"$7"_"$8, $3, $4, $5, $6, "original"}}' > {output}
#         """


# rule format_batch_2:
#     message:
#         """
#         Columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pbs-pval, batch
#         Notes:
#         - tail removes header
#         - vawk prints entire line, adds "novel" batch flag
#         """
#     input:
#         rules.install_datasets.output.batch_2,
#     output:
#         temp(expand(BATCH_FORMAT, batch=2)),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/format_batch_2.stdout",
#         stderr="workflow/logs/format_batch_2.stderr",
#     conda:
#         "../envs/snp-selex.yaml"
#     shell:
#         """
#         tail -n +2 {input} | vawk '{{print $0, "novel"}}' > {output}
#         """


# rule add_oligo_bounds:
#     message:
#         """
#         columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pval, batch, profile, length_bp, oligo_chrm, oligo_lbound, oligo_rbound
#         Note: the arithmetic on the right oligo bound...this is to ensure the variant can be inlucded in the sliding window.
#         subtract 21 from left to address 0-based indexing
#         """
#     input:
#         BATCH_FORMAT,
#     output:
#         temp(BATCH_BOUNDS),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/add_oligo_bounds_{batch}.stdout",
#         stderr="workflow/logs/add_oligo_bounds_{batch}.stderr",
#     shell:
#         """
#         vawk '{{split($2, vid, "_"); print $0, vid[1], vid[2]-21, vid[2]+20}}' {input} > {output}
#         """


# rule combine_batches:
#     message:
#         """
#         Combines data from both batches into single file.
#         """
#     input:
#         batch_1=expand(BATCH_BOUNDS, batch=1),
#         batch_2=expand(BATCH_BOUNDS, batch=2),
#     output:
#         temp(COMBINED_BATCHES),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/combine_batches.stdout",
#         stderr="workflow/logs/combine_batches.stderr",
#     shell:
#         """
#         cat {input.batch_1} {input.batch_2} > {output}
#         """


# rule flag_bound:
#     message:
#         """
#         Adds flag for whether the REF allele is bound. OBS P-value < 0.05 from gvatDB.
#         """
#     input:
#         rules.combine_batches.output,
#     output:
#         temp(FLAG_BOUND),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/flag_pbvar.stdout",
#         stderr="workflow/logs/flag_pbvar.stderr",
#     shell:
#         """
#         vawk '{{if($4 <= 0.05) {{print $0,1}} else {{print $0,0}} }}' {input} > {output}
#         """


# rule flag_pbvar:
#     message:
#         """
#         Adds flag for whether the variant in PB. PBS P-value < 0.01 from gvatDB.
#         """
#     input:
#         rules.flag_bound.output,
#     output:
#         temp(FLAG_PBVAR),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/flag_pbvar.stdout",
#         stderr="workflow/logs/flag_pbvar.stderr",
#     shell:
#         """
#         vawk '{{if($6 <= 0.01) {{print $0,1}} else {{print $0,0}} }}' {input} > {output}
#         """


# rule make_oligo_fasta:
#     message:
#         """
#         Note ID formatting for 4th column in BED
#         """
#     input:
#         selex=rules.flag_pbvar.output,
#         genome=GENOME,
#     output:
#         temp(OLIGO_FASTA),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/make_oligo_fasta.stdout",
#         stderr="workflow/logs/make_oligo_fasta.stderr",
#     shell:
#         """
#         vawk '{{print $8, $9, $10, $2}}' {input.selex} |
#         bedtools getfasta -fi {input.genome} -bed stdin -nameOnly > {output}
#         """


# rule add_oligo_sequence:
#     message:
#         """
#         Adds oligo sequence from FASTA to matrix
#         """
#     input:
#         selex=rules.flag_pbvar.output,
#         fasta=rules.make_oligo_fasta.output,
#     output:
#         main=temp(FINAL_OUTPUT),
#         seqs=temp(SEQUENCES),
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/add_oligo_sequence.stdout",
#         stderr="workflow/logs/add_oligo_sequence.stderr",
#     shell:
#         """
#         grep -v '>' {input.fasta} > {output.seqs} &&
#         paste {input.selex} {output.seqs} -d $'\t' > {output.main}
#         """


# rule final_format:
#     message:
#         """
#         Formats results to BED+
#         """
#     input:
#         rules.add_oligo_sequence.output,
#     output:
#         FINAL_FORMAT,
#     conda:
#         "../envs/snp-selex.yaml"
#     log:
#         stdout="workflow/logs/final_format.stdout",
#         stderr="workflow/logs/final_format.stderr",
#     script:
#         "../scripts/format.py"


# rule liftup_variants:
#     message:
#         """
#         Lifts up variants to hg38
#         """
#     input:
#         selex=rules.final_format.output,
#         chain=LIFTCHAIN,
#     output:
#         mapped=LIFTUP_OUTPUT,
#         unmapped=temp("unMapped"),
#     conda:
#         "../envs/liftover.yaml"
#     log:
#         stdout="workflow/logs/final_format.stdout",
#         stderr="workflow/logs/final_format.stderr",
#     shell:
#         """
#         liftOver {input.selex} {input.chain} -bedPlus=3 {output.mapped} {output.unmapped}
#         """
# #liftOver snp-selex.combined.final.hg19.bed resources/data/snp-selex/hg19ToHg38.over.chain.gz -bedPlus=3 tst.hg38.bed unMapped
