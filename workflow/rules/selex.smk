from snakemake.utils import min_version


# Configuration
configfile: "config/config.yaml"

# Parameters
BATCH_1_URL = config["SNP-SELEX"]["FILE_URLS"]["BATCH_1"]
BATCH_2_URL = config["SNP-SELEX"]["FILE_URLS"]["BATCH_2"]
PROFILES = [i.split("|")[1] for i in config["TFBS-SCAN"]["TARGETS"]]

# Settings
min_version("7.32.4")


rule all:
    input:
        "results/snp-selex/targets.sorted.txt",
        "results/snp-selex/snp-selex.combined.tsv",
    default_target: True


rule install_datasets:
    message:
        """
        Download and upack selex data across both batches. 
        Batch_1 is original and batch_2 is novel.
        """
    output:
        batch_1="resources/data/snp-selex/GSE118725_pbs.obs_pval05.tsv",
        batch_2="resources/data/snp-selex/GSE118725_pbs.novel_batch.tsv",
    params:
        batch_1_url=BATCH_1_URL,
        batch_2_url=BATCH_2_URL,
        batch_1_out="resources/data/snp-selex/GSE118725_pbs.obs_pval05.tsv.gz",
        batch_2_out="resources/data/snp-selex/GSE118725_pbs.novel_batch.tsv.gz",
    log:
        stdout="workflow/logs/install_datasets.stdout",
        stderr="workflow/logs/install_datasets.stderr",
    threads: 1
    shell:
        """
        curl {params.batch_1_url} -o {params.batch_1_out} && gunzip {params.batch_1_out}
        curl {params.batch_2_url} -o {params.batch_2_out} && gunzip {params.batch_2_out}
        """

rule echo_targets:
    message:
        """
        Saves target list as tab delimited file
        """
    output:
        temp("results/snp-selex/targets.txt"),
    params:
        targets=config["TFBS-SCAN"]["TARGETS"],
    log:
        stdout="workflow/logs/echo_targets.stdout",
        stderr="workflow/logs/echo_targets.stderr",
    threads: 1
    run:
        with open(output[0], 'w') as outp:
            for target in params.targets:
                info = target.split("|")
                line = f"{info[0]}\t{info[1]}\t{info[2]}\n"
                outp.write(line)

rule sort_targets:
    message:
        """
        Sorts target list under JOIN specs.
        """
    input:
        rules.echo_targets.output
    output:
        "results/snp-selex/targets.sorted.txt",
    log:
        stdout="workflow/logs/sort_targets.stdout",
        stderr="workflow/logs/sort_targets.stderr",
    threads: 1
    shell:
        """
        sort -t $'\t' -k 1b,1 {input} > {output}
        """

rule format_batch_1:
    message:
        """
        Columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pbs-pval, batch
        Notes:
        - tail removes header
        - tr replaces :
        - vawk first split grabs TF from oligo ID
        - vawk gsub replace - in column 2
        - last vawk split formats oligo coords into vid to match novel batch format 
        """
    input:
      rules.install_datasets.output.batch_1
    output:
        temp("results/snp-selex/format_batch_1.tsv"),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/format_batch_1.stdout",
        stderr="workflow/logs/format_batch_1.stderr",
    threads: 1
    shell:
        """
        tail -n +2 {input} |
        tr ":" "_" |
        vawk '{{split($1, id, "."); print id[1], $2, $3, $4, $9, $10, $5, $6}}' |
        vawk '{{gsub(/\-/, "_", $2); print $0}}' |
        vawk '{{split($2, oligo, "_"); print $1, oligo[1]"_"oligo[2]+20"_"$7"_"$8, $3, $4, $5, $6, "original"}}' > {output}
        """

rule format_batch_2:
    message:
        """
        Columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pbs-pval, batch
        Notes:
        - tail removes header
        - vawk prints entire line, adds "novel" batch flag
        """
    input:
        rules.install_datasets.output.batch_2
    output:
        temp("results/snp-selex/format_batch_2.tsv"),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/format_batch_2.stdout",
        stderr="workflow/logs/format_batch_2.stderr",
    threads: 1
    shell:
        """
        tail -n +2 {input} |
        vawk '{{print $0, "novel"}}' > {output}
        """

rule sort_batch:
    message:
        """
        Alphabetical sort on TF common name, necessary for merging laters.
        Note: tr is to ensure proper tab delimited. 
        Sort param -k 1b,1 is from recommendation from join people
        """
    input:
         "results/snp-selex/format_batch_{batch}.tsv",
    output:
        temp("results/snp-selex/sorted_batch_{batch}.tsv"),
    log:
        stdout="workflow/logs/sort_batch_{batch}.stdout",
        stderr="workflow/logs/sort_batch_{batch}.stderr",
    threads: 1
    shell:
        """
        cat {input} | tr -s '\t' | sort -t $'\t' -k 1b,1 > {output}
        """

rule add_profile_info:
    message:
        """
        columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pval, batch, profile, length_bp
        note inner merge;  to left merge add -a1 -e "." flag for join
        tr is to convert output to tab delim
        """
    input:
        batch=rules.sort_batch.output,
        targets="results/snp-selex/targets.sorted.txt",
    output:
        temp("results/snp-selex/add_profile_info_{batch}.tsv"),
    log:
        stdout="workflow/logs/add_profile_info_{batch}.stdout",
        stderr="workflow/logs/add_profile_info_{batch}.stderr",
    threads: 1
    shell:
        """
        join -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3 {input.batch} {input.targets} |
        tr -s ' ' '\t' > {output}
        """


rule filter_long_motifs:
    message:
        """
        removes motifs that have PWM > 20
        """
    input:
        selex_data=rules.add_profile_info.output,
    output:
        temp("results/snp-selex/filter_long_motifs_{batch}.tsv"),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/filter_long_motifs_{batch}.stdout",
        stderr="workflow/logs/filter_long_motifs_{batch}.stderr",
    threads: 1
    shell:
        """
        vawk '{{if($9 <= 20) print $0}}' {input} > {output}
        """


rule add_oligo_bounds:
    message:
        """
        columns will be: tf, snp, oligo_auc, oligo_pval, pbs, pval, batch, profile, length_bp, oligo_chrm, oligo_lbound, oligo_rbound
        Note: the arithmetic on the right oligo bound...this is to ensure the variant can be inlucded in the sliding window.
        """
    input:
        selex_data=rules.filter_long_motifs.output,
    output:
        "results/snp-selex/add_oligo_bounds_{batch}.tsv",
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/add_oligo_bounds_{batch}.stdout",
        stderr="workflow/logs/add_oligo_bounds_{batch}.stderr",
    threads: 1
    shell:
        """
        vawk '{{split($2, vid, "_"); print $0, vid[1], vid[2]-$9, vid[2]+$9-1}}' {input} > {output}
        """


rule combine_batches:
    message:
        """
        Combines data from both batches into single file.
        """
    input:
        batch_1="results/snp-selex/add_oligo_bounds_1.tsv",
        batch_2="results/snp-selex/add_oligo_bounds_2.tsv",
    output:
        temp("results/snp-selex/combine_batches.tsv"),
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/combine_batches.stdout",
        stderr="workflow/logs/combine_batches.stderr",
    threads: 1
    shell:
        """
        cat {input.batch_1} {input.batch_2} > {output}
        """ 

rule flag_pbvar:
    message:
        """
        Adds flag for whether the variant in PB. PBS P-value < 0.01 from gvatDB.
        """
    input:
        rules.combine_batches.output,
    output:
        "results/snp-selex/snp-selex.combined.tsv",
    conda:
        "../envs/snp-selex.yaml"
    log:
        stdout="workflow/logs/flag_pbvar.stdout",
        stderr="workflow/logs/flag_pbvar.stderr",
    threads: 1
    shell:
        """
        vawk '{{if($6 <= 0.01) {{print $0,1}} else {{print $0,0}} }}' {input} > {output}
        """