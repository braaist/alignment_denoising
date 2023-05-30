import os
from Bio import AlignIO

names = [x.split(".")[0] for x in os.listdir(f"/data/ruslan_gumerov/{config['group_name']}/{config['group_name']}_fasta_aligned/")]

rule denoised_run:
    input:
        expand("{sample}_denoised.tmp", sample=names)

rule original_run:
    input:
        expand("{sample}_original.tmp", sample=names)

rule noisy:
    input:
        "../" + config["group_name"] + "_fasta_aligned/{sample}.fasta",
    output:
        "{sample}_out.fas",
        temp("{sample}_sta.gr"),
        temp("{sample}_typ.eps"),
        temp("{sample}_idx.txt"),
    shell:
        "/opt/bin/noisy {input}"

rule phylip_converter:
    input:
        "{sample}_out.fas" if "denoised" in config["type"] else "../" + config["group_name"] + "_fasta_aligned/{sample}.fasta",
    output:
        "{sample}_{type}.phy",
    params:
        type = lambda wildcards: "denoised" if "denoised" in wildcards.type else "original",
    run:
        alignments = AlignIO.parse(input[0], "fasta")
        AlignIO.write(alignments, output[0], "phylip-sequential")
        
rule run_fastme:
    input:
        "{sample}_{type}.phy",
    output:
        "{sample}_{type}.phy_fastme_tree.nwk",
        temp("{sample}_{type}.phy_fastme_stat.txt"),
    params:
        type = lambda wildcards: "denoised" if "denoised" in wildcards.type else "original",
    shell:
##Attention! Files with encoding error will be skipped
        "/opt/bin/fastme -i {input} -p -z 5051 -T 4 || true"

rule rf_dist:
    input:
        denoised = "{sample}_{type}.phy_fastme_tree.nwk",
    output:
        "{sample}_{type}.tmp"
    params:
        type = lambda wildcards: "denoised" if "denoised" in wildcards.type else "original",
        log_file = lambda wildcards: "_denoised.log" if "denoised" in wildcards.type else "_original.log",
    shell:
        """
        /data/Soft/rf_dist_n {input} /data/ruslan_gumerov/reference_trees/{config[group_name]}_{config[taxonomy_level]}_tree.nwk &> {output} || true;
        echo {output} >> {output};
        cat {output} >> {config[group_name]}{params.log_file};
        """
