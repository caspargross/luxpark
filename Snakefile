import pathlib
import re

all_sampledirs = [ x.name for x in pathlib.Path().iterdir() \
    if x.is_dir() and x.name.startswith("Sample_") ]
all_sampledirs.sort()
all_sampleids = [ re.sub(r'^Sample_', '', x) for x in all_sampledirs ]
if 'unknown' in all_sampleids:
  all_sampleids.remove('unknown')

threads_max = 8

rule all:
  input:
    expand('Sample_{s}/{s}.ccs.aligned.bam', s = all_sampleids),
    expand('Sample_{s}/genome_results.txt', s = all_sampleids),
    expand('Sample_{s}/{s}.dv.vcf.gz', s = all_sampleids),
    expand('Sample_{s}/{s}.dv.tsv', s = all_sampleids),
    expand('Sample_{s}/{s}.pbsv.vcf', s = all_sampleids),
    expand('Sample_{s}/{s}.pbsv.tsv', s = all_sampleids),
    "multiqc_report.html"

rule minimap_align:
  input:
    ref = srcdir('ref/GRCh38.mmi'),
    ccs = 'Sample_{s}/20010E_PacBio_Pilot.{s}.consensusreadset.xml'
  output:
    aligned = 'Sample_{s}/{s}.ccs.aligned.bam',
    aligned_index = 'Sample_{s}/{s}.ccs.aligned.bam.bai'
  threads: threads_max
  log: "logs/{s}_minimap.log"
  conda: "env.yaml"
  shell:
    """
    pbmm2 align \
      -j 8 \
      --sort \
      --sort-memory 2G \
      --sort-threads 4 \
      --preset CCS \
      --sample {wildcards.s} \
      --unmapped \
      {input.ref} {input.ccs} {output.aligned} > {log} 2>&1
    samtools index {output.aligned}
    """

rule qualimap:
  input: rules.minimap_align.output.aligned
  output: "Sample_{s}/genome_results.txt"
  conda: "env.yaml"
  threads: 8
  log: "logs/{s}_qualimap.log"
  params:
    roi = "--feature-file ref/roi_GBA_GRCh38.bed"
  shell:
    """
    qualimap bamqc -bam {input} -c -nt {threads} {params.roi} -outdir Sample_{wildcards.s} -nw 2000 --java-mem-size=12G > {log} 2>&1
    """

rule call_snp_singularity:
  input: "Sample_{s}/{s}.ccs.aligned.bam"
  output: "Sample_{s}/{s}.dv.vcf.gz"
  threads: 6
  params: genome = "/mnt/share/data/genomes/GRCh38.fa"
  log: "logs/{s}_dv.log"
  shell:
    """
    set +u;
    ln -f -s {params.genome} {wildcards.f}/ref.fasta
    ln -f -s {params.genome}.fai {wildcards.f}/ref.fasta.fai
    singularity run  \
      --bind "${{PWD}}/{wildcards.f}":"/input"  \
      --bind "${{PWD}}/{wildcards.f}":"/output" \
      --bind "${{PWD}}/ref":"/ref"    \
      docker://google/deepvariant \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=PACBIO \
      --ref=/ref/GRCh38.fa \
      --reads=/input/{wildcards.s}.ccs.aligned.bam \
      --output_vcf=/output/{wildcards.s}.dv.vcf.gz \
      --num_shards={threads} \
      > {log} 2>&1
    """

rule vcfstats:
  input: rules.call_snp_singularity.output
  output: "Sample_{s}/{s}_stats_var.qcML"
  threads: 1
  log: "logs/{s}_vcfstats.log"
  shell:
    """
    /mnt/storage1/share/opt/ngs-bits-2020_06-9-g409d0101/VariantQC -in {input} -out {output}
    """

rule multiqc:
  input: 
    expand("Sample_{s}/{s}_stats_var.qcML", s = all_sampleids),
    expand("Sample_{s}/genome_results.txt", s = all_sampleids)
  output:
    "multiqc_report.html"
  threads: 1
  params:
    conf = srcdir("config_multiqc.yml")
  shell:
    """
    multiqc -f -c {params.conf} {input}
    """

rule vcf2tsv:
  input:
    rules.call_snp_singularity.output
  output:
   "Sample_{s}/{s}.dv.tsv"
  shell:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%GQ\t%DP\t%AD\t%VAF\t%PL] \n'  {input} > {output}
    """

rule pbsv_discover:
  input:
   rules.minimap_align.output.aligned
  output:
   "Sample_{s}/{s}.svsig.gz"
  conda:
    "env.yaml"
  shell:
    """
    pbsv discover {input} {output}
    """

rule pbsv:
  input:
    svsig = rules.pbsv_discover.output,
    ref = srcdir("ref/GRCh38.fa")
  output:
    "Sample_{s}/{s}.pbsv.vcf"
  conda:
    "env.yaml"
  shell:
    """
    pbsv call --ccs {input.ref} {input.svsig} {output}
    """

rule vcf2tsv_pbsv:
  input:
    rules.pbsv.output
  output:
   "Sample_{s}/{s}.pbsv.tsv"
  conda:
    "env.yaml"
  shell:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%AD\t%DP\t%SAC\t%CN] \n'  {input} > {output}
    """