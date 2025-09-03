rule mkref:
    """Run cellranger mkref to add construct to reference genome"""
    input:
        fasta = join(RESULTSDIR,"custom_reference","genome.fa"),
        gtf = join(RESULTSDIR,"custom_reference","genes.gtf")   
    output:
        ref = join(RESULTSDIR,"custom_reference",config['custom_genome'])
    envmodules:
        TOOLS["cellranger"]
    shell:
        """
        TEMP=/lscratch/$SLURM_JOB_ID

        cellranger mkref --genome={output.ref} \
            --fasta={input.fasta} \
            --genes={input.gtf} \
            --memgb=48 --nthreads=12 \
            --localcores=$SLURM_CPUS_PER_TASK \
            --localmem=92 \
            --jobmode=slurm --maxjobs=20
        """
