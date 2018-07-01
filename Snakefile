#TODO: de novo 

configfile: "config.json"

GENOME_IDS = str(config["GENOME_IDS"]).split()
BLASTDBNAME = config["BLASTDBNAME"]
rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        'results/jobgraph.png',
        'rulegraph.png',
        f'data/BLASTDB/{BLASTDBNAME}.fa.nsq',
    params:
        mem = '1gb'
    

rule download_genomes:
    """
    Download genomes from ncbi ftp server
    """
    output:
        "data/genomes/{GENOME_ID}.fa.gz",
    params:
        mem = '1gb'
    run:
        from Bio import Entrez
        from Bio import SeqIO
        import shutil
        import requests

        Entrez.email = 'wytamma.wirth@me.com'

        with Entrez.esummary(db="assembly", id=wildcards.GENOME_ID, report="full") as esummary_handle:
            esummary_record = Entrez.read(esummary_handle)
        
        AssemblyName = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq'] or esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        if url[-1] == '/':
            url = url[:-1]
        url = f"{url}/{url.split('/')[-1]}_genomic.fna.gz"

        print("Downloading", AssemblyName)

        shell(f'curl -L {url} -o {output}')
        
        
rule combine_and_unzip_genomes:
    """
    Combines the downloaded genomes in prep for baslts databse
    """
    input:
        expand("data/genomes/{GENOME_ID}.fa.gz", GENOME_ID = GENOME_IDS)
    output:
        "data/BLASTDB/{BLASTDBNAME}.fa",
    params:
        mem = '1gb'
    shell:
        """
        gunzip -c {input} > {output}
        """

rule make_blast_db:
    """
    Make blast db
    """
    input:
        makeblastdb = 'tools/ncbi-magicblast-1.3.0/bin/makeblastdb',
        combined_genomes = "data/BLASTDB/{BLASTDBNAME}.fa",
    output:
        'data/BLASTDB/{BLASTDBNAME}.fa.nhd',
        'data/BLASTDB/{BLASTDBNAME}.fa.nhi',
        'data/BLASTDB/{BLASTDBNAME}.fa.nhr',
        'data/BLASTDB/{BLASTDBNAME}.fa.nin',
        'data/BLASTDB/{BLASTDBNAME}.fa.nog',
        'data/BLASTDB/{BLASTDBNAME}.fa.nsd',
        'data/BLASTDB/{BLASTDBNAME}.fa.nsi',
        'data/BLASTDB/{BLASTDBNAME}.fa.nsq',
    params:
        mem = '1gb',
    shell:
        """
        {input.makeblastdb} -in {input.combined_genomes} \
        -dbtype nucl \
        -title {wildcards.BLASTDBNAME} \
        -parse_seqids \
        -hash_index
        """


rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    params:
        mem = '1gb'
    output:
        "rulegraph.png"
    shell:
        """
        rm -f {output}
        snakemake --rulegraph | dot -Tpng > {output}
        """

rule generate_jobgraph:
    """
    Generate a rulegraph for the workflow.
    """
    params:
        mem = '1gb'
    output:
        "results/jobgraph.png"
    shell:
        """
        rm -f {output}
        snakemake --dag | dot -Tpng > {output}
        """

