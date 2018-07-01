#TODO: de novo 

configfile: "config.json"

GENOME_IDS = str(config["GENOME_IDS"]).split()

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        'results/jobgraph.png',
        'rulegraph.png',
        expand("data/genomes/{GENOME_ID}.fa.gz", GENOME_ID = GENOME_IDS)
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
        "data/combined_genomes/{config["BLASTDBNAME"]}.fa",
    params:
        mem = '4gb'
    shell:
        """
        gunzip -c {input} > {output}
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

