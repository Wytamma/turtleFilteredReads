#TODO: de novo 

configfile: "config.json"

GENOME_IDS = config["GENOME_IDS"].split()

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
    

download_genomes:
    """
    Download genomes from ncbi ftp server
    """
    output:
        expand("data/genomes/{GENOME_ID}.fa.gz", GENOME_ID = GENOME_IDS)
    params:
        mem = '1gb'
    run:
        from Bio import Entrez
        from Bio import SeqIO
        import shutil
        import requests

        Entrez.email = 'wytamma.wirth@me.com'

        with Entrez.esummary(db="assembly", id={wildcards.GENOME_ID}, report="full") as handle:
            esummary_record = Entrez.read(esummary_handle)
        
        AssemblyName = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        url = url + url.split('/')[-1] + '_genomic.fna.gz'

        print("Downloading", AssemblyName)
        r = requests.get(url, stream=True)
        with open('data/genomes/{wildcards.GENOME_ID}+".fa.gz"', 'wb') as f:
            shutil.copyfileobj(r.raw, f)
        

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

