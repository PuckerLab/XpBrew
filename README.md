

# XpBrew - an end-to-end pipeline for RNA-seq data processing

XpBrew is a simple Python-based pipeline for fetching RNA-seq reads from SRA and processing them to get the expression TPM files for further downstream analyses. Following are the salient features of the pipeline:

-> Optional BUSCO analysis of the PEP file obtained from the input CDS file

-> Reattempting capability to handle network disruptions during SRA prefetch and fasterq-dump

-> Iterative removal of SRA files after fetching and kallisto quantification of every accession to free-up storage space

-> Gene expression quantification by pseudopalignment using kallisto

-> Removal of invalid RNA-seq samples

-> Optional removal of isoforms

<h2>📢 News & Updates</h2>

-> **v0.5.1 is out and is the latest stable release!**

-> **v0.5.1 major changes include**:

  -> Metadata config file is mandatory to use the XpBrew pipeline. This is done to ensure
     proper metadata collection right at the start!

  -> The metadata information is used to generate a json file that later helps in building a 
     metadata reference tracking sheet

  -> A helper script has been added to the repo to facilitate obtaining CDS file(s) from GFF3 & assembly files

-> **v0.5.0 is out and is the latest stable release!**

-> **v0.5.0 major changes include**:

  -> Xpression_collector is now XpBrew!

  -> Removes kallisto files containing greater than 40% of nans in the TPMs column from failed kallisto runs

  -> Checks for duplicate samples and retains the first occurrence of a duplicate sample during the merging step 

-> **v0.4.0 is out and is the latest stable release!**

-> **v0.4.0 major changes include**:

  -> Xpression_collect now provides support for processing multiple references at once

  -> Options provided to choose ENA-FTP based fetch of accessions by integrating kingfisher-download

  -> Improved pandas-based merge provisions for easy updates of the expression files from multiple runs


-> **v0.3.0 is out and is the latest stable release!**

-> **v0.3.0 major changes include**:

  -> The storage space available at your end will be checked 

  -> Based on the storage and the number of parallel prefetches you specify the accessions queueing up after prefetch will be determined

  -> By default the storage cutoff puts the queue to freeze if the storage exceeds 75% of your available disk space; 
     This can be tuned with the --storage_percentage flag

-> **v0.2.1. is out and is the latest stable release!**

-> **v0.2.1 major changes include**:

  -> Modified isoform removal to handle multi-line FASTA files

  -> Added compressed file processing support to multi-line FASTA file handling function

-> **v0.2 is out and is the latest stable release!**

-> **v0.2 updates include**:

  -> buscolineage flag option entry has been simplified

  -> The pipeline can be restarted now after disruptions
  
  -> SRA file fetching and kallisto runs take place in an interleaved manner improving processing time
  
  -> Better error handling and network issue tolerant fall backs in SRA file fetching
  
  -> You can now merge already exisiting TPM/ repr_TPM files with the ones produced in the current run of XpBrew

## Workflow

<img width="7224" height="8192" alt="Xpression_collector drawio" src="https://github.com/user-attachments/assets/ea2ded2f-f774-422e-9d17-5c897d5eea93" />


## Installation and dependencies

### (1) Manual installation

```
git clone https://github.com/PuckerLab/XpBrew
```
**Mandatory dependencies**
- **Tools** - kallisto, sratoolkit
- **Python libraries** - pandas (v2.3.1 or greater), matplotlib (v3.10.5 or greater)

**Optional dependencies**
- **Tools** - BUSCO, kingfisher-download 
- **Python libraries** - pyfiglet (v1.0.2 or greater), rich (v15.0.0 or greater)

### (2) Docker installation

```
docker pull shakunthalan/xpbrew:latest
```

It is recommended to run the docker image as a user and not root:

```
 docker run --rm -u $(id -u) -v /path/to/data:/path/to/data shakunthalan/xpbrew:latest
```

**Note:** If you are using a docker image of the tool, you need not specify the dependencies' full paths while running the tool; 
The dependencies are built-in in the docker image and that makes it simple to use the tool across systems without the need for local installations.

### (3) Installation in a conda environment

This method of installation installs all the dependencies in a conda environment using the environment.yml file in this repository

```
git clone https://github.com/PuckerLab/XpBrew

cd XpBrew

conda env create -f environment.yml

conda activate XpBrew
```

If you have a fresh install of conda, please follow the below steps:

```
git clone https://github.com/PuckerLab/XpBrew

cd XpBrew

conda update -n base conda-libmamba-solver libmamba libmambapy

conda env create -f environment.yml

conda activate XpBrew
```
If your package metadata collection or environment solving steps are taking long, follow the below steps:

```
git clone https://github.com/PuckerLab/XpBrew

cd XpBrew

conda update -n base conda-libmamba-solver libmamba libmambapy

#Check for the solver
conda config --show solver

#If it shows solver: classic
#Do the following:

conda config --set solver libmamba

#Check again. It should show:

solver: libmamba

#Then continue with the installation as follows:

conda env create -f environment.yml

conda activate XpBrew
```

## Running the pipeline

```
Usage:

python3 XpBrew.py --input_config <CONFIG_FILE>
                               --sra <TXT_FILE_WITH_SRA_ACCESSIONS_LIST>
                               --out <OUTPUT_DIR>
MANDATORY:

--input_config         STR     Full path to tab separated config file specifying the reference name,
                               full path to CDS file, and full path to GFF file

--metadata_config	   STR     Full path to tab separated metadata config file;
                               More details provided further below

--sample_name          STR     Genus_species name of the species being processed


--sra                  STR     Full path to TXT file with one SRA accession per line


--out                  STR     Full path to output directory


OPTIONAL:

--sample_name          STR     Name of the species; default is sample

--fetch				   STR	   Method of fetching the accessions; sratoolkit | kingfisher;
                               default is sratoolkit

--fetch_fallback	   STR	   Option to fallback to sratoolkit in case kingfisher is chosen as
                               the fetch method; yes | no default is no

--kingfisher		   STR	   Full path to kingfisher executable

--gff_config           STR     Full path to GFF config file specifying the reference name in
                               the first column followed by the different GFF attributes -
                               child_attribute, child_parent_linker, parent_attribute;
                               default values in the config file are:
                               default	ID	Parent	ID

--min_sra_file_size    INT     Minimum file size cutoff in MB to check prefetched SRA
                               file sizes to cath sralite files; default cutoff is 1MB

--annotation_qc        STR     yes or no for BUSCO-based QC of the PEP file; default is yes

--threads              STR     Total number of cores for running the pipeline; default is 4

--parallel_prefetch    STR     Number of SRA accessions to be prefetched parallely; default is 2

--storage_percentage   FLOAT   Storage threshold percentage cutoff to adjust queueing of
                               parallel prefetched accessions according to available
                               disk space

--attempts             STR     Number of attempts at prefetching and accession from SRA;
                               default is just 1 attempt

--wait                 STR     Base wait time for sleep in case of network delays for prefetch;
                               Increases exponentially with a base of 2 for each reattempt;
                               default is 20 seconds

--remove_isoforms      STR     Optional step to remove isoforms; yes or no;
                               default is yes

--merge_tpms           STR     Full path to config file containing the full paths to the
                               filtered_tpm, and/or repr_filtered_tpms to be merged

                               If multiple references are present, provide the full path
                               to a folder containing the config files per reference

                               Every config file must end with an extension .merge_config

                               Sample config file name: tair10.merge_config

                               First column will have paths to filtered_tpm files one per line;
                               Second column will have paths to filtered_repr_tpm files,
                               one per line

--only_merge'		   STR	   Inserting this flag triggers only tpm and repr tpm file merging
                               and not the entire expression analysis run

--min                  INT     Minimum percentage expression of top 100 genes

--max                  INT     Maximum percentage expression of top 100 genes

--black                STR     SRA IDs to be removed or blacklisted in a TXT file with
                               one SRA accession ID per line

--busco                STR     Full path to BUSCO; Specify busco_docker if BUSCO is
                               installed via docker

--busco_lineage        STR     BUSCO lineage of the specific organism

--busco_version        STR     Version of BUSCO; default is v6.1.0

--container_version    STR     Container version of the BUSCO docker image;
                               default is cv1

--docker_host_path     STR     Host path to be mounted on for running the docker image
 
-docker_container_path STR     Container path of the docker image

--organism_type        STR     eukaryote or prokaryote for BUSCO analysis;
                               default is eukaryote

--kallisto             STR     Full path to the kallisto executable

--prefetch             STR     Full path to the prefetch executable

--fasterq-dump         STR     Full path to the fasterq-dump executable

--uncompressed         STR     Provide this flag if your read files are uncompressed

--clean_up             STR     Clean up the TMP folder after a successful run;
                               default is yes

```
## More details

-> The --metadata_config file is being included as a mandatory parameter to enforce proper metadata 
   collection at the start of the pipeline; Following are the detailed instructions to make this file:
   One line per reference with 9 tab separated columns in this order:
				1. Reference name - must match the reference name in the input config file
				2. Assembly version
				3. Annotation version
				4. Database access link or accession of the assembly/annotation
				5. Reference - publication or DOI of the assembly/annotation
				6. Family of the organism
				7. Order of the organism
				8. CDS processing method - how the CDS file was generated, e.g. get_cds_from_assembly.py_v0.3 Make sure to enter version number!
				9. Responsible person - name of the user running XpBrew
				Fields may contain spaces within but NOT tabs; All 9 fields must be complete for the tool to run


-> The pipeline has the ability to restart after disruptions. Simply repeat the command you used before the disruption

-> The TMP folder gets cleaned up after a successful run of the pipeline. You can retain it in case you want to take a deeper look

-> --parallel_prefetch flag determines the number of parallel prefetch processes at a time

  The optimal number of parallel prefetches cannot be defined precisely as it is influenced by the network bandwidth and the disk space available to store  .sra files

  But some performance optimizations have shown that a prefetch process uses anywhere between 30-45% of a CPU

  Based on this the pipeline runs two parallel prefetches by default with a dynamic queueing mechanism
  
  In the dynamic queueing and pushing design, once an accession's prefetch is done, it is pushed to the next step of fasterq-dump+pigz 
  and kallisto which function serially since both of these are CPU-bound tasks showing maximum performance when utilizing all available 
  cores at a time. The vacant spot in the prefetch block is now taken by the next accession in line making it a combinatorial scheme 
  utilizing parallelization and serialization as necessary

  Hence it is recommended to use the default parameter for --parallel_prefetch flag unless you are sure about abundance of storage space 
  and excellent network bandwidth

-> Alternative isoforms removal needs the GFF file to be given along with the CDS file and it is important for the CDS FASTA headers and the transcript identifiers in the GFF file to match to process without errors

-> Detailed guidelines on preparing the config file for --input_config flag are:

   The input config file is a simple tab-separated TXT file that has 3 columns in the order:

                    (i)   Reference name

					(ii)  Full path to CDS file

					(iii) Full path to GFF file (optional)

-> Detailed guidelines on preparing the config file for --gff_config flag are:

    The GFF config TXT file is a simple tab separated TXT file that:
					
          Needs 4 columns in the order

		  			(i) Reference name - this needs to be the same as the reference name in the first column of the input config file

					(ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID

					(iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
					respective parent feature like Parent - Note: base assumption by the tool is that all child levels
					have the same child-parent linker attribute fields. 
          
                    For eg., if Parent is the child-parent linker in the mRNA feature line,
                    then Parent will be the child-parent linker for all other child-level feature lines in the GFF

					(iv) parent_attribute: attribute field of the gene feature like ID

-> In case the headers in your CDS file are formatted differently when compared to the transcript identifiers in the GFF file supplied, you can use a helper called Fasta_fix.py to tackle this issue. The helper script and its detailed usage can be found in https://github.com/ShakNat/DupyliCate

-> --fetch allows you to choose the repository for fetching the accessions - ENA-FTP fetch via kingfisher or fetching from SRA via sratoolkit; While choosing the kingfisher fetch method, there is an option to activate fall-back to sratoolkit fetching in case the accession is not found in ENA with --fetch_fallback. In such cases, the run time can increase since the first trial is done with ENA-FTP and in case of a failure, sratoolkit is attempted as fallback. Hence it is first recommended to ensure presence of all accessions in your list in one of these databases and choose the fetch method accordingly.

-> To allow for easy merge and future updates to existing expression files, the --merge_tpms and --only_merge flags are provided. The difference between --merge_tpms and --only_merge is that the former can carry out accession fetching and expression analysis run and do the merging of the config specified files with the results produced in that run. Whereas with --only_merge, the fetching of accessions and expression analysis steps are skipped, the script checks if the output files sample_unfiltered.tpms.tsv, sample_unfiltered.repr.tpms.tsv, sample.tpms.tsv, sample.repr.tpms.tsv, sample.pep.fasta, sample.repr.pep.fasta, sample.repr.cds.fasta are already present in the putput folder specified and if so, activates the merge of the tpm files. Hnece it is important to note that the standalone merge can work out only if the output folder specified already has a set of output files from a previous XpBrew run.

## Helper scripts:

**(1) get_cds_from_assembly.py** - Python script to obtain CDS file from assembly and GFF files

Usage instructions of **get_cds_from_assembly.py**:

```
python3 get_cds_from_assembly.py
    --assembly <full path to file or folder of assembly files>
    --gff <full path to file or folder of gff files>
    --out <full path to output directory>

    NOTE_1: Assembly files should have extension of .fa, .fna or .fasta
			(optionally followed by .gz)
    NOTE_2: GFF files should have extension of .gff or .gff3
			(optionally followed by .gz)
    NOTE_3: An assembly file and a GFF file are paired when their names are
			identical after removing the extension (case ignored),
            e.g. species1.fasta.gz and species1.gff3
    NOTE_4: GFF files without CDS features are skipped with a warning,
			because exon-based sequences would include UTRs
    NOTE_5: GFF seqids (column 1) must match the FASTA headers (first word only).
			Transcripts on contigs that are missing in the FASTA are skipped;
			if no seqid matches at all, the species is skipped with
            an error message

    optional:
        --gff_config <full path to txt file containing the gff parameters to be used
			separated by spaces or tabs>
        Needs 4 columns in the order
        (i) base file name - same as the base name you use for the gff file
			(without extension) | all in case all the gff have the same
			attribute pattern
        (ii) child_attribute: attribute field of the mRNA or transcript
			 feature in the file like ID
        (iii) child_parent_linker: attribute field of the mRNA or transcript,
              CDS, exon features that link them with their respective parent
              feature like Parent - Note: base assumption by the script is that
              all child levels have the same child-parent linker attribute fields.
              For eg., if Parent is the child-parent linker in the mRNA feature line,
              then Parent will be the child-parent linker for all other child-level
              feature lines in the GFF
        (iv) parent_attribute: attribute field of the gene feature like ID
             Default if no config is given (or no matching row): ID, Parent, ID

        --process_pseudos <yes|no>
        yes (default): CDS of pseudogenes are extracted like all others
        no: CDS of pseudogenes are skipped

        --transsplice <warn|skip|process>
        Controls transcripts that look trans-spliced (CDS segments on different strands
        or carrying part= attributes)
        warn (default): print a warning and extract using the standard single-strand rules
        (sequence might be wrong)
        skip: print a warning and leave these transcripts out
        process: print a warning and extract them, ordering segments by their part= numbers
        and handling the strand of each segment separately; transcripts without part= numbers
        on all segments are skipped

        --contig_boundary_violation_process <yes|no>
        Controls CDS segments whose GFF coordinates run past the end of the contig.
        For instance, genes that cross the origin of circular genomes (mitochondria,
        plastids) this way, e.g. a CDS segment at 394625-395019 on a 367808 bp contig

        yes: such segments wrap around to the start of the contig (circular contig assumed)
        and a warning is printed for each
        no (default): coordinates are used as given, so such segments come out empty
        or cut short without a warning; check the summary for sequences whose length is not a multiple of 3

    Transcripts with unstranded features (strand other than + or -) are always
    skipped with a warning.
    CDS lines with several comma-separated parents are assigned to each of
    the parents.

    Output files in the output directory:
        <name>.cds.fasta - CDS sequences, one per transcript
        <name>.cds_summary.txt - counts of sequences, length not a multiple of 3,
        missing ATG start,missing terminal stop, internal stops,
        non-ACGT characters, duplicate IDs
        <name>.cds_problems.tsv - IDs and lengths of the flagged sequences
        with the issues found
        cds_summary_all_species.tsv - one summary row per species

```


## Third party tool references

-> https://github.com/ncbi/sra-tools

-> Bray N, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology 34, 525–527 (2016). https://doi.org/10.1038/nbt.3519

-> Delaney K. Sullivan, Kyung Hoi (Joseph) Min, Kristján Eldjárn Hjörleifsson, Laura Luebbert, Guillaume Holley, Lambda Moses, Johan Gustafsson, Nicolas L. Bray, Harold Pimentel, A. Sina Booeshaghi, Páll Melsted, Lior Pachter. kallisto, bustools and kb-python for quantifying bulk, single-cell and single-nucleus RNA-seq. Nature Protocols 20, 587-607 (2025). https://doi.org/10.1038/s41596-024-01057-0

-> Fredrik Tegenfeldt, Dmitry Kuznetsov, Mosè Manni, Matthew Berkeley, Evgeny M Zdobnov, Evgenia V Kriventseva, OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D516–D522, https://doi.org/10.1093/nar/gkae987

-> Woodcroft, B. J., Cunningham, M., Gans, J. D., Bolduc, B. B., & Hodgkins, S. B. (2023). Kingfisher: A utility for procurement of public sequencing data. [Computer software]. DOI: 10.5281/zenodo.10525085.


## Citation:

Natarajan, S. et al. XpBrew and PanXpresso - automatic RNA-seq processing workflow and comprehensive collection of gene expression data. 2026.08.27.747620 Preprint at https://doi.org/10.64898/2026.08.27.747620 (2026).





