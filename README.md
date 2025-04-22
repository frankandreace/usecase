## How to handle omics data 

Data is either PROVIDER generated (from in-house tissue) or HOSPITAL generated (no tissue?).  
It is crucial to have an associated metadata tracking table.
A database (simple SQLite, manageable with Python) might be overkill for the beginning, but it would ensure that different users can access different kind of data and not erroneously modify or delete the data. There is overhead in maintaining it and make people know how to query. Best future-proof solution.
A CSV is an easier solution fot the beginning but it seems to me a too big single point of failure.
The goal is to register for each sample useful information to reconstruct WHERE, WHEN, WHO and other specifics for ANALYSIS (batch effect) or QC.
| SAMPLE_ID | SOURCE_PROVIDER | DATA_TYPE | PLATFORM | SEQUENCING RUN | PATIENT_ID | TISSUE_SPECIFICS | CREATION_DATE | PATIENT_SEX | OTHER_PATIENT_DATA ... |
| --------- | ------- |------- |------- |------- |------- |------- |------- |------- |------- |
| HG001  | HOSPITAL    | WGS | ILLUMINA NEXTSEQ | WGS001 | 1  |  COLON   |   25/12/2024  |   M  |  ...   |
| PR002  | PRIVATE_IN_HOUSE    |  BULK_RNA_SEQ | ILLUMINA NOVASEQ | RNS001 | 2   |  COLON   |   12/01/2025  |   F  |  ...   |

### Manage Large Files
Using a VCF example
1. Sort the VCF by position (using BCFTOOLS or PICARD (or awk));
2. Split the VCF by chromosome
    * It is not essential, it depends on the actual needs (if VOI are in a subset of chromosomes);
    * Files are more tractable as they are smaller -> lower memory footprint for tools that load file in memory;
    * Faster lookups;
3. Block gzip it (bgzip) if not done in the previous steps;
4. Index the VFC(s) with TABIX for faster region retrieval;

For other sequencing data
* BAM or CRAM files instead of FASTQ? Not a fan, it might offer good compression (Illumina proprietary compression is based on alignment to a reference) and it is one step less when running the pipeline but the fact that it is aligned to a particular reference makes recalculations more costly in case you want to change the reference (see changes between GRCh37, GRCh38, CHM13).
* Other data can be just compressed, if not in binary format already.

### Future proof computational pipelines
1. Use a WORKFLOW MANAGER: either Snakemake or Nextflow:
    * The goal is to ensure reproducibility of the work that is being done so it can be audited, checked by collaboratos, easily recomputed if needed.
    * Depending on the workflow needs (Snakemake takes a lot of time to compute the DAG if the workflow is complex and with tons of data, Nextflow generates a LOT of 'ghost' files and can be an issue to the system it is operating on (max files it can handle)).

2. Use VERSION CONTROL (GIT)
    * Use github repository to store the workflow rules (functions) and other useful data (constants, parameters). Not the actual data.
    * This also ensures that we can control the evolution of the pipeline, revert modifications, check who changed what.

3. Use a modular structure to ensure easiness of understanding and reusability:
    * 1 file per funtionality (eg. 1 for READS MAPPING, 1 for VARIANT CALLING, 1 for EACH QC step), with inside thr RULE() to achieve it, e.g. different functions to map long or short reads, ecc...
    * Constants and other parameters in a separated file;
    * Software and library dependencies described in a separated file, with used version;

4. Use CONTAINERS for software dependencies (there are nuances)
    * The goal is to have a consistent environment across different platforms, as it might be needed.
    * DOCKERS provide an 'offline' copy of the softwares and the libraries needed. It is usually considered better because it stores the version and does not require internet downloading on the fly (and errors related to relying on something outside your possession). They might be of large size if there are a lot of software in use.
    * (BIO)CONDA is online and runtime. You just need to specify the requirements in a YAML file (e.g. Snakemake)

5. CLEAR DOCUMENTATION  
    * of what each function does;
    * why a function is necessary;
    * what the parameters mean;
    * why the constants have been set to a certain value;
    * MODIFICATION LOGBOOK; 

6. LIMITATIONS:
    * Software dependencies might be a problem if some versions are discontinued  (e.g. version locking in Bioconda);
    * If tools are not maintained anymore (common in bioinfo, especially research tools);
    * If softwares need different version of the same dependency, there is more confusion and additional storage for Docker Images, it should be avoided at all costs;
    * Change in HPC provider can be an issue (some functionalities, if not very common, are not supported);
    * Documentations and logging are a burden, need to find the sweet spot between recording all useful info and not make people hate it.

### QC - RNA SEQ EXAMPLE
PLOT ALL METRICS for VISUAL INSPECTION and DISCUSSION with MANAGERS.

1. PRE-ALIGNMENT control of reads coming from the sequencer. 
    * GOAL: Ensure sample quality and consistency  and flag low quality samples.
    * Take out low quality reads (Q_SCORE < 20 (0.99), 30);
    * Remove adapter contamination if present;
    * Tools: **Fastqc**, fastp, trimgalore.

2. POST-ALIGNMENT control of reads after the alignment.
    * GOAL: Verify alignment quality (alginment try to align all the sequences), remove duplicates, coverage uniformity of genes, sequencing bias,...
    * Discard reads that have either:
        * poor alignment score;
        * are not uniquely mapped;
        * overalp with >1 gene;
        * Pair end do not overlap together;
    * Picard allows to verify alignment metrics, mark duplicates
    * Rseqc for PCR bias, coverage uniformity, ...
3. BIGGEST RISKS.
    1. POOR SAMPLE QUALITY -> either less power from analysis, either have to discard;
        * Perform QC before and after the alignment. Resequence if necessary.
    2. BATCH EFFECT in ORIGIN of SAMPLE, SEQUENCING RUN, CHEMESTRY, POPULATION:
        * Keep meatadat in a table, normalize count values (e.g. RNA-SEQ), add covariates to the differential analysis;
    3. Operations on data are POORLY recorded, leads to data processing errors:
        * Log all operations either with Workflow manager or in a separate table.


### WGS PREPROCESSING STEPS
1.  Remove (germline) variants that appear in healthy pop or that are flagged as common in healthy population databases (gnomAD);
2.  Check variations on several known databases and flag or remove them.
    * see dbSNP, gnomAD (see pt. 1), COSMIC (cancer), ClinVar;
3.  Vairant-effect prediction with ANNOVAR (not done but should be interesting to see).
4. T-testing with 1Kgenomes population data (no info on health status but should be generally healthy individuals. Since they are > 3k individuals, I would filter low frequency variants as some will be pathogenic).


### BULK RNA-SEQ PREPROCESSING STEPS
1. See above for QC RNA example.
2. Quantify expression of transcript and produce a GENE MATRIX (recording COUNTS)
    * counting tools are SALMON, KALLISTO, (I co-developed MUSET for reference-free counting)
3. Normalize counts (see TPM or transcript per million formula).
4. Differential expression analysis -> apply statistical measure to identify quantitative changes.