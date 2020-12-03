# summaRygwasqc

This R script is designed to automate the processing of GWAS summary statistics en mass and produce standardised output according to the same QC parameters where possible.

This was inspired by Richard Anney's stata script designed to do a similar function (https://github.com/ricanney/stata/blob/master/code/s/summaryqc.ado)

############################### REQUIRED R PACKAGES

library(data.table)

library(R.utils)     

library(optparse)

library(BiocManager) 

#BiocManager::install("rtracklayer").   

library(rtracklayer). 

############################### FULL LIST OF OPTIONS

--GWAS Path and name of GWAS Summary Statistics (REQUIRED)

--HRC Location and prefix of HRC reference panel information (REQUIRED)

--PRS Derive polygenic risk scores? (SET TO T IF YOU WANT TO CALCULATE PRS OR F IF NOT. REQUIRED)

--MHC Remove the MHC locus (Set to T if you want to remove, F if not. Default is T)

--ithres INFO score threshold (IF FLAG NOT USED DEFAULT IS 0.8)

--mthres MAF threshold (IF FLAG NOT USED DEFAULT IS 0.01)

--chr CHR column name (OPTIONAL - AUTOMATICALLY DETECTED but can be specified manually with this flag)

--bp BP column name (OPTIONAL - AUTOMATICALLY DETECTED)

--effectallele Column name for allele corresponding to or/beta/z (OPTIONAL - AUTOMATICALLY DETECTED)

--nonEA Column name for allele NOT corresponding to or/beta/z (OPTIONAL - AUTOMATICALLY DETECTED)

--beta BETA column name (OPTIONAL - AUTOMATICALLY DETECTED)

--or OR column name (OPTIONAL - AUTOMATICALLY DETECTED)

--z Z-Score column name (OPTIONAL - AUTOMATICALLY DETECTED)

--pvalue Column name for GWAS p-value (OPTIONAL - AUTOMATICALLY DETECTED)

--se Column name for standard error (e.g SE) (OPTIONAL - AUTOMATICALLY DETECTED)

--info Column name for info or imputation quality (OPTIONAL - AUTOMATICALLY DETECTED)

--snp Column name for SNP or rsID (OPTIONAL - AUTOMATICALLY DETECTED)

--frq Column name for Allele 1 frequency (OPTIONAL - AUTOMATICALLY DETECTED)

--prsice2 Full path to PRSice2 binary (OPTIONAL only use when --PRS T)

--target Full path and prefix for target bed/bim/fam (OPTIONAL only use when --PRS T)

--plink1.9 Full path to plink1.9 (OPTIONAL only use when --PRS T)

--prsice2R Full path and prefix for PRSice.R (OPTIONAL only use when --PRS T)

############################### Brief order of operations:
############################### 1) PARSE GWAS AND HRC DATA USING COMMAND LINE INPUT
*** Please run using the command like this (ASSUMING ALL DATA IS IN THE SAME DIRECTORY)
*** All output is stored in a separate subdirectory containing the GWAS file and SGWQC suffix.
*** If you run this more than once on datasets of the same name, the existing log and other files will be automatically deleted.

############################### 2 AUTOMATED DETECTION AND PROCESSING OF GWAS COLUMN NAMES
*** Identifies columns with the following headers or similar:
SNP (ideally contains rsIDS)
CHR* (Chromosome)
BP (Base position)
A1* (A1 - reference allele. Must correspond with the OR)
A2 (A2 - non-reference allele).
FRQ (Allele frequency - should correspond with A1)
OR/BETA/Z* (Effect column)
SE (Standard error)
P* (SNP P-Value)
INFO (Imputation quality score)
Assumes that the GWAS is tab, space or comma delimited
Column detection is based on pre-defined column names. Ie, the P-value column might be named P, P_VAL, P-Value etc. It searches pre-defined matches against the GWAS column names. Columns that cannot be matched are flagged.
Columns with * are considered essential. For now CHR is essential because the script breaks up the GWAS by CHR for memory and speed.
Some GWAS are formatted so CHR and BP are in the same column, usually with the format CHR:BP. The script will automatically search for columns with this format if no CHR or BP column is found. If present, it will create a new CHR and BP column based on the information in this column. However, if the format is not CHR:BP then the script will stop.
Future functionality may include support for adding in allele information if the format in a single column is CHR:BP:A1:A2.
If there is no SNP column, or there are low numbers of matches (<10000 per CHR) with the HRC reference panel based on rsIDs it will cross-reference physical SNP coordinates (assuming CHR and BP exist) with the HRC and extract the matchig rsID. If there are still low numbers of SNPs matched with the HRC (<10000 per CHR) after matching on RS and BP columns, the script will end. This indicates probable build issues.
Update V1.3 - added build check. If coordinate matches between the GWAS and HRC are low, it may be the GWAS is in an older (or newer build). This now automatically updates genomic coordinates to hg19 (from hg18) if this is found.

############################### 3) QC AND SNP REMOVAL
Checks probable GWAS build (hg18 or 19, not currently compatable with hg17 or 20)
Converts letters in allele columns to be upper case.
Checks and removes non-numeric values in columns where there is only expected to be

**SNPs are removed if:

Not autosomal (only keeps chrs 1-22)
MAF < 1%
INFO < 0.8 or > 2
Non-numeric values are found in rows where there is only expected to be numeric input (ie OR/BETA/INFO/Frequency)
SNPs have multiple alleles in A1 or A2 (indels)
Removes rows with any non ATGC SNP. Note this, and the previous step are theoretically different but produce the same number of SNPs being removed. This is not an error but these numbers should not be considered independent.

Removes rows with p-value <0 or >1 that are outside of the expected range.
If both A1 and A2 are present, removes any SNPs where both alleles are discordant with the HRC. Order of alleles here is irrelvant (ie GWAS A1 can equal HRC A1 or A2 and vice versa)
If only A1 (or if there is a single column with A2), this allele is checked against whether it exists in HRC A1 or A2. If it doesn't match the SNP is removed.


Removes SNPs where the BP differs between the GWAS and HRC (after checking for build issues)

############################### 4) CONVERTS OR TO BETA (logOR)
If the GWAS effect column is OR or similar, then this is automatically converted log(OR) and renamed as beta.
Columns named as BETA or Z or similar will not be changed.

############################### 5) CLEANUP
If all columns are present, renames the headers to those below and re-orders the data into the following column order.
Label Column_Number
SNP   1
CHR   2
BP   3
A1   4
A2   5
FRQ   6
BETA   7
SE   8
P   9
INFO  10

Removes any other columns in the GWAS that are not as described above

############################### 6) WRITES CLEANED GWAS TO FILE
*** The cleaned summary stats file will be written using the existing GWAS filename with the suffix *summaRyQC.txt.gz
*** This is automatically gz compressed

############################### 7) LOOK AT THE LOG FILE 
The log file will be written using the existing GWAS filename with the suffix *summaRyQC.log
This contains information about all column detection and QC steps.
Please note the total number of SNPs after QC may not match completely with the numbers reported (ie SNPs with a MAF <1% and INFO < 0.8 will be counted twice).
The script should take 2-5 minutes depending on the size of the GWAS
Email hubbardl@cardiff.ac.uk to report errors/bugs/suggestions.
