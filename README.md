# summaRygwasqc

This R script is designed to automate the processing of GWAS summary statistics, aligns SNP information against the HRC reference panel, performs QC and outputs results in a standardised format. This will take approximately 2-4 minutes depending on the size of the GWAS.  

This was inspired by Richard Anney's stata script designed to do a similar function (https://github.com/ricanney/stata/blob/master/code/s/summaryqc.ado)

# Required R Packages

library(data.table)  
library(R.utils)  
library(optparse)  
library(BiocManager)  
BiocManager::install("rtracklayer")  
library(rtracklayer)

# Full list of options

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
--plink1.9 Full path to plink1.9 (OPTIONAL only use when --PRS T) (BETA)
--prsice2R Full path and prefix for PRSice.R (OPTIONAL only use when --PRS T) (BETA)
--N (OPTIONAL - Column name for total sample size)
--outdelim (OPTIONAL - Column delimiter in output file)
--prscs (OPTIONAL - Select T if you want output restricted to SNP, A1,A2 and BETA, uncompressed in tab delim format)


# 2 Automated detection of column names
Identifies columns with the following headers or similar:  
SNP (ideally contains rsIDS)  
CHR* (Chromosome)  
BP (Base position)  
A1* (A1 - reference allele. Must correspond with the OR)  
A2 (A2 - non-reference allele)    
FRQ (Allele frequency - should correspond with A1)  
OR/BETA/Z* (Effect column)  
SE (Standard error)  
P* (SNP P-Value)  
INFO (Imputation quality score)  

Please note that if the script fails, this is most likely due to an error recognising the column IDs. Please use the appropriate flags as described above if you need to specify the column name.    

Assumes that the GWAS is tab, space or comma delimited    
  
Column detection is based on pre-defined column names. Ie, the P-value column might be named P, P_VAL, P-Value etc. It searches pre-defined matches against the GWAS column names. Columns that cannot be matched are flagged.
Columns with * are considered essential. For now CHR is essential because the script breaks up the GWAS by CHR for memory and speed.
Some GWAS are formatted so CHR and BP are in the same column, usually with the format CHR:BP. The script will automatically search for columns with this format if no CHR or BP column is found. If present, it will create a new CHR and BP column based on the information in this column. However, if the format is not CHR:BP then the script will stop.  

Future functionality may include support for adding in allele information if the format in a single column is CHR:BP:A1:A2.  

If there is no SNP column, or there are low numbers of matches (<10000 per CHR) with the HRC reference panel based on rsIDs it will cross-reference physical SNP coordinates (assuming CHR and BP exist) with the HRC and extract the matchig rsID. If there are still low numbers of SNPs matched with the HRC (<10000 per CHR) after matching on RS and BP columns, the script will end. This indicates probable build issues.  

Update V1.3 - added build check. If coordinate matches between the GWAS and HRC are low, it may be the GWAS is in an older (or newer build). This now automatically updates genomic coordinates to hg19 (from hg18) if this is found.  

# 3) QC and SNP removal  
Checks probable GWAS build (hg18 or 19, not currently compatable with hg17 or 20)  
Converts letters in allele columns to be upper case (if lower)  
  
SNPs are removed if:  
Not autosomal (only keeps chrs 1-22)  
MAF < Specified threshold
INFO < Specified threshold or > 2
Non-numeric values are found in rows where there is only expected to be numeric input (ie OR/BETA/INFO/Frequency)  
SNPs have multiple alleles in A1 or A2 (indels)  
Removes rows with any non ATGC SNP. Note this, and the previous step are theoretically different but produce the same number of SNPs being removed. This is not an error but these numbers should not be considered independent.  
Removes rows with p-value <0 or >1  
If both A1 and A2 are present, removes any SNPs where both alleles are discordant with the HRC. Order of alleles here is irrelvant (ie GWAS A1 can equal HRC A1 or A2 and vice versa)  
If only A1 (or if there is a single column with A2), this allele is checked against whether it exists in HRC A1 or A2. If it doesn't match the SNP is removed.  
Removes SNPs where the BP differs between the GWAS and HRC (after checking for build issues)  

# 4) Converts OR TO BETA (logOR)
If the GWAS effect column is OR or similar, then this is automatically converted log(OR) and renamed as beta. Columns named as BETA or Z will not be changed.  

# 5) Cleanup
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
  
# 6) Writes clean GWAS summary statistics to file. 
The cleaned summary stats file will be written using the existing GWAS filename with the suffix *summaRyQC.txt.gz. This is automatically gz compressed.
  
# 7) Look at the log file
The log file will be written using the existing GWAS filename with the suffix *summaRyQC.log
This contains information about all column detection and QC steps.
Please note the total number of SNPs after QC may not match completely with the numbers reported (ie SNPs with a MAF <1% and INFO < 0.8 will be counted twice).
The script should take 2-5 minutes depending on the size of the GWAS
Email hubbardl@cardiff.ac.uk to report errors/bugs/suggestions.

# 8) Tutorial  
We will run summaRygwasqc on the PGC2 schizophrenia summary statistics file using a MAF of 1%, an INFO threshold of 0.8 and remove the MHC.    
Rscript summaRygwasqc.V1.3.R --GWAS ckqny.scz2snpres.gz --HRC HRC.Chr --PRS F --mthres 0.01 --ithres 0.8 --MHC T

The log file produced can be seen below:

Running summaRygwasqc v1.3  
Processing GWAS file ckqny.scz2snpres.gz  
Interpreting p-value to GWAS column 9: p  
Interpreting effect size to GWAS column 7: or  
Interpreting allele 1 to GWAS column 3 : a1. This is assumed to be the effect allele that corresponds to the effect size.  
Interpreting allele 2 to GWAS column 4 : a2  
Interpreting Chromosome to GWAS column 1: hg19chrc  
Interpreting Base Position to GWAS column 5: bp  
Interpreting SNP to GWAS column 2: snpid  
Interpreting INFO or imputation quality to GWAS column 6: info  
Interpreting standard error to GWAS column 8: se  
GWAS does not have an interpretable allele frequency column  
Number of SNPs in GWAS summary statistics file before QC = 9444230  
Number of SNPs in GWAS after removing rows with unexpected character values in numeric rows OR/BETA/Z/FRQ/INFO columns = 9444230  
8624491/9444230 (91.32%) of SNPs have an rsID.   
Performing QC and cross-referencing SNP BP and allele codes with the HRC v1.1 reference panel  
Checking Genome Build of ckqny.scz2snpres.gz    
100 % of SNPs had a concordant BP with the HRC on Chromosome 1. The genome build is likely to be HG19/b37   
Performing QC for chromosome 1   
Performing QC for chromosome 2   
Performing QC for chromosome 3   
Performing QC for chromosome 4   
Performing QC for chromosome 5   
Performing QC for chromosome 6   
Removing 39408 SNPs from the MHC Region (chr6:28477797-33448354)   
Performing QC for chromosome 7   
Performing QC for chromosome 8   
Performing QC for chromosome 9   
Performing QC for chromosome 10   
Performing QC for chromosome 11   
Performing QC for chromosome 12   
Performing QC for chromosome 13   
Performing QC for chromosome 14   
Performing QC for chromosome 15   
Performing QC for chromosome 16   
Performing QC for chromosome 17   
Performing QC for chromosome 18   
Performing QC for chromosome 19   
Performing QC for chromosome 20   
Performing QC for chromosome 21   
Performing QC for chromosome 22   
Changing effect size column from OR to BETA (logOR)  
Number of SNPs not present in the HRC: 1923424  
Number of SNPs with duplicated rsIDs: 0  
Number of SNPs with non ATGC alleles: 786309  
Number of SNPs multiple alleles in a single column: 786309  
Number of SNPs removed with MAF < 1 %:  0   
Number of SNPs removed with INFO < 0.8  or > 2:  1153851   
Number of SNPs removed with p-values outside of 0 and 1:  0   
Number of SNPs removed with discordant BP to HRC:  0   
Number of SNPs removed with discordant alleles to HRC:  10461   
Number of SNPs after QC:  6809917   
Updating GWAS headers and re-ordering columns  
Old SNP column ID: snpid. New SNP column ID:SNP  
Old chromosome column ID: hg19chrc. New chromosome column ID: CHR  
Old base position column ID: bp. New SNP column ID: BP  
Old allele1 column ID: a1. New allele1 column ID: A1  
Old allele2 column ID: a2. New allele2 column ID: A2  
Old OR column ID: or. New OR column ID: BETA. Note, OR has been log transformed.  
Old standard error column ID: se. New standard error column ID: SE  
Old p-value column ID: p. New p-value column ID: P  
Old imputation quality column ID: info. New imputation quality column ID: INFO  
Writing cleaned GWAS summary statistics file to ckqny.scz2snpres.gz.summaRyQC.txt.gz  
Output can be found in this directory:  ckqny.scz2snpres_SGWQC   
gwassummaRyQC complete.  







