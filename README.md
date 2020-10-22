# summaRygwasqc

This R script is designed to automate the processing of GWAS summary statistics en mass and produce standardised output according to the same QC parameters where possible. This was inspired by Richard Anney's stata script designed to do a similar function (https://github.com/ricanney/stata/blob/master/code/s/summaryqc.ado)

############################### REQUIRED PACKAGES

** "data.table"
** "R.utils" - this must be installed to allow data.table to read gz extensions 

*Both packages are required as GWAS are read in using gz. If your system cannot read in gz files then you will need to uncompress beforehand 


Brief order of operations:

###############################

############################### 1) PARSE GWAS AND HRC DATA USING COMMAND LINE INPUT 
*** Please run using the command like this: 

RScript --no-restore summaRygwasqc.R "name and path of gwas" "~/HRC.Chr"

*** Alternatively to run within R or RStudio uncomment lines 16 and 17 of the script and replace with your own values. 

*** The first option is the full path and name of the GWAS summary statistics file. 
*** The second is the path of the HRC data. The full name is HRC.Chr*.txt.gz. Please only include the directory path and "HRC.Chr" prefix

*** All output is stored in the folder containing the GWAS data. 

*** If you run this more than once, the existing log and other files will be automatically deleted.  
###############################

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

*** Assumes that the GWAS is tab, space or comma delimited

*** Column detection is based on pre-defined column names. Ie, the P-value column might be named P, P_VAL, P-Value etc. It searches pre-defined matches against
*** the GWAS column names. Columns that cannot be matched are flagged. 

*** Columns with * are considered essential. For now CHR is essential because the script breaks up the GWAS by CHR for memory and speed.
*** Future functionality will include the user being able to manually specify column names but this is not currently impemented.  

*** Some GWAS are formatted so CHR and BP are in the same column, usually with the format CHR:BP. 
*** The script will automatically search for columns with this format if no CHR or BP column is found. 
*** If present, it will create a new CHR and BP column based on the information in this column. However, if the format is not CHR:BP then the script will stop. 
*** Future functionality may include support for adding in allele information if the format in a single column is CHR:BP:A1:A2. 

*** If there is no SNP column, or there are low numbers of matches (<10000 per CHR) with the HRC reference panel based on rsIDs
*** it will cross-reference physical SNP coordinates (assuming CHR and BP exist) with the HRC and extract the matchig rsID. 
*** If there are still low numbers of SNPs matched with the HRC (<10000 per CHR) after matching on RS and BP columns, the script will end. 
*** This indicates probable build issues. 

###############################

############################### 3) QC AND SNP REMOVAL

*** Converts letters in allele columns to be upper case. 
*** Checks and removes non-numeric values in columns where there is only expected to be 
**  SNPs are removed if:
* Not autosomal (only keeps chrs 1-22)
* MAF < 1%
* INFO < 0.8 or > 2
* Non-numeric values are found in rows where there is only expected to be numeric input (ie OR/BETA/INFO/Frequency)

* SNPs have multiple alleles in A1 or A2 (indels)
* Removes rows with any non ATGC SNP. Note this, and the previous step are theoretically different but produce the same number of SNPs being removed. 
* This is not an error but these numbers should not be considered independent. 

* Removes rows with p-value <0 or >1 that are outside of the expected range. 

* If both A1 and A2 are present, removes any SNPs where both alleles are discordant with the HRC. Order of alleles here is irrelvant (ie GWAS A1 can equal HRC A1 or A2 and vice versa)
* If only A1 (or if there is a single column with A2), this allele is checked against whether it exists in HRC A1 or A2. If it doesn't match the SNP is removed. 

* Removes SNPs where the BP differs between the GWAS and HRC. 

###############################

############################### 4) CONVERTS OR TO BETA (logOR)

*** If the GWAS effect column is OR or similar, then this is automatically converted log(OR) and renamed as beta.
*** Columns named as BETA or Z or similar will not be changed. 

###############################
 
############################### 5) CLEANUP

*** Removes any other columns in the GWAS that are not present in 2)
*** If all columns are present, renames the headers to those below and re-orders the data into the following column order.

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

###############################

############################### 6) WRITES CLEANED GWAS TO FILE

*** The cleaned summary stats file will be written using the existing GWAS filename with the suffix *summaRyQC.txt.gz
*** This is automatically gz compressed

############################### 7) LOOK AT THE LOG FILE - NOTHING IS PRINTED TO THE SCREEN AS STANDARD

*** The log file will be written using the existing GWAS filename with the suffix *summaRyQC.log
*** This contains information about all column detection and QC steps. 
*** Please note the total number of SNPs after QC may not match completely with the numbers reported (ie SNPs with a MAF <1% and INFO < 0.8 will be counted twice).
*** Future versions may tidy this up if requested. 
*** The script should take 3-10 minutes depending on the size of the GWAS

*** Email hubbardl@cardiff.ac.uk to report errors/bugs/suggestions.  
