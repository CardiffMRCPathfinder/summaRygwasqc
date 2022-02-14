### summaRygwasqc V1.3.3
### Leon Hubbard, Richard Anney, Cardiff MRC Pathfinder WS1

### OPTIMISATIONS NEEDED:
### 1) ADD FUNCTION TO UPDATE BP TO DIFFERENT BUILD IF ALLELE CODES MATCH
### 2) Update output files for SNPs removed with the snp id and column for exclusion


### Load required libraries
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(optparse))
suppressMessages(library(BiocManager))
#BiocManager::install("rtracklayer")
suppressMessages(library(rtracklayer))
#BiocManager::install("liftOver")
#library(liftover)

option_list <- list(
  make_option("--GWAS", action="store_true", type="character", help="Location of GWAS Summary Statistics"),
  make_option("--HRC", action="store_true", type="character", help="Location and prefix of HRC reference panel information"),
  make_option("--chr", action="store_true", type="character", help="CHR column name"),
  make_option("--bp", action="store_true", type="character", help="BP column name"),
  make_option("--effectallele", action="store_true",type="character", help="Column name for allele corresponding to or/beta/z"),
  make_option("--nonEA", action="store_true",type="character", help="Column name for allele NOT corresponding to or/beta/z"),
  make_option("--beta", action="store_true", type="character", help="BETA column name"),
  make_option("--or", action="store_true", type="character", help="OR column name"),
  make_option("--z", action="store_true",type="character", help="Z-Score column name"),
  make_option("--pvalue", action="store_true",type="character", help="Column name for GWAS p-value"),
  make_option("--se", action="store_true",type="character", help="Column name for standard error (e.g SE)"),
  make_option("--info", action="store_true",type="character", help="Column name for info or imputation quality (e.g INFO)"),
  make_option("--snp", action="store_true",type="character", help="Column name for SNP or rsID (e.g SNP"),
  make_option("--frq", action="store_true",type="character", help="Column name for Allele 1 frequency (e.g FRQ)"),
  make_option("--PRS", action="store_true",type="character", help="Derive polygenic risk scores?"),
  make_option("--prsice2", action="store_true",type="character", help="Full path to PRSice2 binary"),
  make_option("--target", action="store_true",type="character", help="Full path and prefix for target bed/bim/fam"),
  make_option("--plink1.9", action="store_true",type="character", help="Full path to plink1.9"),
  make_option("--prsice2R", action="store_true",type="character", help="Full path and prefix for PRSice.R"),
  make_option("--ithres", action="store_true",type="character", help="INFO score threshold"),
  make_option("--mthres", action="store_true",type="character", help="MAF score threshold"),
  make_option("--MHC", action="store_true",type="character", help="Remove MHC"),
  make_option("--N", action="store_true",type="character", help="Column name for total sample size"),
  make_option("--outdelim", action="store_true",type="character", help="Column delimiter in output file"),
  make_option("--prscs", action="store_true",type="character", help="Select T if you want output restricted to SNP, A1,A2 and BETA, uncompressed in tab delim format"),
  make_option("--convert", action="store_true",type="character", help="Select T if you have an existing SGWQC gwas file and want to convert to prscs format ONLY"),
  make_option("--lava", action="store_true",type="character", help="Select T if you want output restricted to SNP, A1,A2 and BETA, uncompressed in tab delim format"),
  make_option("--ldsc", action="store_true",type="character", help="Select T if you want output restricted to SNP, A1,A2 and BETA, uncompressed in tab delim format")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
opt <- parse_args(parser)

GWAS.file=opt$GWAS
HRC.PREFIX=opt$HRC
chr.matches=opt$chr
bp.matches=opt$bp
snp.matches=opt$snp
allele1.matches=opt$effectallele
allele2.matches=opt$nonEA
p.val.matches=opt$pvalue
beta.matches=opt$beta
or.matches=opt$or
z.matches=opt$z
se.matches=opt$se
info.matches=opt$info
frq.matches=opt$frq
infothres=opt$ithres
mafthres=opt$mthres
mhc.rm=opt$MHC
N.matches=opt$N
prscs=opt$prscs
lava=opt$lava
ldsc=opt$ldsc



outdir=paste(gsub(".txt","",gsub(".assoc","",gsub(".gz","",GWAS.file))),"_SGWQC",sep="")

suppressWarnings(if(dir.exists(outdir)==F){
  system(paste("mkdir ",outdir,sep=""))
}
)


if(length(mhc.rm)==0){
  mhc.rm=T
} 

if(length(mhc.rm)>0 & mhc.rm=="T"){
  mhc.rm=T
} 

if(length(mhc.rm)>0 & mhc.rm=="F"){
  mhc.rm=F
} 

### PRSCS output
if(length(prscs)==0){
  prscs=F
} 

if(length(prscs)>0 & prscs=="T"){
  prscs=T
} 

if(length(prscs)>0 & prscs=="F"){
  prscs=F
} 

### Lava output
if(length(lava)==0){
  lava=F
} 

if(length(lava)>0 & lava=="T"){
  lava=T
} 

if(length(lava)>0 & lava=="F"){
  lava=F
} 

### LDSC output
if(length(ldsc)==0){
  ldsc=F
} 

if(length(ldsc)>0 & ldsc=="T"){
  ldsc=T
} 

if(length(ldsc)>0 & ldsc=="F"){
  ldsc=F
}

path = system.file(package="liftOver", "extdata", "hg18ToHg19.over.chain")
ch = import.chain("hg18ToHg19.over.chain")


# GWAS.file="AUDIT_T_CLEAN.txt.gz"
#HRC.PREFIX="HRC.Chr"
#chr.matches=NULL
#bp.matches=NULL
#allele1.matches=NULL
#snp.matches=NULL
#p.val.matches=NULL
#allele2.matches=NULL
#beta.matches=NULL
#or.matches=NULL
#z.matches=NULL
#se.matches=NULL
#info.matches=NULL
#frq.matches=NULL
#infothres=0.8
#mafthres=0.01
#mhc.rm=T
#outdir=NULL

if(length(infothres)==0){
  infothres=0.8
}

if(length(mafthres)==0){
  mafthres=0.01
}

####
logmestart=function(logmestart){write.table(paste(Sys.time(),": ",logmestart,sep=""),file=paste(outdir,"/",GWAS.file,".summaRyQC.log",sep=""),col.names=F,quote=F,row.names=F)}
logme=function(logme){write.table(paste(Sys.time(),": ",logme,sep=""),file=paste(outdir,"/",GWAS.file,".summaRyQC.log",sep=""),col.names=F,quote=F,row.names=F,append=T)}

exist=function(exist)(if(exists(exist)==T){
  get(exist)
})
### Removing previous versions of log and temporary files for the GWAS if they exist
if(file.exists(paste(outdir,"/",GWAS.file,".summaRyQC.log",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".summaRyQC.log",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.nonstandard.alleles",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.nonstandard.alleles",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.maf.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.maf.snps",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.info.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.info.snps",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.multiple.alleles",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.multiple.alleles",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.p.val.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.p.val.snps",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.nonhrc.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.nonhrc.snps",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.hrcdiscbp.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.hrcdiscbp.snps",sep=""))
}

if(file.exists(paste(outdir,"/",GWAS.file,".remove.hrcdiscalleles.snps",sep=""))==T){
  system(paste("rm ",outdir,"/",GWAS.file,".remove.hrcdiscalleles.snps",sep=""))
}

logmestart("Running summaRygwasqc v1.3.4")
cat(paste("Running summaRygwasqc v1.3.4\n"))

############## Step 1 - Auto detect column names and stop pipeline if essential columns are missing

### Start by parsing column names. 
### It is quicker and less memory intensive to just read the header row with the column names
### If there are any major issues here then the pipeline stops before reading the actual data in
logmestart(paste("Processing GWAS file ",GWAS.file,sep=""))
cat(paste("Processing GWAS file ",GWAS.file,"\n",sep=""))

#logme(paste("Removing SNPs with INFO <",infothres,sep=""))
#cat(paste("Removing SNPs with INFO <",infothres,"\n",sep=""))

#logme(paste("Removing SNPs with MAF <",mafthres,sep=""))
#cat(paste("Removing SNPs with MAF <",mafthres,"\n",sep=""))

gwas.header=as.character(readLines(GWAS.file,n=1))

### I am assuming here the GWAS files will either be tab, space, or comma delimited 
if(length(grep("\t",gwas.header,value=T))>0){
  names.header=unlist(strsplit(gwas.header,split = "\t"))
  gwas.header=as.data.frame(as.matrix(do.call(rbind,strsplit(gwas.header,split = "\t"))))
  names(gwas.header)=names.header
} else if(length(grep(" ",gwas.header,value=T))>0){
  names.header=unlist(strsplit(gwas.header,split = " "))
  gwas.header=as.data.frame(as.matrix(do.call(rbind,strsplit(gwas.header,split = " "))))
  names(gwas.header)=names.header
} else if(length(grep(",",gwas.header,value=T))>0){
  names.header=unlist(strsplit(gwas.header,split = ","))
  gwas.header=as.data.frame(as.matrix(do.call(rbind,strsplit(gwas.header,split = ","))))
  names(gwas.header)=names.header
}

### Names sanity check - some files have white space delimiter but R 
if(length(which(names.header == ""))>0){
  names.header=names.header[which(names.header != "")]
  names(gwas.header)=names.header
}

# Read in GWAS to count NA values
INPUTERROR <- FALSE

# Note that print(b) fails since b doesn't exist
suppressWarnings(tryCatch(GWAS=as.data.frame(fread(input=GWAS.file,header=T,stringsAsFactors = F,fill=F)), error = function(e) { INPUTERROR <<- TRUE}))
if(INPUTERROR==T){
  GWAS=as.data.frame(fread(input=GWAS.file,header=T,stringsAsFactors = F,fill=T))
}

## Remove columns where all values are NA

GWAS.1 <- GWAS
GWAS <- GWAS.1[, unlist(lapply(GWAS.1, function(x)!all(is.na(x))))]
get_common_cols <- function(df1, df2)  intersect(names(df1), names(df2))
common_cols <- get_common_cols(GWAS, GWAS.1)
rm(GWAS.1)
`%ni%` <- Negate(`%in%`)
removed.na.cols=names.header[names.header %ni% common_cols]
if(length(removed.na.cols)>=1){
  for (i in colnames(removed.na.cols)){
    logme(paste("Removed column: ",print(removed.na.cols)," as this was identified as having missing values for all rows",sep=""))
    cat(paste("Removed column: ",print(removed.na.cols)," as this was identified as having missing values for all rows\n",sep=""))
  }}

# update names.header
names.header=names.header[names.header %in% common_cols]

### Check format of column headers to ensure there are values for all essential columns 

### P-Value checks

p.val.labels=c("P","p","pvalue","p-value","P-value","P-Value","P_VAL","PVAL","PVALUE","GC_PVALUE","P-val","Pval","Pvalue","P_LINREG","P_LOGREG","p_value","PVAL","P.value","P_BOLT_LMM","p_wald","p_T")

if(length(p.val.matches)==0){
  p.val.matches=names.header[names.header %in% p.val.labels]
}

if(length(p.val.matches)>0){
  p.val.col=which(names.header==p.val.matches)
} else if(length(p.val.matches)==0){
  stop("GWAS does not have a p-value column")
}
p.val.text=paste("Interpreting p-value to GWAS column ",p.val.col,": ",names.header[p.val.col],sep="")
cat(paste(p.val.text,"\n",sep=""))
###

### N checks

Nlabels=c("TotalSampleSize","N","Nca","Nco","TotalN","NCAS","NCON","n","Nsum","N_analyzed","Total_N")
N.matches=names.header[names.header %in% Nlabels]

suppressWarnings(if(length(N.matches)==0 & ldsc==F){
  N.text=paste("No N Samples column requested",sep="")
} else if(length(N.matches)==0 & ldsc==T){
  N.text=paste("You requested LDSC output files, but no N column has been identified in the GWAS.")
  cat(paste(N.text,"\n",sep=""))
} else if(length(N.matches)==0 & lava==T){
  N.text=paste("You requested LAVA output files, but no N column has been identified in the GWAS.")
  cat(paste(N.text,"\n",sep=""))
} else if(length(N.matches)==1){
  N.col=which(names.header==N.matches)
  N.text=paste("Interpreting N Samples to GWAS column ",N.col,": ",names.header[N.col],sep="")
} else if(length(N.matches)==2){
  N.col1=which(names.header==N.matches[1])
  N.col2=which(names.header==N.matches[2])
  N.text=paste("Interpreting N Samples to GWAS column ",N.col1,": ",names.header[N.col1]," and ", N.col2,": ",names.header[N.col2],sep="")
  N.text1=paste("Two N columns have been identified. The max value in each column will be used and summed to give a total N.")
  cat(paste(N.text1,"\n",sep=""))
}
)

cat(paste(N.text,"\n",sep=""))

### Effect size checks

beta.labels=c("BETA","Beta","beta","LOG_OR","Effect","Effect_Beta","b","LogOR","est","stdBeta","beta_T")

if(length(beta.matches)==0){
  beta.matches=names.header[names.header %in% beta.labels]
}

if(length(or.matches)==0){
  or.labels=c("OR","or","oddsratio","OddsRatio","OR(A1)","odds_ratio")
  or.matches=names.header[names.header %in% or.labels]
}

if(length(z.matches)==0){
  z.labels=c("Z","z","Z-Score","zscore","Z-score","z_score","Z_score","Z_Score","Z_SCORE","Zscore")
  z.matches=names.header[names.header %in% z.labels]
}

if(length(beta.matches)>0){
  effect.col=which(names.header==beta.matches)
} else if(length(or.matches)>0){
  effect.col=which(names.header==or.matches)
} else if(length(z.matches)>0){
  effect.col=which(names.header==z.matches)
} else if(length(beta.matches)==0 & length(or.matches)==0 & length(z.matches)==0){
  logme("GWAS does not have an interpretable effect size column")
  stop("GWAS does not have an interpretable effect size column")
}

effect.text=paste("Interpreting effect size to GWAS column ",effect.col,": ",names.header[effect.col],sep="")
cat(paste(effect.text,"\n",sep=""))

###

### Allele label checks

allele1.labels=c("a1","A1","REF","ref","ref_allele","allele","Allele","ALLELE1","ALLELE_1","EA","EFFECT_ALLELE","Effect-allele","Effect_Allele","Allele1","allele2","Effect_allele","EFF_ALLELE","effect_allele","allele1","a_0","testedAllele")
allele2.labels=c("a2","A2","ALT","alt","alt_allele","OA","ALLELE2","ALLELE_2","OTHER_ALLELE","NON_EFFECT_ALLELE","NEA","Other-allele","Non_Effect_Allele","ALLELE0","Allele2","Other_allele","NONEFF_ALLELE","other_allele","allele2","neffect_allele","a_1","otherAllele")

if(length(allele1.matches)==0){
  allele1.matches=names.header[names.header %in% allele1.labels]
}

if(length(allele2.matches)==0){
  allele2.matches=names.header[names.header %in% allele2.labels]
}

if(length(allele1.matches)+length(allele2.matches)==2){
  allele1.cols=which(allele1.matches==names.header)
  allele2.cols=which(allele2.matches==names.header)
  
  allele.text.a1=paste("Interpreting allele 1 to GWAS column ",paste(allele1.cols,sep=":")," : ",names.header[allele1.cols],". This is assumed to be the effect allele that corresponds to the effect size.",sep="")
  allele.text.a2=paste("Interpreting allele 2 to GWAS column ",paste(allele2.cols,sep=":")," : ",names.header[allele2.cols],sep="")
  
} else if(length(allele1.matches)+length(allele2.matches)==1){
  allele1.cols=which(allele1.matches==names.header)
  allele2.cols=which(allele2.matches==names.header)
  
  a1.length=length(allele1.cols)
  a2.length=length(allele2.cols)
  
  if(length(a1.length)==1){
    allele.text.a1=paste("Interpreting allele 1 to GWAS column ",paste(allele1.cols,sep=":")," : ",names.header[allele1.cols],sep="")
  } else if (length(a2.length)==1){
    allele.text.a1=paste("Interpreting allele 1 to GWAS column ",paste(allele2.cols,sep=":")," : ",names.header[allele2.cols],sep="")
  }
} else if(length(a1.length)==0 & length(a2.length)==0){
  logme("GWAS does not have any allele columns")
  stop("GWAS does not have any allele columns")
}
cat(paste(allele.text.a1,"\n",sep=""))
cat(paste(allele.text.a2,"\n",sep=""))


### Chromosome label checks

chr.labels=c("hg19chrc","CHR","chr","Chr","chromosome","CHROMOSOME","Chromosome","chrom","CHROM")

if(length(chr.matches)==0){
  chr.matches=names.header[names.header %in% chr.labels]
}

if(length(chr.matches)>0){
  chr.col=which(names.header==chr.matches)
  chr.text=paste("Interpreting Chromosome to GWAS column ",chr.col,": ",names.header[chr.col],sep="")
} else if(length(chr.matches)==0){
  chr.text=paste("GWAS does not have an interpretable CHR column")
}
cat(paste(chr.text,"\n",sep=""))


bp.labels=c("BP","bp","baseposition","Bp","POSITION","Position_hg19","POS","pos","Pos","Position(hg19)","PosGRCh37","base_pair_location")

if(length(bp.matches)==0){
  bp.matches=names.header[names.header %in% bp.labels]
}

if(length(bp.matches)>0){
  bp.col=which(names.header==bp.matches)
  bp.text=paste("Interpreting Base Position to GWAS column ",bp.col,": ",names.header[bp.col],sep="")
} else if(length(bp.matches)==0){
  bp.text=paste("GWAS does not have an interpretable BP column")
}
cat(paste(bp.text,"\n",sep=""))

### If there are no separate chr and base position columns, this might be in a single column with CHR and BP separated by a colon.
if(length(chr.matches)==0 & length(bp.matches)==0){
  
  logme("No CHR or BP column recognised. Attempting to identify whether CHR and BP are in a single column.")
  cat("No CHR or BP column recognised. Attempting to identify whether CHR and BP are in a single column.\n")
  if(length(grep(".gz",GWAS.file))==1){
    system(paste("zcat ",GWAS.file," | head -n100 | grep : > ",outdir,"/",GWAS.file,".header",sep=""))
  } else if(length(grep(".gz",GWAS.file))==0){
    system(paste("cat ",GWAS.file," | head -n100 | grep : > ",outdir,"/",GWAS.file,".header",sep=""))
  }
  
  gwas.lite=read.table(file=paste(outdir,"/",GWAS.file,".header",sep=""),header=F,stringsAsFactors = F)
  colon.finder=table(as.vector(do.call(rbind,lapply(1:nrow(gwas.lite),function(x){
    grep(":",gwas.lite[x,])
  }))))
  if(length(colon.finder)==0){
    stop("Unable to identify a column containing CHR and BP separated by :")
  } else if(length(colon.finder)>0){
    likely.chr.bp.col=which(colon.finder == max(colon.finder))
    gwas.lite.chrbp=gwas.lite[,likely.chr.bp.col]
    chr.bp=as.data.frame(do.call(rbind,strsplit(gwas.lite.chrbp,split = ":")))
    if(length(unique(chr.bp[,1]))<=26 & length(unique(chr.bp)[,2])==nrow(chr.bp)){
      chr.col=likely.chr.bp.col
      chr.text=paste("Interpreting Chromosome to GWAS column ",likely.chr.bp.col,": ",names.header[likely.chr.bp.col],sep="")
      bp.col=likely.chr.bp.col
      bp.text=paste("Interpreting Base Position to GWAS column ",bp.col,": ",names.header[bp.col],sep="")
      logme(paste("Interpreting CHR and BP as being in the same column: ",names(gwas.header)[likely.chr.bp.col],sep=""))
      cat(paste("Interpreting CHR and BP as being in the same column: ",names(gwas.header)[likely.chr.bp.col],"\n",sep=""))
      
    }
  }
}

snp.labels=c("snpid","SNP","SNPID","snp","rsid","RSID","rsID","MARKERNAME","RS_NUMBER","RS_NUMBERS","MarkerName","marker","MARKER","Marker","Marker_Name","SNPNAME","snpname","markername","ID","variant_id")

if(length(snp.matches)==0){
  snp.matches=names.header[names.header %in% snp.labels]
}

if(length(snp.matches)>0){
  snp.col=which(names.header==snp.matches)
  snp.text=paste("Interpreting SNP to GWAS column ",snp.col,": ",names.header[snp.col],sep="")
  cat(paste(snp.text,"\n",sep=""))
  
} else if(length(snp.matches)==0){
  logme("GWAS does not have an interpretable SNP column. Will attempt to match SNP information using CHR and BP.")
  cat("GWAS does not have an interpretable SNP column. Will attempt to match SNP information using CHR and BP.\n")
  
}

info.labels=c("INFO","info","R2","r2","INFO-score","Info-score","info-score","IMPINFO","minINFO","INFO_UKB")

if(length(info.matches)==0){
  info.matches=names.header[names.header %in% info.labels]
}

if(length(info.matches)>0){
  info.col=which(names.header==info.matches)
  info.text=paste("Interpreting INFO or imputation quality to GWAS column ",info.col,": ",names.header[info.col],sep="")
} else if(length(info.matches)==0){
  info.text=paste("GWAS does not have an interpretable INFO or imputation quality column")
}
cat(paste(info.text,"\n",sep=""))


se.labels=c("SE","se","standarderror","standard_error","Standard_error","Standard_Error","STANDARD_ERROR","STANDARDERROR","StdErr","StdErrLogOR","SEBETA","stderr","se_T")
if(length(se.matches)==0){
  se.matches=names.header[names.header %in% se.labels]
}
if(length(se.matches)>0){
  se.col=which(names.header==se.matches)
  se.text=paste("Interpreting standard error to GWAS column ",se.col,": ",names.header[se.col],sep="")
} else if(length(se.matches)==0){
  se.text=paste("GWAS does not have an interpretable standard error column")
}
cat(paste(se.text,"\n",sep=""))

### Note the allele frequency column can be a pain, especially if this header includes the number of samples in the study (PGC study)
### Use grep to also search for these terms before throwing an error warning. 
### The user can also manually specifiy the number of the frequency column. 

frq.labels=c("freq","FREQ","FREQ_A1","FREQ_ref","FREQ_alt","FREQ_U","FRQ","FRQ_U","frq","frequency","F_U","EAF","Effect-allele-frequency","Freq_European_1000Genomes","eaf_hapmap_CEU","Freq1","FREQ1","Freq_A1","Freq.A1","A1FREQ","Freq_HapMap","Freq_1000g","Freq","EFF_ALLELE_FREQ","AF1","FCON","effect_allele_frequency","EAF_HRC","EAF_UKB","EAF_A1","A1FREQ","FREQA1_HRC")
if(length(frq.matches)==0){
  frq.matches=names.header[names.header %in% frq.labels]
}

pgc.freq=grep("FRQ_U",names.header,value=T)
if(length(pgc.freq)>0){
  frq.matches=pgc.freq
}

if(length(frq.matches)>0){
  frq.col=which(names.header==frq.matches)
  frq.text=paste("Interpreting allele frequency to GWAS column ",frq.col,": ",names.header[frq.col],sep="")
} else if(length(frq.matches)==0){
  frq.text=paste("GWAS does not have an interpretable allele frequency column")
}

cat(paste(frq.text,"\n",sep=""))

column.checker=rbind(exist("chr.text"),exist("bp.text"),exist("snp.text"),exist("allele.text.a1"),exist("allele.text.a2"),exist("frq.text"),exist("effect.text"),exist("se.text"),exist("p.val.text"),exist("info.text"))

if(ncol(column.checker)>1){
  warning("More than one column has been assigned the same GWAS label. This is most likely the SNP ID column")
  unique.labels=which(table(c(column.checker))==1)
  if(length(grep("SNP",names(unique.labels),value=T))>0){
    GWAS.SNP.COLS=as.data.frame(fread(input=GWAS.file,header=T,stringsAsFactors = F,select=snp.col))
    SNP.length.col1=length(grep("rs",GWAS.SNP.COLS[,1]))
    SNP.length.col2=length(grep("rs",GWAS.SNP.COLS[,2]))
    if(SNP.length.col1>SNP.length.col2){
      names(GWAS.SNP.COLS)[1]
      snp.matches=names.header[names.header %in% names(GWAS.SNP.COLS)[1]]
      
      if(length(snp.matches)>0){
        snp.col=which(names.header==snp.matches)
      } else if(length(snp.matches)==0){
        stop("GWAS does not have a SNP column")
      }
      snp.text=paste("Interpreting SNP to GWAS column ",snp.col,": ",names.header[snp.col],sep="")
    } else if(SNP.length.col2>SNP.length.col1){
      
      snp.matches=names.header[names.header %in% names(GWAS.SNP.COLS)[2]]
      
      if(length(snp.matches)>0){
        snp.col=which(names.header==snp.matches)
      } else if(length(snp.matches)==0){
        stop("GWAS does not have a SNP column")
      }
      snp.text=paste("Interpreting SNP to GWAS column ",snp.col,": ",names.header[snp.col],sep="")
    }
  } 
  
  else if(length(grep("SNP",names(unique.labels),value=T))==0){
    stop("Check GWAS to see if column names have very similar names, e.g SNP and rsID")
  }}

column.checker=rbind(exist("chr.text"),exist("bp.text"),exist("snp.text"),exist("allele.text.a1"),exist("allele.text.a2"),exist("frq.text"),exist("effect.text"),exist("se.text"),exist("p.val.text"),exist("info.text"))
logme(column.checker[,1])
#########################################################
#########################################################
#########################################################

### Stage 2 - check format of GWAS file and remove or filter SNPs based on eronious or unexpected values. 

logme(paste("Number of SNPs in GWAS summary statistics file before QC = ",nrow(GWAS),sep=""))
cat(paste("Number of SNPs in GWAS summary statistics file before QC = ",nrow(GWAS),"\n",sep=""))

if(length(N.matches)==1){
  MAX.N=max(GWAS[,N.col])
  logme(paste("Highest N= ",MAX.N,sep=""))
  cat(paste("Highest N= ",MAX.N,"\n",sep=""))
}

if(length(N.matches)==2){
  MAX.N=max(GWAS[,N.col1],na.rm = T)+max(GWAS[,N.col2],na.rm = T)
  logme(paste("Highest N= ",MAX.N,sep=""))
  cat(paste("Highest N= ",MAX.N,"\n",sep=""))
}

### Doubling back from earlier, if there is likely a single chr/bp column in the GWAS this will parse that column into separate CHR and BP columns
if(length(exist("likely.chr.bp.col"))>0){
  GWAS$CHR=do.call(rbind,strsplit(GWAS[,likely.chr.bp.col],split = ":"))[,1]
  GWAS$BP=do.call(rbind,strsplit(GWAS[,likely.chr.bp.col],split = ":"))[,2]
  
  names.header=names(GWAS)
  chr.matches=names.header[names.header %in% chr.labels]
  if(length(chr.matches)>0){
    chr.col=which(names.header==chr.matches)
    chr.text=paste("Interpreting Chromosome to GWAS column ",chr.col,": ",names.header[chr.col],sep="")
  } else if(length(chr.col)==0){
    chr.text=paste("GWAS does not have a CHR column")
    stop("GWAS does not have a CHR column")
  }
  
  bp.labels=c("BP","bp","baseposition","Bp","POSITION","Position_hg19","POS","pos","Pos","BASEPOSITION","base_position")
  bp.matches=names.header[names.header %in% bp.labels]
  if(length(bp.matches)>0){
    bp.col=which(names.header==bp.matches)
    bp.text=paste("Interpreting Base Position to GWAS column ",bp.col,": ",names.header[bp.col],sep="")
  } else if(length(bp.matches)==0){
    bp.text=paste("GWAS does not have an interpretable BP column")
    warning("GWAS does not have a BP column")
  }
}

### Check classes of expected numeric columns. 

### Pvalue class
if(class(GWAS[,p.val.col])=="character" | class(GWAS[,p.val.col])=="numeric"){
  p.val.character.vals=which(!grepl('^[0-9]',GWAS[,p.val.col]))
  logme(paste("Identified and removed ",length(p.val.character.vals)," rows from GWAS with an unexpected non-numeric value in the P-Value column",sep=""))
  cat(paste("Identified and removed ",length(p.val.character.vals)," rows from GWAS with an unexpected non-numeric value in the P-Value column\n",sep=""))
  GWAS[,p.val.col]=as.numeric(GWAS[,p.val.col])
  
  if(length(p.val.character.vals)>0){
    GWAS=GWAS[-p.val.character.vals,]
  }
}

### Frq class
if(length(frq.matches)>0){
  if(class(GWAS[,frq.col])=="character" | class(GWAS[,frq.col])=="numeric"){
    frq.character.vals=which(!grepl('^[0-9]',GWAS[,frq.col]))
    logme(paste("Identified and removed ",length(frq.character.vals)," rows from GWAS with an unexpected non-numeric value in the allele frequency column",sep=""))
    cat(paste("Identified and removed ",length(frq.character.vals)," rows from GWAS with an unexpected non-numeric value in the allele frequency column\n",sep=""))
    GWAS[,frq.col]=as.numeric(as.character(GWAS[,frq.col]))
    
    if(length(frq.character.vals)>0){
      GWAS=GWAS[-frq.character.vals,]
    }
  }
}

### Info class
if(length(info.matches)>0){
  if(class(GWAS[,info.col])=="character" | class(GWAS[,info.col])=="numeric"){
    info.character.vals=which(!grepl('^[0-9]',GWAS[,info.col]))
    logme(paste("Identified and removed ",length(info.character.vals)," rows from GWAS with an unexpected non-numeric value in the INFO column",sep=""))
    cat(paste("Identified and removed ",length(info.character.vals)," rows from GWAS with an unexpected non-numeric value in the INFO column\n",sep=""))
    GWAS[,info.col]=as.numeric(GWAS[,info.col])
    
    if(length(info.character.vals)>0){
      GWAS=GWAS[-info.character.vals,]
    }
  }
}

### Effect size class
if((length(beta.matches)>0 | length(or.matches)>0 | length(z.matches)>0) & class(GWAS[,effect.col])=="character" | class(GWAS[,effect.col])=="numeric"){
  beta.character.vals=which(is.na(GWAS[,effect.col]))
  logme(paste("Identified and removed ",length(beta.character.vals)," rows from GWAS with an unexpected non-numeric value in the effect size column",sep=""))
  cat(paste("Identified and removed ",length(beta.character.vals)," rows from GWAS with an unexpected non-numeric value in the effect size column\n",sep=""))
  GWAS[,effect.col]=as.numeric(GWAS[,effect.col])
  
  if(length(beta.character.vals)>0){
    GWAS=GWAS[-beta.character.vals,]
  }
}

if(length(beta.matches)==0 | length(or.matches)>0){
  min.or=min(GWAS[,effect.col])
  
  if(min.or<0){
    stop(paste("You have negative values in your OR column. Please check if this is mislabelled"))
  }
}        

logme(paste("Number of SNPs in GWAS after removing rows with unexpected character values in numeric rows OR/BETA/Z/FRQ/INFO columns = ",nrow(GWAS),sep=""))
cat(paste("Number of SNPs in GWAS after removing rows with unexpected character values in numeric rows OR/BETA/Z/FRQ/INFO columns = ",nrow(GWAS),"\n",sep=""))

### If the total number of snps with an rsID is less than 50% of the total number of SNPs, check this against a reference panel. 

if(length(exist("snp.col"))>0){
  rs.id=length(grep("rs",GWAS[,snp.col]))
  rs.txt=paste(rs.id,"/",nrow(GWAS)," (",round(rs.id/nrow(GWAS)*100,digits=2),"%) of SNPs have an rsID.",sep="")
} else if(length(exist("snp.col"))==0){
  rs.txt=paste("There is no SNP column. Attempting to get rsIDs by cross-checking chr/bp in the HRC.",sep="")
  nosnp=T
}

if(length(exist("snp.col")==0) & length(exist("chr.col"))==0 & length(exist("bp.col"))==0){
  logme("GWAS contains no identifiable CHR, BP or SNP column")
  stop("GWAS contains no identifiable CHR, BP or SNP column")
}

logme(rs.txt)
cat(rs.txt,"\n")
### Chromosome check
### Remove trailing "chr" from values if applicable
chr.values=unique(GWAS[,chr.col])
if(length(grep("chr",chr.values))>0){
  GWAS[,chr.col]=as.numeric(as.character(gsub("chr","",GWAS[,chr.col])))
}
### Only keep autosomal SNPs
GWAS=GWAS[GWAS[,chr.col] %in% 1:22,]


#HRC=as.data.frame(fread(input="~/Downloads/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz",header=T,stringsAsFactors = F,select=1:5))
#names(HRC)=c("CHR","BP","SNP","REF","ALT")
#lapply(1:22,function(x){
#hrc.chr=HRC[which(HRC$CHR==x),]
#gz.output=gzfile(paste("HRC.Chr",x,".txt.gz",sep=""))
#write.table(hrc.chr,gz.output,col.names=T,row.names=F,quote=F)
#})

### For speed and memory usage we are looping over by chr
### This will automatically flag GWAS that are not build hg19
### Where possible the BP and CHR will be updated to newer physical coordinates. 
### No attempt is made to update allele codes
### Where alleles differ from the HRC these will be removed. The order of A1 and A2 is ignored relative to the HRC.
### Loop over each chromosome to perform checks for base position, 

var.get=paste(grep("matches",ls(),value=T))

logme(paste("Performing QC and cross-referencing SNP BP and allele codes with the HRC v1.1 reference panel"))
cat(paste("Performing QC and cross-referencing SNP BP and allele codes with the HRC v1.1 reference panel\n"))

all.chrs.cleaned=as.data.frame(do.call(rbind,lapply(1:22,function(chrom){
  if(chrom==1){
    logme(paste("Checking Genome Build of ",GWAS.file,sep=""))
    cat(paste("Checking Genome Build of",GWAS.file,"\n"),sep="")
    
    
    GWAS_chr=GWAS[which(GWAS[,chr.col]==chrom),]
    HRC.chr=as.data.frame(fread(input=paste(HRC.PREFIX,chrom,".txt.gz",sep=""),header=T,stringsAsFactors = F))
    names(HRC.chr)=c("CHR.HRC","BP.HRC","SNP.HRC","REF.HRC","ALT.HRC")
    
    if(length(exist("nosnp"))==0){
      
      merged.chr=merge(GWAS_chr,HRC.chr,by.x=snp.matches,by.y="SNP.HRC",sort=F)
      if(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)>0.8){
        write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
        logme(paste(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)*100,"% of SNPs had a concordant BP with the HRC on Chromosome 1. The genome build is likely to be HG19/b37",sep=""))
        cat(paste(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)*100,"% of SNPs had a concordant BP with the HRC on Chromosome 1. The genome build is likely to be HG19/b37","\n"))
        assign("update.build",FALSE,envir = .GlobalEnv)
      } else if(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)<0.8){
        logme(paste(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)*100,"% of SNPs had a concordant BP with the HRC on Chromosome 1. The genome build is unlikely to be HG19/b37",sep=""))
        cat(paste(nrow(merged.chr[which(merged.chr[,names(merged.chr) %in% bp.matches] == merged.chr$BP.HRC),])/nrow(merged.chr)*100,"% of SNPs had a concordant BP with the HRC on Chromosome 1.The genome build is unlikely to be HG19/b37.","\n",sep=""))
        assign("update.build",TRUE,envir = .GlobalEnv)
      }
    }
  }
  
  
  if(update.build==T){
    if(chrom==1){
      logme(paste("Updating genome build using Liftover and chain file: hg18ToHg19.over.chain"))
      cat(paste("Updating genome build using Liftover and chain file: hg18ToHg19.over.chain","\n"))
    }
    
    GWAS_chr=GWAS[which(GWAS[,chr.col]==chrom),]
    
    GWAS_chr$group=seq(1:nrow(GWAS_chr))
    
    GETTOBED=GWAS_chr[,c(chr.col,bp.col)]
    GETTOBED=cbind(paste("chr",GETTOBED[,1],sep=""),paste(GETTOBED[,2],"-",GETTOBED[,2],sep=""))
    
    grObject <- GRanges(seqnames=GETTOBED[,1], ranges=GETTOBED[,2])
    
    # run liftOver
    results <- as.data.frame(liftOver(grObject, ch))
    
    GWAS_chr.LO=merge(GWAS_chr,results,by="group",sort=F)
    GWAS_chr.LO[,chr.col+1]=as.numeric(as.character(gsub("chr","",GWAS_chr.LO$seqnames)))
    GWAS_chr.LO[,bp.col+1]=GWAS_chr.LO$end
    logme(paste("On Chr",chrom,": ",nrow(GWAS_chr.LO),"/",nrow(GWAS_chr)," SNPs were mapped to HG19/b37 coordinates",sep=""))
    cat(paste("On Chr",chrom,": ",nrow(GWAS_chr.LO),"/",nrow(GWAS_chr)," SNPs were mapped to HG19/b37 coordinates","\n"))
    GWAS_chr=GWAS_chr.LO[,c(2:(ncol(GWAS_chr.LO)-6))]
    assign("GWAS_chr",GWAS_chr,envir = .GlobalEnv)
    
  }
  
  
  
  if(update.build==F){
    GWAS_chr=GWAS[which(GWAS[,chr.col]==chrom),]
    logme(paste("Performing QC for chromosome",chrom))
    cat(paste("Performing QC for chromosome",chrom,"\n"))
  }
  
  if(mhc.rm==T & chrom==6){
    GWAS_chr.mhc=GWAS_chr[which(GWAS_chr[,chr.col]==6 & GWAS_chr[,bp.col]>28477797 & GWAS_chr[,bp.col]<33448354),]
    
    logme(paste("Removing",nrow(GWAS_chr.mhc),"SNPs from the MHC Region (chr6:28477797-33448354)"))
    cat(paste("Removing",nrow(GWAS_chr.mhc),"SNPs from the MHC Region (chr6:28477797-33448354)","\n"))
    if(nrow(GWAS_chr)>0){
      GWAS_chr=GWAS_chr[-which(GWAS_chr[,chr.col]==6 & GWAS_chr[,bp.col]>28477797 & GWAS_chr[,bp.col]<33448354),]
    }
  }  
  
  HRC.chr=as.data.frame(fread(input=paste(HRC.PREFIX,chrom,".txt.gz",sep=""),header=T,stringsAsFactors = F))
  names(HRC.chr)=c("CHR.HRC","BP.HRC","SNP.HRC","REF.HRC","ALT.HRC")
  
  
  
  ### If GWAS has a SNP ID column, merges with the HRC using this col
  if(length(exist("nosnp"))==0){
    merged.chr=merge(GWAS_chr,HRC.chr,by.x=snp.matches,by.y="SNP.HRC",sort=F)
    names.merged=names(merged.chr)
    if(nrow(merged.chr)>=10000){
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
    }
    
    if(nrow(merged.chr)<10000){
      logme(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column. Attempting match on base position.",sep=""))
      cat(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column. Attempting match on base position.\n",sep=""))
      
      merged.chr.bp=merge(GWAS_chr,HRC.chr,by.x=bp.matches,by.y="BP.HRC",sort=F)
      merged.chr.bp$SNP=merged.chr.bp$SNP.HRC
      names.merged=names(merged.chr.bp)
      
      merged.chr=merged.chr.bp
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
      
      logme(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.",sep=""))
      cat(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base posiiton.\n",sep=""))
      
      if(length(exist("likely.chr.bp.col"))>0){
        merged.chr[,names(merged.chr) %in% snp.matches]=NULL
        snp.matches="SNP"
      }
      
      if(nrow(merged.chr)<10000){
        logme(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column and base position. Please check the formatting of these columns or genome build.",sep=""))
        cat(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column and base position. Please check the formatting of these columns or genome build.\n",sep=""))
        stop("See log file for error")
      }
    }} else if(length(exist("nosnp"))==1){
      
      ### If GWAS has no SNP col, but has a BP col, merges with the HRC using this 
      merged.chr.bp=merge(GWAS_chr,HRC.chr,by.x=bp.matches,by.y="BP.HRC",sort=F)
      merged.chr.bp$SNP=merged.chr.bp$SNP.HRC
      names.merged=names(merged.chr.bp)
      merged.chr=merged.chr.bp
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
      snp.matches="SNP"
      
      logme(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.",sep=""))
      cat(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.\n",sep=""))
    }
  
  ### Check to make sure alleles are upper case otherwise this will cause issues merging with the HRC later
  ### Check and flag SNPs with non standard alleles - not strictly essential as these will be thrown out later anyway
  ### Remove rows with multiple alleles
  
  Alleles=c("C","T","A","G")
  
  if(length(allele1.cols)==1 &  length(allele2.cols)==0){
    
    GWAS_chr[,allele1.cols[1]]=toupper(GWAS_chr[,allele1.cols])
    remove.nonstandard.alleles=as.integer(row.names(GWAS_chr[!GWAS_chr[,allele1.cols] %in% unlist(Alleles),]))
    if(length(remove.nonstandard.alleles)==0){
      remove.nonstandard.alleles=NULL
    }
    
    multiple.alleles.a1=which(nchar(GWAS_chr[,allele1.cols])>1)
    if(length(multiple.alleles.a1)==0){
      remove.multiple.alleles=NULL
    } else if(length(multiple.alleles.a1)>0){
      remove.multiple.alleles=multiple.alleles.a1
    }
    
  } else if(length(allele1.cols)==1 &  length(allele2.cols)==1){
    
    GWAS_chr[,allele1.cols]=toupper(GWAS_chr[,allele1.cols])
    GWAS_chr[,allele2.cols]=toupper(GWAS_chr[,allele2.cols])
    
    remove.nonstandard.alleles=as.integer(row.names(GWAS_chr[!GWAS_chr[,allele1.cols] %in% unlist(Alleles) | !GWAS_chr[,allele2.cols] %in% unlist(Alleles),]))
    write.table(length(remove.nonstandard.alleles),file=paste(outdir,"/",GWAS.file,".","remove.nonstandard.alleles",sep=""),col.names=F,row.names=F,quote=F,append=T)
    
    
    if(length(remove.nonstandard.alleles)==0){
      remove.nonstandard.alleles=NULL
    }
    
    multiple.alleles.a1=which(nchar(GWAS_chr[,allele1.cols])>1)
    multiple.alleles.a2=which(nchar(GWAS_chr[,allele2.cols])>1)
    remove.multiple.alleles=unique(c(multiple.alleles.a1,multiple.alleles.a2))
    write.table(length(remove.multiple.alleles),file=paste(outdir,"/",GWAS.file,".","remove.multiple.alleles",sep=""),col.names=F,row.names=F,quote=F,append=T)
  }
  
  ### Allele frequency check. Remove SNPs with a MAF < 1%. 
  if(length(frq.matches)>0){
    
    lowest.maf=min(min(GWAS_chr[,frq.col]),1-max(GWAS_chr[,frq.col]))
    if(lowest.maf<as.numeric(as.character(mafthres))){
      remove.maf.snps=which(GWAS_chr[,frq.col]<as.numeric(as.character(mafthres)) | GWAS_chr[,frq.col]>(1-as.numeric(as.character(mafthres))))
      
    } else if(lowest.maf>=mafthres){
      remove.maf.snps=NULL
    }
  } else if(length(frq.matches)==0){
    remove.maf.snps=NULL
  }
  
  write.table(length(remove.maf.snps),file=paste(outdir,"/",GWAS.file,".","remove.maf.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
  
  
  ### Info check (if applicable). Remove SNPs with an INFO < INFO Threshold
  if(length(info.matches)>0){
    remove.info.snps=which(GWAS_chr[,info.col]<as.numeric(as.character(infothres)) | GWAS_chr[,info.col]>2)
    
    if(length(remove.info.snps)==0){
      remove.info.snps=NULL
    }} else if(length(info.matches)==0){
      remove.info.snps=NULL
    }
  write.table(length(remove.info.snps),file=paste(outdir,"/",GWAS.file,".","remove.info.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
  
  ### P-value check. Remove SNPs with P < 0 or > 1. 
  if(length(p.val.matches)>0){
    remove.p.val.snps=which(GWAS_chr[,p.val.col]<0 | GWAS_chr[,p.val.col]>1)
    if(length(remove.p.val.snps)==0){
      remove.p.val.snps=NULL
    }} else if(length(p.val.matches)==0){
      remove.p.val.snps=NULL
    }
  write.table(length(remove.p.val.snps),file=paste(outdir,"/",GWAS.file,".","remove.p.val.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
  
  rm.snps=unique(c(remove.nonstandard.alleles,remove.maf.snps,remove.info.snps,remove.multiple.alleles,remove.p.val.snps))
  
  
  if(length(rm.snps)>0){
    GWAS_chr=GWAS_chr[-rm.snps,]
  }
  
  HRC.chr=as.data.frame(fread(input=paste(HRC.PREFIX,chrom,".txt.gz",sep=""),header=T,stringsAsFactors = F))
  names(HRC.chr)=c("CHR.HRC","BP.HRC","SNP.HRC","REF.HRC","ALT.HRC")
  
  ######## ISSUES HERE. NOT WORKIG AS INTENDED - PUTTING OUT BLANK ROWS
  ### TO DO WITH THE SNPCOL IF IT DOESNT EXIST (IE NEW SNP FROM CHR AND BP)
  
  
  ### If GWAS has a SNP ID column, merges with the HRC using this col
  if(length(exist("nosnp"))==0){
    merged.chr=merge(GWAS_chr,HRC.chr,by.x=snp.matches,by.y="SNP.HRC",sort=F)
    names.merged=names(merged.chr)
    if(nrow(merged.chr)>=10000){
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
    }
    
    if(nrow(merged.chr)<10000){
      logme(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column. Attempting match on base position.",sep=""))
      cat(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column. Attempting match on base position.\n",sep=""))
      
      merged.chr.bp=merge(GWAS_chr,HRC.chr,by.x=bp.matches,by.y="BP.HRC",sort=F)
      merged.chr.bp$SNP=merged.chr.bp$SNP.HRC
      names.merged=names(merged.chr.bp)
      
      merged.chr=merged.chr.bp
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
      
      logme(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.",sep=""))
      cat(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base posiiton.\n",sep=""))
      
      if(length(exist("likely.chr.bp.col"))>0){
        merged.chr[,names(merged.chr) %in% snp.matches]=NULL
        snp.matches="SNP"
      }
      
      if(nrow(merged.chr)<10000){
        logme(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column and base position. Please check the formatting of these columns or genome build.",sep=""))
        cat(paste("Fewer than 10000 SNPs have been mapped to the HRC on Chr",chrom," using IDs in the SNP column and base position. Please check the formatting of these columns or genome build.\n",sep=""))
        stop("See log file for error")
      }
    }} else if(length(exist("nosnp"))==1){
      
      ### If GWAS has no SNP col, but has a BP col, merges with the HRC using this 
      merged.chr.bp=merge(GWAS_chr,HRC.chr,by.x=bp.matches,by.y="BP.HRC",sort=F)
      merged.chr.bp$SNP=merged.chr.bp$SNP.HRC
      names.merged=names(merged.chr.bp)
      merged.chr=merged.chr.bp
      write.table(nrow(GWAS_chr)-nrow(merged.chr),file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
      snp.matches="SNP"
      
      logme(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.",sep=""))
      cat(paste(nrow(merged.chr)," SNPs have been mapped to the HRC on Chr",chrom," using base position.\n",sep=""))
    }
  
  if(nrow(merged.chr)<10000){
    logme("Fewer than 10000 SNPs have been mapped to the HRC. This is likely a build issue.")
    stop("Fewer than 10000 SNPs have been mapped to the HRC. This is likely a build issue.")
  }
  
  ### If the number of SNPs that match using RS IDs are low, this suggests that there are not rsIDs in the GWAS SNP column.
  ### Instead we will attempt to match on physical coordinates (base position).
  ### If the number of matches is still low then this will cause the pipeline to stop
  
  
  ### Base position check
  
  names.merged=names(merged.chr)
  
  remove.discordant.bp=which(merged.chr[,which(names.merged==bp.matches)]!=merged.chr$BP.HRC)
  base.position.discordant=length(remove.discordant.bp)
  write.table(base.position.discordant,file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscbp.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
  
  #if(length(remove.discordant.bp)>50000){
  #  logme("More than 10000 SNPs have discordant BP with the HRC. This is likely a build issue.")
  #  stop("More than 10000 SNPs have discordant BP with the HRC. This is likely a build issue.")
  #}
  
  if(length(remove.discordant.bp)>0){
    merged.chr=merged.chr[-remove.discordant.bp,]
  }
  
  ### Allele overlap check
  if(length(allele1.cols)==1 & length(allele2.cols)==0){
    allele.overlap.prop=nrow(merged.chr)-length(which(merged.chr[,which(names.merged==allele1.matches)]==merged.chr$REF.HRC | merged.chr[,which(names.merged==allele1.matches)]==merged.chr$ALT.HRC))
    
    if(allele.overlap.prop>0){
      merged.chr=merged.chr[which(merged.chr[,which(names.merged==allele1.matches)]==merged.chr$REF.HRC | merged.chr[,which(names.merged==allele1.matches)]==merged.chr$ALT.HRC),]
    }
    
  } else if(length(allele1.cols)==1 & length(allele2.cols)==1){
    
    allele.overlap.prop=nrow(merged.chr)-length(which(
      (merged.chr[,which(names.merged==allele1.matches)]==merged.chr$REF.HRC & merged.chr[,which(names.merged==allele2.matches)]==merged.chr$ALT.HRC) |
        (merged.chr[,which(names.merged==allele2.matches)]==merged.chr$REF.HRC & merged.chr[,which(names.merged==allele1.matches)]==merged.chr$ALT.HRC)))
    
    
    if(allele.overlap.prop>0){
      merged.chr=merged.chr[which(
        (merged.chr[,which(names.merged==allele1.matches)]==merged.chr$REF.HRC & merged.chr[,which(names.merged==allele2.matches)]==merged.chr$ALT.HRC) |
          (merged.chr[,which(names.merged==allele2.matches)]==merged.chr$REF.HRC & merged.chr[,which(names.merged==allele1.matches)]==merged.chr$ALT.HRC)),]
    }
  }
  
  write.table(allele.overlap.prop,file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscalleles.snps",sep=""),col.names=F,row.names=F,quote=F,append=T)
  
  ### Remove HRC colums
  merged.chr=merged.chr[,-grep("HRC",names(merged.chr))]
  names.header=names(merged.chr)
  ### Get new column numbers in case this has changed since merging
  merged.chr=merged.chr[,names(merged.chr) %in% c(names.header)]
  
  if(length(exist("nosnp"))==1){
    var.get=c(var.get,"snp.matches")
  }
  
  cols.keep=unique(as.vector(do.call(rbind,lapply(1:length(var.get),function(x){
    if(length(get(var.get))>0){
      var.get.use=get(var.get[x])
      var.get.use
    } 
  }))))
  
  merged.chr.use=merged.chr[,names(merged.chr) %in% cols.keep]
  
  ### Effect size processing - transform OR to logOR
  
  if(length(or.matches)>0){
    merged.chr.use$BETA=log(merged.chr.use[,which(names(merged.chr.use)==or.matches)])
    merged.chr.use[,which(names(merged.chr.use)==or.matches)]=NULL
  }
  merged.chr.use
})))

#####################################
#####################################

if(length(or.matches)>0){
  logme(paste("Changing effect size column from OR to BETA (logOR)"))
  cat(paste("Changing effect size column from OR to BETA (logOR)\n"))
}

if(length(exist("likely.chr.bp.col"))>0){
  snp.matches="SNP"
}

### Remove any duplicated rs numbers
if(length(exist("nosnp"))==0){
  rs.dups=which(duplicated(all.chrs.cleaned[,which(names(all.chrs.cleaned)==snp.matches)]))
  if(length(rs.dups)>0){
    all.chrs.cleaned=all.chrs.cleaned[-rs.dups,]
  }} else if(length(exist("nosnp"))==1){
    rs.dups=which(duplicated(all.chrs.cleaned[,which(names(all.chrs.cleaned)=="SNP")]))
    if(length(rs.dups)>0){
      all.chrs.cleaned=all.chrs.cleaned[-rs.dups,]
    }}

logme(paste("Number of SNPs not present in the HRC: ",abs(sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),header=F))),sep=""))
cat(paste("Number of SNPs not present in the HRC: ",abs(sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.nonhrc.snps",sep=""),header=F))),"\n",sep=""))

logme(paste("Number of SNPs with duplicated rsIDs: ",length(rs.dups),sep=""))
cat(paste("Number of SNPs with duplicated rsIDs: ",length(rs.dups),"\n",sep=""))

logme(paste("Number of SNPs with non ATGC alleles: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.nonstandard.alleles",sep=""),header=F)),sep=""))
cat(paste("Number of SNPs with non ATGC alleles: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.nonstandard.alleles",sep=""),header=F)),"\n",sep=""))

logme(paste("Number of SNPs multiple alleles in a single column: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.multiple.alleles",sep=""),header=F)),sep=""))
cat(paste("Number of SNPs multiple alleles in a single column: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.multiple.alleles",sep=""),header=F)),"\n",sep=""))

logme(paste("Number of SNPs removed with MAF <",as.numeric(as.character(mafthres))*100,"%: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.maf.snps",sep=""),header=F))))
cat(paste("Number of SNPs removed with MAF <",as.numeric(as.character(mafthres))*100,"%: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.maf.snps",sep=""),header=F)),"\n"))

logme(paste("Number of SNPs removed with INFO <",as.numeric(as.character(infothres))," or > 2: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.info.snps",sep=""),header=F))))
cat(paste("Number of SNPs removed with INFO <",as.numeric(as.character(infothres))," or > 2: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.info.snps",sep=""),header=F)),"\n"))

logme(paste("Number of SNPs removed with p-values outside of 0 and 1: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.p.val.snps",sep=""),header=F))))
cat(paste("Number of SNPs removed with p-values outside of 0 and 1: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.p.val.snps",sep=""),header=F)),"\n"))

logme(paste("Number of SNPs removed with discordant BP to HRC: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscbp.snps",sep=""),header=F))))
cat(paste("Number of SNPs removed with discordant BP to HRC: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscbp.snps",sep=""),header=F)),"\n"))

logme(paste("Number of SNPs removed with discordant alleles to HRC: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscalleles.snps",sep=""),header=F))))
cat(paste("Number of SNPs removed with discordant alleles to HRC: ",sum(read.table(file=paste(outdir,"/",GWAS.file,".","remove.hrcdiscalleles.snps",sep=""),header=F)),"\n"))

logme(paste("Number of SNPs after QC: ",nrow(all.chrs.cleaned)))
cat(paste("Number of SNPs after QC: ",nrow(all.chrs.cleaned),"\n"))

### Rename and reorder columns with standardised header names:
### SNP CHR BP A1 A2 FRQ BETA SE P INFO 

logme(paste("Updating GWAS headers and re-ordering columns"))
cat(paste("Updating GWAS headers and re-ordering columns\n"))

names.template=data.frame(c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","P","INFO"),1:10)
names(names.template)=c("Label","Col")

####### Rename 
if(length(snp.matches)>0 & length(exist("likely.chr.bp.col"))==0){
  logme(paste("Old SNP column ID: ",snp.matches,". New SNP column ID:", "SNP",sep=""))
  cat(paste("Old SNP column ID: ",snp.matches,". New SNP column ID:", "SNP\n",sep=""))
  
  snp.matches=c(snp.matches,1)
} else if(length(snp.matches)>0 & length(exist("likely.chr.bp.col"))==1){
  logme(paste("Old SNP column ID: ",snp.matches,". New SNP column ID:", "SNP",sep=""))
  cat(paste("Old SNP column ID: ",names.header[likely.chr.bp.col],". New SNP column ID:", "SNP\n",sep=""))
} else if(length(exist("nosnp"))==1){
  logme(paste("Old SNP column ID: No SNP column. New SNP column ID:", "SNP"))
  cat(paste("Old SNP column ID: No SNP column. New SNP column ID:", "SNP\n"))
  snp.matches=c("SNP",1)
}

if(length(chr.matches)>0){
  logme(paste("Old chromosome column ID: ",chr.matches,". New chromosome column ID:"," CHR",sep=""))
  cat(paste("Old chromosome column ID: ",chr.matches,". New chromosome column ID:"," CHR\n",sep=""))
  chr.matches=c(chr.matches,2)
}

if(length(bp.matches)>0){
  logme(paste("Old base position column ID: ",bp.matches,". New SNP column ID:", " BP",sep=""))
  cat(paste("Old base position column ID: ",bp.matches,". New SNP column ID:", " BP\n",sep=""))
  
  bp.matches=c(bp.matches,3)
}

if(length(allele1.matches)>0){
  logme(paste("Old allele1 column ID: ",allele1.matches,". New allele1 column ID:", " A1",sep=""))
  cat(paste("Old allele1 column ID: ",allele1.matches,". New allele1 column ID:", " A1\n",sep=""))
  
  allele1.matches=c(allele1.matches,4)
}

if(length(allele2.matches)>0){
  logme(paste("Old allele2 column ID: ",allele2.matches,". New allele2 column ID:", " A2",sep=""))
  cat(paste("Old allele2 column ID: ",allele2.matches,". New allele2 column ID:", " A2\n",sep=""))
  
  allele2.matches=c(allele2.matches,5)
}

if(length(frq.matches)>0){
  logme(paste("Old allele frequency column ID: ",frq.matches,". New allele1 column ID:", " FRQ",sep=""))
  cat(paste("Old allele frequency column ID: ",frq.matches,". New allele1 column ID:", " FRQ\n",sep=""))
  
  frq.matches=c(frq.matches,6)
}

if(length(or.matches)>0){
  logme(paste("Old OR column ID: ",or.matches,". New OR column ID:", " BETA. Note, OR has been log transformed.",sep=""))
  cat(paste("Old OR column ID: ",or.matches,". New OR column ID:", " BETA. Note, OR has been log transformed.\n",sep=""))
  
  or.matches=c(or.matches,7)
}

if(length(beta.matches)>0){
  logme(paste("Old BETA column ID: ",beta.matches,". New BETA column ID:", "BETA",sep=""))
  cat(paste("Old BETA column ID: ",beta.matches,". New BETA column ID:", "BETA\n",sep=""))
  
  beta.matches=c(beta.matches,7)
}

if(length(z.matches)>0){
  logme(paste("Old Z column ID: ",z.matches,". New Z column ID:", " BETA",sep=""))
  cat(paste("Old Z column ID: ",z.matches,". New Z column ID:", " BETA\n",sep=""))
  
  z.matches=c(z.matches,7)
}

if(length(se.matches)>0){
  logme(paste("Old standard error column ID: ",se.matches,". New standard error column ID:", " SE",sep=""))
  cat(paste("Old standard error column ID: ",se.matches,". New standard error column ID:", " SE\n",sep=""))
  
  se.matches=c(se.matches,8)
}

if(length(p.val.matches)>0){
  logme(paste("Old p-value column ID: ",p.val.matches,". New p-value column ID:", " P",sep=""))
  cat(paste("Old p-value column ID: ",p.val.matches,". New p-value column ID:", " P\n",sep=""))
  
  p.val.matches=c(p.val.matches,9)
}

if(length(info.matches)>0){
  logme(paste("Old imputation quality column ID: ",info.matches,". New imputation quality column ID:", " INFO",sep=""))
  cat(paste("Old imputation quality column ID: ",info.matches,". New imputation quality column ID:", " INFO\n",sep=""))
  
  info.matches=c(info.matches,10)
}

cols.keep=na.omit(unique(as.vector(do.call(rbind,lapply(1:length(var.get),function(x){
  if(length(get(var.get))>0){
    var.get.use=get(var.get[x])
    as.data.frame(cbind(var.get.use[1],var.get.use[2]))
  } 
})))))

names.template=merge(names.template,cols.keep,by.x="Col",by.y="V2")
col.order=c(as.character(names.template$V1))

if(length(grep(or.matches[1],col.order))>0){
  col.order[grep(or.matches[1],col.order)]="BETA"
}

all.chrs.cleaned=all.chrs.cleaned[,col.order]
names(all.chrs.cleaned)=names.template$Label

if(length(N.matches)>0){
  all.chrs.cleaned$N=MAX.N
}

logme(paste("Writing cleaned GWAS summary statistics file to ",GWAS.file,".summaRyQC.txt.gz",sep=""))
cat(paste("Writing cleaned GWAS summary statistics file to ",GWAS.file,".summaRyQC.txt.gz\n",sep=""))

gz.output=gzfile(paste(outdir,"/",GWAS.file,".summaRyQC.txt.gz",sep=""))
write.table(all.chrs.cleaned,gz.output,col.names=T,row.names=F,quote=F,sep="\t")


### Clean header col search
names.cleaned=names(all.chrs.cleaned)
SNP.COL=names.cleaned[names.cleaned %in% "SNP"]
A1.COL=names.cleaned[names.cleaned %in% "A1"]
A2.COL=names.cleaned[names.cleaned %in% "A2"]
BETA.COL=names.cleaned[names.cleaned %in% "BETA"]
if(length(N.matches)>0){
  N.COL=names.cleaned[names.cleaned %in% "N"]
}
P.COL=names.cleaned[names.cleaned %in% "P"]



if(prscs==T){
  logme(paste("Writing cleaned GWAS summary statistics in PRSCS format to ",GWAS.file,".summaRyQC.PRSCS.txt",sep=""))
  cat(paste("Writing cleaned GWAS summary statistics in PRSCS format to ",GWAS.file,".summaRyQC.PRSCS.txt\n",sep=""))
  prscs.output=paste(outdir,"/",GWAS.file,".summaRyQC.PRSCS.txt",sep="")
  write.table(all.chrs.cleaned[,c(SNP.COL,A1.COL,A2.COL,BETA.COL,P.COL)],prscs.output,col.names=T,row.names=F,quote=F,sep="\t")
}

if(lava==T){
  logme(paste("Writing cleaned GWAS summary statistics in Lava format to ",GWAS.file,".summaRyQC.lava.txt",sep=""))
  cat(paste("Writing cleaned GWAS summary statistics in Lava format to ",GWAS.file,".summaRyQC.lava.txt\n",sep=""))
  prscs.output=paste(outdir,"/",GWAS.file,".summaRyQC.lava.txt",sep="")
  write.table(all.chrs.cleaned[,c(SNP.COL,A1.COL,A2.COL,N.COL,BETA.COL,P.COL)],prscs.output,col.names=T,row.names=F,quote=F,sep=" ")
}

if(ldsc==T & length(N.matches)>0){
  logme(paste("Writing cleaned GWAS summary statistics in LDSC format to ",GWAS.file,".summaRyQC.ldsc.txt.gz",sep=""))
  cat(paste("Writing cleaned GWAS summary statistics in LDSC format to ",GWAS.file,".summaRyQC.ldsc.txt.gz\n",sep=""))
  gz.output.ldsc=gzfile(paste(outdir,"/",GWAS.file,".summaRyQC.ldsc.txt.gz",sep=""))
  write.table(all.chrs.cleaned[,c(SNP.COL,A1.COL,A2.COL,N.COL,P.COL,BETA.COL)],gz.output.ldsc,col.names=T,row.names=F,quote=F,sep=" ")
}

if(ldsc==T & length(N.matches)==0){
  logme(paste("Writing cleaned GWAS summary statistics in LDSC format to ",GWAS.file,".summaRyQC.ldsc.txt.gz",sep=""))
  cat(paste("Writing cleaned GWAS summary statistics in LDSC format to ",GWAS.file,".summaRyQC.ldsc.txt.gz\n",sep=""))
  gz.output.ldsc=gzfile(paste(outdir,"/",GWAS.file,".summaRyQC.ldsc.txt.gz",sep=""))
  write.table(all.chrs.cleaned[,c(SNP.COL,A1.COL,A2.COL,P.COL,BETA.COL)],gz.output.ldsc,col.names=T,row.names=F,quote=F,sep=" ")
}

logme(paste("Output can be found in this directory: ",outdir,sep=""))
cat(paste("Output can be found in this directory: ",outdir,"\n"))

logme(paste("gwassummaRyQC complete.",sep=""))
cat(paste("gwassummaRyQC complete.\n\n"))



### PRS Functionality
PRS=opt$PRS
if(PRS=="T"){
  logme(paste("Running PRS module",sep=""))
  cat(paste("Running PRS module\n",sep=""))
  
  plinkexe=opt$plink1.9 
  TARGET.prefix=opt$target
  prsice2=opt$prsice2
  prsice2r=opt$prsice2R
  
  logme(paste("Reading target bim file: ",TARGET.prefix,".bim",sep=""))
  cat(paste("Reading target bim file: ",TARGET.prefix,".bim\n",sep=""))
  target.bim=as.data.frame(fread(input=paste(TARGET.prefix,".bim",sep=""),header=F,stringsAsFactors = F))
  
  ### Only retain SNPs present in both discovery and target datasets
  all.chrs.cleaned.overlap=all.chrs.cleaned[all.chrs.cleaned$SNP %in% unlist(target.bim$V2),]
  
  logme(paste("Number of SNPs overlapping discovery and target datasets: ",nrow(all.chrs.cleaned.overlap),sep=""))
  cat(paste("Number of SNPs overlapping discovery and target datasets: ",nrow(all.chrs.cleaned.overlap),"\n",sep=""))
  
  if(nrow(all.chrs.cleaned.overlap)<100000){
    logme("Fewer than 100000 SNPs overlap between discovery and target datasets. This is likely a build issue.")
    stop("Fewer than 100000 SNPs overlap between discovery and target datasets. This is likely a build issue.")
  }
  
  gz.output=gzfile(paste(outdir,"/",GWAS.file,".",TARGET.prefix,".summaRyQC.txt.gz",sep=""))
  write.table(all.chrs.cleaned.overlap,gz.output,col.names=T,row.names=F,quote=F,sep="\t")
  write.table(all.chrs.cleaned.overlap$SNP,paste(outdir,"/",GWAS.file,".",TARGET.prefix,".summaRyQC.snplist",sep=""),col.names=T,row.names=F,quote=F)
  
  logme(paste("Creating new bed/bim/fam file for snps overlapping discovery and target datasets",sep=""))
  cat(paste("Creating new bed/bim/fam file for snps overlapping discovery and target datasets\n",sep=""))
  
  ### Create new target bed/bim/fam for overlapping SNPs with discovery 
  system(paste(plinkexe," --bfile ",TARGET.prefix," --extract ",GWAS.file,".",TARGET.prefix,".summaRyQC.snplist --make-bed --out ",GWAS.file,".",TARGET.prefix,".forPRS",sep=""))
  
  logme(paste("Running PRSice2 and deriving PRS in target dataset",sep=""))
  cat(paste("Running PRSice2 and deriving PRS in target dataset\n",sep=""))
  system(paste("Rscript ",prsice2r," --prsice ",prsice2," --base ",GWAS.file,".",TARGET.prefix,".summaRyQC.txt.gz ","--target ",GWAS.file,".",TARGET.prefix,".forPRS ","--bar-levels 0.05,1 --fastscore --no-regress --interval 0.0005 --out ",GWAS.file,".",TARGET.prefix,".PRSICEOUT",sep=""))
  
  logme(paste("PRS have been created in target dataset",sep=""))
  cat(paste("PRS have been created in target dataset\n",sep=""))
}


#dir.create(paste(outdir,"/",GWAS.file,"_summaryqc",sep=""))system(paste("mv ",GWAS.file,".summaRyQC.txt.gz ",GWAS.file,"_summaryqc/",sep=""))
#system(paste("mv *.remove.* ",GWAS.file,"_summaryqc/",sep=""))

#logme(paste("Moving summaRyQC files to new directory: ",GWAS.file,"_summaryqc",sep=""))
#cat(paste("Moving summaRyQC files to new directory: ",GWAS.file,"_summaryqc\n",sep=""))


#if(PRS=="T"){
#dir.create(paste(outdir,"/",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_summaryqc",sep=""))
#system(paste("mv ",GWAS.file,".",TARGET.prefix,".summaRyQC.txt.gz ",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_PRS/",sep=""))
#system(paste("mv *forPRS* ",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_PRS/",sep=""))
#system(paste("mv *PRSICEOUT* ",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_PRS/",sep=""))

#logme(paste("Moving summaRyQC PRS files to new directory: ",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_PRS/",sep=""))
#cat(paste("Moving summaRyQC PRS files to new directory: ",GWAS.file,"_summaryqc/",GWAS.file,".",TARGET.prefix,"_PRS/\n",sep=""))
#}
