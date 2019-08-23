Cleaning Up RNASeq and Sample Attribute Dataframes
================
Michael Kesling
8/20/2019

### Comparing Normalized RNA-Seq Data to Sample/Subject Attribute Data

I am in the process of creating design matrices for performing machine learning on a number of breast samples, both healthy and tumor.

I've assembled a bunch of sample attributes as well as data that have been heavily normalized from the [Wang ComBat batch normalization paper](https://www.nature.com/articles/sdata201861.pdf).

### Importing Wang Dataset RNASeq IDs

We start by importing the 3 dataframes of Wang data, and extract the RNASeq sample IDs and place them into a list

``` r
wangTCGAnormBreast <- read.table("../data/Wang_Norm_Data/brcarsemfpkmtcga.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
colnames(wangTCGAnormBreast) <- gsub("[.]", "-", colnames(wangTCGAnormBreast))
wangTCGAtumorBreast <- read.table("../data/Wang_Norm_Data/brcarsemfpkmtcgat.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
colnames(wangTCGAtumorBreast) <- gsub("[.]", "-", colnames(wangTCGAtumorBreast))
wangGTEXnormBreast <- read.table("../data/Wang_Norm_Data/breastrsemfpkmgtex.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
colnames(wangGTEXnormBreast) <- gsub("[.]", "-", colnames(wangGTEXnormBreast))
RNASeq_IDs <- c(colnames(wangGTEXnormBreast), colnames(wangTCGAnormBreast), colnames(wangTCGAtumorBreast))
```

### Importing Attribute Dataframes

We'll next import the sample attribute data into data frames:

``` r
GTEX_AttributesDF <- read.table("GTEX_Attributes.txt", sep="\t", header=TRUE,
                              stringsAsFactors = FALSE)
TCGA_Attributes_fullDF <- read.table("TCGA_Attributes_full.txt", sep="\t", header=TRUE,
                              stringsAsFactors = FALSE)
print(c(as.character(dim(GTEX_AttributesDF)), as.character(dim(TCGA_Attributes_fullDF))))
```

    ## [1] "15438" "11"    "1097"  "18"

We can see that GTEX has many more instances, but that's because each donor had many organs harvested. TCGA only deals with breast tissue, tumorous or otherwise. Let's see how we'd subset the GTEx Attribute data frame by seeing what tissue types there are.

``` r
gtex_tissue_counts <- tibble(GTEX_AttributesDF$SMTSD) %>% group_by(GTEX_AttributesDF$SMTSD) %>% tally()
head(gtex_tissue_counts, 25)
```

    ## # A tibble: 25 x 2
    ##    `GTEX_AttributesDF$SMTSD`                    n
    ##    <chr>                                    <int>
    ##  1 Adipose - Subcutaneous                     523
    ##  2 Adipose - Visceral (Omentum)               363
    ##  3 Adrenal Gland                              205
    ##  4 Artery - Aorta                             315
    ##  5 Artery - Coronary                          183
    ##  6 Artery - Tibial                            543
    ##  7 Bladder                                     11
    ##  8 Brain - Amygdala                           124
    ##  9 Brain - Anterior cingulate cortex (BA24)   148
    ## 10 Brain - Caudate (basal ganglia)            199
    ## # … with 15 more rows

We see that there are 306 *Breast - Mammary Tissue* tissue samples in all of the GTEX\_AttributesDF, which represent the attributes on all GTEx samples as of August 2019. The Wang publication was submitted in the Fall of 2017, so the work would have been late 2016-early 2017.

### Subsetting based on Tissue = Breast

Next, we'll subset the GTEX\_AttributesDF data frame by those whose detailed tissue is *Breast - Mammary Tissue*, and we see if there are any cases where more than 1 sample is derived from one donor.

``` r
GTEX_Attributes_Breast <- GTEX_AttributesDF %>% filter(SMTSD=="Breast - Mammary Tissue")
tibble(GTEX_Attributes_Breast$SUBJID) %>% group_by(GTEX_Attributes_Breast$SUBJID) %>% tally() %>% arrange(desc(n))
```

    ## # A tibble: 304 x 2
    ##    `GTEX_Attributes_Breast$SUBJID`     n
    ##    <chr>                           <int>
    ##  1 GTEX-145ME                          2
    ##  2 GTEX-ZDTT                           2
    ##  3 GTEX-1117F                          1
    ##  4 GTEX-111YS                          1
    ##  5 GTEX-1122O                          1
    ##  6 GTEX-117XS                          1
    ##  7 GTEX-117YX                          1
    ##  8 GTEX-1192X                          1
    ##  9 GTEX-11DXW                          1
    ## 10 GTEX-11DXY                          1
    ## # … with 294 more rows

We see that 2 subjects, GTEX-145ME and GTEX-ZDTT, have 2 samples each. However, those may not be relevant for the Wang dataset.

### Cross-Referencing RNASeq IDs by GTEx Attribute Dataframe

We need to ensure that all GTEx RNA samples from the Wang dataset have corresponding entried in the GTEx Attribute dataframe.

``` r
fracGTEX <- dim(inner_join(tibble(RNASeq_IDs), GTEX_Attributes_Breast, by=c("RNASeq_IDs" = "SAMPID")))[1] / sum(grepl("^GTEX", colnames(wangGTEXnormBreast)))
print(paste0("The fraction of Wang RNASeq samples having entries in the GTEx attribute dataframe is ", fracGTEX, " ."))
```

    ## [1] "The fraction of Wang RNASeq samples having entries in the GTEx attribute dataframe is 1 ."

That is perfect. Let's check on the 2 GTEx donors with duplicate samples:

``` r
print(paste0("RNASeq replicate number for GTEX-145ME is ", sum(grepl("GTEX-145ME", RNASeq_IDs)),
             " and the RNASeq replicate number for GTEX-ZDTT is ", sum(grepl("GTEX-ZDTT", RNASeq_IDs)), "."))
```

    ## [1] "RNASeq replicate number for GTEX-145ME is 1 and the RNASeq replicate number for GTEX-ZDTT is 0."

So GTEX-145ME is in the Wang dataset, but only one sample of it, and GTEX-ZDTT isn't present at all.

### Modifying TCGA Attributes Full Dataframe

When the TCGA Attributes file was created, it was indexed by something called *subject sample id*, which is entirely distinct from the TCGA IDs that are used here. Further, each patient usually has several TCGA IDs, as it's not unusual for both healthy and tumorous tissue to be harvested. That being said, the *primary\_diagnosis* will refer to the tumor, and not the healthy sample. That means some of the TCGA Attributes will apply to some RNASeq samples and some will not. Towards this end, we'll modify the TCGA Attributes Full Dataframe so that it's indexed by TCGA ID, each with the right attributes.

``` r
createIncrementalDF <- function(df){
   tcgas <- unlist(strsplit(df$tcgaID, ";"))
   tmpDF <- do.call("rbind",replicate(length(tcgas),df, 
                                      simplify=FALSE))
   i = 1
   for(tcga in tcgas){
      if(grepl("-01.$", tcga)==FALSE){
         tmpDF$primary_diagnosis[i] = "healthy"
         tmpDF$tumor_stage[i] = "non-cancer"
      }
      tmpDF$tcgaID[i] = tcga
      i <- i + 1
   }
   return(tmpDF)
}

TCGA_Attributes_expandedDF <- createIncrementalDF(TCGA_Attributes_fullDF[1,]) 
for(j in 2:dim(TCGA_Attributes_fullDF)[1]){
   TCGA_Attributes_expandedDF <- rbind(TCGA_Attributes_expandedDF, 
                                       createIncrementalDF(TCGA_Attributes_fullDF[j,]))
}  

sprintf("Number rows for initial TCGA sample attribute df is %d, and the number of rows after putting a unique TCGA-ID per row is %d", dim(TCGA_Attributes_fullDF)[1], dim(TCGA_Attributes_expandedDF)[1])
```

    ## [1] "Number rows for initial TCGA sample attribute df is 1097, and the number of rows after putting a unique TCGA-ID per row is 3353"

### Verifying that all Wang TCGA RNASeq samples have Attribute Data

We'll split the TCGA Attribute Data by healthy/tumor and then determine for the Wang RNASeq IDs what fraction of healthy/tumor in the TCGA Attribute Dataframes. We also need a shortening of the TCGA-ID field while matching, as the RNASeq samples have more information (longer ID) than the ID associated with a patient.

``` r
TCGA_Attr_HealthyDF <- TCGA_Attributes_expandedDF %>% filter(primary_diagnosis == "healthy")
TCGA_Attr_TumorDF <- TCGA_Attributes_expandedDF %>% filter(primary_diagnosis != "healthy")

RNASeq_SampleIDs <- gsub("(TCGA-[^-]+[-][^-]+[-][^-]+)[-].*", "\\1", RNASeq_IDs)

fracTCGAhealthy <- dim(inner_join(tibble(RNASeq_SampleIDs), TCGA_Attr_HealthyDF, by=c("RNASeq_SampleIDs" = "tcgaID")))[1] / sum(grepl("^TCGA", colnames(wangTCGAnormBreast)))

fracTCGAtumor <- dim(inner_join(tibble(RNASeq_SampleIDs), TCGA_Attr_TumorDF, by=c("RNASeq_SampleIDs" = "tcgaID")))[1] / sum(grepl("^TCGA", colnames(wangTCGAtumorBreast)))
sprintf("The fraction of Wang RNASeq healthy samples having entries in the TCGA attribute dataframe is %f and the fraction of Wang tumor samples having entries is %f", fracTCGAhealthy, fracTCGAtumor)
```

    ## [1] "The fraction of Wang RNASeq healthy samples having entries in the TCGA attribute dataframe is 1.000000 and the fraction of Wang tumor samples having entries is 0.997963"

So it appears that only a single tumor RNASeq sample is missing attribute information. We'll identify it and delete it from the data frame.

``` r
RNAseqIDtumor <- gsub("(TCGA-[^-]+[-][^-]+[-][^-]+)[-].*", "\\1",colnames(wangTCGAtumorBreast))
full_join(tibble(RNAseqIDtumor), TCGA_Attr_TumorDF, by=c("RNAseqIDtumor" = "tcgaID")) %>% filter(is.na(vital_status))
```

    ## # A tibble: 4 x 18
    ##   RNAseqIDtumor sampleId gender race  days_to_birth ethnicity vital_status
    ##   <chr>         <chr>    <chr>  <chr>         <int> <chr>     <chr>       
    ## 1 Hugo_Symbol   <NA>     <NA>   <NA>             NA <NA>      <NA>        
    ## 2 Entrez_Gene_… <NA>     <NA>   <NA>             NA <NA>      <NA>        
    ## 3 TCGA-BH-A0B2… <NA>     <NA>   <NA>             NA <NA>      <NA>        
    ## 4 TCGA-BH-A0B2… <NA>     <NA>   <NA>             NA <NA>      <NA>        
    ## # … with 11 more variables: age_at_index <int>, days_to_death <int>,
    ## #   year_of_birth <int>, year_of_death <int>, primary_diagnosis <chr>,
    ## #   tumor_stage <chr>, age_at_diagnosis <int>, morphology <chr>,
    ## #   prior_treatment <chr>, tissue_or_organ_of_origin <chr>,
    ## #   prior_malignancy <chr>

And it appears that there are 2 *TCGA-BH-A0B2-01A*-derived RNA Samples that have no corresponding sample attribute data. We will therefore remove them from the dataframe containing the RNASeq data, *wangTCGAtumorBreast*.

``` r
dimBefore <- dim(wangTCGAtumorBreast)
toRemove <- c(grep("TCGA-BH-A0B2-01A", colnames(wangTCGAtumorBreast)))
wangTCGAtumorBreast <- wangTCGAtumorBreast[,-toRemove]
grep("TCGA-BH-A0B2-01A", colnames(wangTCGAtumorBreast))
```

    ## integer(0)

``` r
dimAfter <- dim(wangTCGAtumorBreast)
sprintf("We can see that the *TCGA-BH-A0B2-01A*-derived RNA samples were 2 in number, as the number of columns before removal was %d, and after was %d", dimBefore[2], dimAfter[2])
```

    ## [1] "We can see that the *TCGA-BH-A0B2-01A*-derived RNA samples were 2 in number, as the number of columns before removal was 984, and after was 982"

### Selecting Samples for Analysis

As it stands, there are 980 tumor breast samples and a total of 199 healthy breast samples having sample attribute data. The end goal is to perform machine learning algorithms of these data. However, if we work on the entire dataset, then 982 / (982+199) = 83% of all samples are tumors, and there's only a potential 17% classification improvement to be made over the baseline model (selecting all samples are tumors). In order that we can better discriminate between different machine learning outputs, we'll set 50% as tumors and 50% as healthy. So we should randomly select 199 tumorous samples to add to our healthy sets.

``` r
set.seed(3235321)                                  # random seed
idx <- sample(3:dim(wangTCGAtumorBreast)[2], 199, replace=FALSE)
idx <- c(1,2,idx)                                  # first 2 cols are gene info
wangTCGAtumorBreastsubset <- wangTCGAtumorBreast[,idx]
dim(wangTCGAtumorBreastsubset)
```

    ## [1] 19738   201

We will also want to ensure that the order of the genes is the same in all 3 dataframes.

``` r
all(wangGTEXnormBreast$Entrez_Gene_Id == wangTCGAnormBreast$Entrez_Gene_Id) && 
all(wangGTEXnormBreast$Entrez_Gene_Id == wangTCGAtumorBreast$Entrez_Gene_Id)
```

    ## [1] TRUE

### Creating a Single RNASeq (FPKM) Dataframe and Saving to a File

As everything is now in order, we'll bind together the RNASeq data frames to create one with 398 batch-normalized RNASeq measurements. As mentioned in a different document, the units of measurement are *FPKM* derived from the *RSEM* algorithm that was quantile normalized.

``` r
wangBreastFPKM398 <- cbind(wangGTEXnormBreast, wangTCGAnormBreast[3:dim(wangTCGAnormBreast)[2]],
      wangTCGAtumorBreastsubset[3:dim(wangTCGAtumorBreastsubset)[2]])
dim(wangBreastFPKM398)
```

    ## [1] 19738   400

The *wangBreastFPKM398* dataframe seems to be the correct size, including the 2 gene information columns. We write it to a file.

``` r
write.table(wangBreastFPKM398, "wangBreastFPKM398.txt", sep="\t", row.names = FALSE)
```

### Creating the Design Matrix

Some of the sample attributes can be used the create a design matrix X. \#\#\# Exploring Predictor Variables and Cross-Mapping IDs

``` r
paste("TCGA: ", colnames(TCGA_Attr_TumorDF), " GTEx: ", colnames(GTEX_AttributesDF))
```

    ##  [1] "TCGA:  sampleId  GTEx:  SAMPID"                   
    ##  [2] "TCGA:  gender  GTEx:  SMRIN"                      
    ##  [3] "TCGA:  race  GTEx:  SMTS"                         
    ##  [4] "TCGA:  days_to_birth  GTEx:  SMTSD"               
    ##  [5] "TCGA:  ethnicity  GTEx:  SME2MPRT"                
    ##  [6] "TCGA:  vital_status  GTEx:  SMMAPRT"              
    ##  [7] "TCGA:  age_at_index  GTEx:  SMGNSDTC"             
    ##  [8] "TCGA:  days_to_death  GTEx:  SUBJID"              
    ##  [9] "TCGA:  year_of_birth  GTEx:  SEX"                 
    ## [10] "TCGA:  year_of_death  GTEx:  AGE"                 
    ## [11] "TCGA:  primary_diagnosis  GTEx:  DTHHRDY"         
    ## [12] "TCGA:  tumor_stage  GTEx:  SAMPID"                
    ## [13] "TCGA:  age_at_diagnosis  GTEx:  SMRIN"            
    ## [14] "TCGA:  morphology  GTEx:  SMTS"                   
    ## [15] "TCGA:  prior_treatment  GTEx:  SMTSD"             
    ## [16] "TCGA:  tissue_or_organ_of_origin  GTEx:  SME2MPRT"
    ## [17] "TCGA:  prior_malignancy  GTEx:  SMMAPRT"          
    ## [18] "TCGA:  tcgaID  GTEx:  SMGNSDTC"

GTEx only has about half the sample attribute fields, and most of them have non-informative labels. GTEx label Description (and whether used for design matrix) SAMPID GTEx public sample ID --&gt; Fully maps to RNASeq file IDs SMRIN RIN Number. Not used in design matrix SMTS Tissue Type. Not used. See SMTSD. SMTSD Tissue Type with more detail. Maps to TCGA: tissue\_or\_organ\_of\_origin, but we'll need to standardize the language for creating Factors. SME2MPRT 2-end mapping rate: (\# 2-end mapped reads / \# 1-end mapped reads). Not used. SMMAPRT Mapping rate: Ratio of total mapped reads to total reads. Not used. SMGNSDTC Genes detected. Total \# genes with 5+ exon-mapped reads. Not used SUBJID Truncated SAMPID. Not used. SEX 1=Male, 2=Female. Used, and maps to TCGA: gender. Need confirm same values. AGE Age at time of organ donation (death). Only gives age by decade. Will need to alter TCGA: age DTHHRDY Hardy score describing death. Not used because only GTEx has this info.

There are *2 design matrices* we'll create. One with all the sample attributes we'd *like* to use, and one where the attributes will be populated across the dataset.

#### Ideal Sample Attributes

The ideal attributes, amongst which we have, are: \* Age (age\_at\_index)
\* Ethnicity
\* Race
\* Vital Status
\* primary diagnosis
\* Tumor Stage \* Tissue or Organ of Origin \* Prior Malignancy
I am skipping *Sex* because only 1% of the breast samples are of men (see below). I'm simply going to remove them from the dataset.

With *Age*, it's only given as a range (e.g. 30-39), so rather than lose the granularity of the TCGA datasets, I'll take a range and place it at the midpoint (e.g. 35).

In the case of the GTEx, which has almost no sample info publicly available, I'll leave the missing data as unknown *NA*, which can either be removed or imputed later.

So we need to perform mapping between the TCGA and GTEx sample attribute columns and we need to standardize mapping to the RNASeq file.

### Subsetting Sample Attribute files

Before looking at the ratios of predictor values, we'll subset each of the attribute file by whether there is an RNASeq sample in hand.

``` r
GTEX_RNASeqIDs <- tibble(RNASeq_IDs) %>% filter(grepl("^GTEX", RNASeq_IDs))
GTEX_Attr_Breast_RNA  <- inner_join(GTEX_RNASeqIDs, GTEX_Attributes_Breast, by=c("RNASeq_IDs" = "SAMPID"))


TCGAhealthy_RNASeqIDs <- as.data.frame(cbind(RNASeq_IDs, RNASeq_SampleIDs)) %>%
   filter(grepl("^TCGA", RNASeq_IDs)) %>% filter(grepl("-11[AB]-", RNASeq_IDs))
# now need to cbind abbreviated ID and join with that new id.
TCGAhealthy_Attr_RNA <- inner_join(TCGAhealthy_RNASeqIDs, TCGA_Attr_HealthyDF, by=c("RNASeq_SampleIDs"="tcgaID"))
```

    ## Warning: Column `RNASeq_SampleIDs`/`tcgaID` joining factor and character
    ## vector, coercing into character vector

``` r
# with the TCGA tumors, we need to first grab the subsetted RNA_SeqIDs from
# the wangTCGAtumorBreastsubset file and from that, we'll generate the RNA_SampleSeqIDs
TCGAtumorSelected_RNASeqIDs <-tibble(RNASeqTumor=colnames(wangTCGAtumorBreastsubset)) %>%
   filter(grepl("^TCGA", RNASeqTumor))
TCGAtumorSelected_RNASeqIDs <- cbind(TCGAtumorSelected_RNASeqIDs, RNASeqTumorSample= gsub("(TCGA-[^-]+[-][^-]+[-][^-]+)[-].*", "\\1", TCGAtumorSelected_RNASeqIDs$RNASeqTumor))
# now we perform the inner join to gather the relevant tumor sample attributes
TCGAtumorSelected_Attr_RNA <- inner_join(TCGAtumorSelected_RNASeqIDs, TCGA_Attr_TumorDF, by=c("RNASeqTumorSample"="tcgaID"))
```

    ## Warning: Column `RNASeqTumorSample`/`tcgaID` joining factor and character
    ## vector, coercing into character vector

Let's start with tissue\_of\_origin. The TCGA has 6 values:

``` r
paste(levels(factor(TCGAhealthy_Attr_RNA$tissue_or_organ_of_origin)),
levels(factor(TCGAtumorSelected_Attr_RNA$tissue_or_organ_of_origin)))
```

    ## [1] "Breast, NOS Breast, NOS"                                      
    ## [2] "Lower-inner quadrant of breast Upper-inner quadrant of breast"
    ## [3] "Overlapping lesion of breast Upper-outer quadrant of breast"  
    ## [4] "Upper-outer quadrant of breast Breast, NOS"

"Breast, NOS" means invasive ductal carcinoma, which is from the mammary tissue and not the surrounding adipose tissue. It should be noted that the healthy TCGA samples are all derived from patients that also had cancer, but at this point, we don't know from what part of the breast the healthy breast samples from TCGA came from, though we could dig farther into the GDC later to find out. In that case, we'd be able to use *breast subsection* as a potential confounder / predictor variable. So for now, I've decided to leave this out.

and GTEx has 1 level:

``` r
levels(factor(GTEX_Attr_Breast_RNA$SMTSD))
```

    ## [1] "Breast - Mammary Tissue"

So since all samples are from the breast, it won't be included as a variable in the design matrix.

In GTEX, the sex variable is called 'sex' and has values 1=Male, 2=Female In TCGA, the sex variable is called 'gender' and has values 'female' and 'male'.

``` r
TCGAhealthy_Attr_RNA %>% group_by(gender) %>% summarize(Count = n())
```

    ## # A tibble: 2 x 2
    ##   gender Count
    ##   <chr>  <int>
    ## 1 female   109
    ## 2 male       1

For TCGA, the ratio of females to males in healthy tissue is about 100:1

``` r
TCGAtumorSelected_Attr_RNA %>% group_by(gender) %>% summarize(Count = n())
```

    ## # A tibble: 2 x 2
    ##   gender Count
    ##   <chr>  <int>
    ## 1 female   197
    ## 2 male       2

And for TCGA cancer samples the female to male ratio is still about 100:1

``` r
GTEX_Attr_Breast_RNA %>% group_by(SEX) %>% summarise(Count = n())
```

    ## # A tibble: 1 x 2
    ##     SEX Count
    ##   <int> <int>
    ## 1     2    89

GTEx has all samples from females.

#### Modifying GTEX Attribute Dataframe

So in the analysis, I'll probably filter out the *males*. For now, I'll convert the GTEX *SEX* category using the *female* value rather than *2*. I'm also going to change the age range to a midpoint age, and fill in other missing fields.

``` r
GTEX_Attr_Breast_RNA$SEX <- rep("female", dim(GTEX_Attr_Breast_RNA)[1])
names(GTEX_Attr_Breast_RNA)[names(GTEX_Attr_Breast_RNA) == "SEX"] <- "gender"
GTEX_Attr_Breast_RNA$AGE <- as.numeric(gsub("-[0-9]{2}", "", GTEX_Attr_Breast_RNA$AGE)) + 5
names(GTEX_Attr_Breast_RNA)[names(GTEX_Attr_Breast_RNA) == "AGE"] <- "age"
names(GTEX_Attr_Breast_RNA)[names(GTEX_Attr_Breast_RNA) == "RNASeq_IDs"] <- "RNASeq_ID"
GTEX_Attr_Breast_RNA$primary_diagnosis <- "healthy"
GTEX_Attr_Breast_RNA$tumor_stage <- -1  # numeric variable whose value means something
GTEX_Attr_Breast_RNA$ethnicity <- NA
GTEX_Attr_Breast_RNA$race <- NA
GTEX_Attr_Breast_RNA$prior_malignancy <- "no"
GTEX_Attr_Breast_RNA$vital_status <- "Dead"  # death unrelated to cancer
```

#### Modifying and Cleanup of the TCGA Attribute Dataframes

I'm converting the *tumor\_stage* variable into a numeric variable, as the stage of cancer progresses along with the stages -1 -&gt; 0 -&gt; 1 -&gt; 2 -&gt; 3 -&gt; 4. So here I will do that to both dataframes. Grade GX, also *stage x* means that it cannot be determined, so I'll replace those with NA. During this replacement, I found all the Tumor RNA samples to be either Stage 1, 2, 3, or unknown. No Stage 0 or Stage 4 were present.

``` r
TCGAtumorSelected_Attr_RNA$tumor_stage <- gsub("stage iii.*", 3, TCGAtumorSelected_Attr_RNA$tumor_stage) %>% 
   {gsub("stage ii.*", 2, .)} %>% {gsub("stage i.*", 1, .)} %>% 
   {gsub("not reported", NA, .)} %>% {gsub("stage x", NA, .)}
# I will set healthy breast tissue to stage -1:
TCGAhealthy_Attr_RNA$tumor_stage <- -1
# fix a label
names(TCGAhealthy_Attr_RNA)[names(TCGAhealthy_Attr_RNA) == "age_at_index"] <- "age"
names(TCGAtumorSelected_Attr_RNA)[names(TCGAtumorSelected_Attr_RNA) == "age_at_index"] <- "age"
TCGAhealthy_Attr_RNA$ethnicity <- gsub("not reported", NA, TCGAhealthy_Attr_RNA$ethnicity)
TCGAhealthy_Attr_RNA$race <- gsub("not reported", NA, TCGAhealthy_Attr_RNA$race)
TCGAtumorSelected_Attr_RNA$ethnicity <- gsub("not reported", NA, TCGAtumorSelected_Attr_RNA$ethnicity)
TCGAtumorSelected_Attr_RNA$race <- gsub("not reported", NA, TCGAtumorSelected_Attr_RNA$race)
names(TCGAtumorSelected_Attr_RNA)[names(TCGAtumorSelected_Attr_RNA) == "RNASeqTumor"] <- "RNASeq_ID"
names(TCGAhealthy_Attr_RNA)[names(TCGAhealthy_Attr_RNA) == "RNASeq_IDs"] <- "RNASeq_ID"
TCGAtumorSelected_Attr_RNA$prior_malignancy <- gsub("not reported", NA, TCGAtumorSelected_Attr_RNA$prior_malignancy)
```

### Create Design Matrix of Sample Attributes based on RNASeq File Headers

The ordering of samples in the RNASeq dataset is GTEX, then healthy TCGA, then tumorous TCGA. I'm going to merge the relevant sample attribute data into a single dataframe in that order, as well as save it to a file.

``` r
cols <- c("RNASeq_ID", "gender", "age", "race", "ethnicity", "prior_malignancy", "vital_status",
          "primary_diagnosis", "tumor_stage")
sampleAttr398 <- rbind(GTEX_Attr_Breast_RNA[cols], TCGAhealthy_Attr_RNA[cols],
TCGAtumorSelected_Attr_RNA[cols])
write.table(sampleAttr398, "sampleAttr398BRCA.txt", sep="\t", row.names=FALSE)
```

### Final Data Matrix

I'll finish off by creating a data matrix which combines the RNASeq data with the Sample Attribute data.

``` r
tmp <- cbind(colnames(sampleAttr398), rep(NA, 9), t(sampleAttr398))
colnames(tmp) <- tmp[1,]
colnames(tmp)[colnames(tmp)=="RNASeq_ID"] <- "Hugo_Symbol"
colnames(tmp)[2] <- "Entrez_Gene_Id"
wangBreastFPKM389_Attrib <- rbind(wangBreastFPKM398, tmp[2:9,])
write.table(wangBreastFPKM389_Attrib, "wangBreastFPKM389_Attrib.txt", sep="\t", row.names=FALSE)
```

### Summary

This markdown file created 3 output files:
1. wangBreastFPKM398.txt which was the RNASeq (FPKM) batch-normalized data which were selected for further analysis. The 398 samples include 199 cancerous and 199 healthy, all of which had some sample attribute data associated with them.
2. sampleAttr398BRCA.txt which contained 8 different sample attributes for each of the 398 samples.
3. wangBreastFPKM389\_Attrib.txt which contained the combined data from files 1 and 2.
