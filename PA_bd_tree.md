From sequences to trees with Bd Fluidigm data
================

For this workshop I will show you how to take our Bd Fluidigm data and
turn it into a phylogenetic tree. Fun\!

First, some preparation

1.  I’ve included the program
    [Astral](https://github.com/smirarab/ASTRAL) in the workshop github
    folder. Navigate to the workshop folder and unzip the .zip file.

**only for mac, step 2 won’t work for windows** 2. I’ve also included
the program newick utils in the workshop github folder. Go to the
workshop folder and unzip it. Then navigate to the folder that has the
newick utils in terminal and run the following to make the binaries we
will need to run:

> ./configure make

Ok now let’s make sure we have the libraries we need and go ahead and
install them:

Now let’s load the libraries and check to make sure they are all loaded.

Now let’s read in the sequence data contained in the “seq\_data” files.
Each fasta file represents an individual Bd sample. This folder includes
Bd samples collected from around Pymatuning, in addition to some
reference Bd samples to contextualize our resutls. These fasta files
have one amplicon sequence per locus that has multiple alleles coded
using [IUPAC ambiguity
codes](https://droog.gs.washington.edu/parc/images/iupac.html).

``` r
#get a vector of all the file names in the folder
files_PA <- list.files("./seq_data", full.names = T)
```

Now let’s see how many sequences we have for each primer/sample.

``` r
#make empty objects
length_table_PA <- tibble(sample = "", length=0, file="")

#extract sample names
samplenames_PA <- sapply(files_PA, function(x)unlist(strsplit(x,split=".fasta_primerfix.fasta_trim"))[1])
samplenames_PA <- sapply(samplenames_PA, function(x)unlist(strsplit(x,split="./seq_data/Sample."))[2])

#count number of sequences per sample
for (i in 1: length(files_PA)){
    seq <- readDNAStringSet(files_PA[i], format = "fasta")
    length_table_PA <- add_row(length_table_PA, sample=samplenames_PA[i], length=length(seq), file=files_PA[i])
}

#plot the results
hist(as.numeric(length_table_PA$length[-1]), breaks=50,main="PA", xlab="Number of sequences")
```

![](PA_bd_tree_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

We can see that each sample has a variable amount of missing data. At
this stage in the process we would choose what cutoff we would like to
use for including samples.

``` r
#get list of only samples with at least 90 sequences
files_PA_90 <- filter(length_table_PA, length>90)

filenames_PA_90 <- files_PA_90$file
length(filenames_PA_90)
```

    ## [1] 82

``` r
#filters out 1 bad sample
```

Ok now that we have chosen our samples, let’s get our data organized. We
are going to do that by creating a matrix m that we will then populate
with our sequence data. Here each row is a sample and each column is a
single locus

``` r
#make empty matrix of the right size
m_PA <-matrix(NA, nrow=nrow(files_PA_90), ncol=243)
#rename the columns for the loci
colnames(m_PA) <- seq(1,243)
#rename the rows for the samples
rownames(m_PA) <- files_PA_90$sample
#covert
seq_matrix <- as_tibble(m_PA)

#populate the empty matrix with the ambiguities sequences for PA
#takes a minute to run
for (i in 1:length(filenames_PA_90)){
  seq <- readDNAStringSet(filenames_PA_90[i])
  for (j in 1:length(seq)){
    col_match <- as.numeric(names(seq))
    m_PA[i,col_match[j]] <- as.character(seq[j])
     }
}

#here might be a cool place to view the m_PA matrix to get an idea of the data we are working with.

#in Rstudio
#View(m_PA)
```

Now we have a matrix *m\_PA* with samples as rows and amplicons as
columns. First we can eliminate amplicons with no data. Let’s get the
average length of each sequence to determine which has no data. We can
also find the min max and mean sequence length to identify potential bad
sequences.

Now let’s explore the problem amplicons. These are ones missing too much
data that we should trim out.

``` r
#identify amplicons with more than 1/2 missing data - only considering the samples of interest and not the reference samples (hence the +13)
badamps_PA <- which(as.numeric(n_bases_PA[,7]) > (.5*(length(filenames_PA_90)+13)))

#now trim out these amos from the matrix
m_trim_PA <- m_PA[,-badamps_PA]

#check how many we are left with
dim(m_trim_PA) 
```

    ## [1]  82 176

Now that we have a trimmed dataset let’s just double check we don’t have
any errors in our sequence data. One way to do this is to check the
difference in length between the min and max amplicon.

``` r
#create an empty object to count bases in the new trimmed matrix
n_bases_PA_trim <- matrix(NA, nrow=ncol(m_trim_PA), ncol=7)

#characterize each locus
for (i in 1:ncol(m_trim_PA)){
  n_bases_PA_trim[i,2] <- mean(nchar(m_trim_PA[,i]),na.rm=T)
  n_bases_PA_trim[i,3] <- min(nchar(m_trim_PA[,i]),na.rm=T)
  n_bases_PA_trim[i,4] <- max(nchar(m_trim_PA[,i]),na.rm=T)
  n_bases_PA_trim[i,5] <- median(nchar(m_trim_PA[,i]),na.rm=T)
  n_bases_PA_trim[i,6] <- max(nchar(m_trim_PA[,i]),na.rm=T) - min(nchar(m_trim_PA[,i]),na.rm=T)
  n_bases_PA_trim[i,7] <- sum(is.na(m_trim_PA[,i]))
}

#rename columns and add amplicon labels
n_bases_PA_trim[,1] <- colnames(m_trim_PA)
colnames(n_bases_PA_trim) <- c("amp","mean","min","max","median","diff","sum_missing")

n_bases_t <- as_tibble(n_bases_PA_trim)
plot(n_bases_t$diff, ylab="seq length difference", xlab="amplicon")
```

![](PA_bd_tree_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#get the amplicon name for those that have really big differences between min and max amplicon length so we can double check there aren't bad sequences there
n_bases_t[which(as.numeric(n_bases_t$diff) > 15),1]
```

    ## # A tibble: 3 x 1
    ##   amp  
    ##   <chr>
    ## 1 41   
    ## 2 225  
    ## 3 243

For each sample with data, aligns all loci separately. This is done only
for the gene tree to species tree approach.

Now we want to make a tree for each locus. This allows us to collapse
all the information contained in one locus into the relative
relationships between the samples. This is called the “gene tree to
species tree” method. These trees will then be used to create another
tree using the program ASTRAL-III. To do this I typically use the
program Geneious (as I will demonstrate) to run the command line program
raxml to make a maximum liklihood tree with bootstrap support. Raxml can
also be run locally after installing it. This installation can be a bit
tricky and making the trees can be time consuming so I’ll show you how I
do it then provide the trees. Really you can use any program that makes
trees with bootstrap support.

Here is the command line version of the raxml command that I ran in
Geneious:

> raxmlHPC-SSE3-MAC -s input.phy -n output -m GTRCAT -f a -x 1 -N 100 -p
> 1

You can also run raxml [online](https://raxml-ng.vital-it.ch/#/). Lots
of options\!

For the purposes of this workshop I’ve made all the trees for us and
they are in the folder *PA\_locus\_trees*

Next we need to concatenate all trees into one file to be input into
Astral. You can do this with a simple command on the command line. In
the directory with all the trees (and only the trees) with nodes
collapsed I ran the following:

for mac: \> cat \* \> combined\_trees.newick

for windows: \> copy /b \*.newick combined\_trees.newick

Here is code that will run this command right from R.

**run for mac only**

``` r
system("cat PA_locus_trees/* > PA_trees_combined.newick")
```

**the file *PA\_combined\_trees.newick* is also in the folder
“intermediate\_files\_jic” so you can just copy paste this into your
working directory if you can’t get this concatenate command to work**

Then we can use the newick utils tools to collapse all nodes with \<10
bootstrap support. Make sure you have extracted the newick utils zipped
file and run the compile and make commands (see above). Unfortunately
this program does not work on windows so for those with windows we will
use a file I have already collapsed for us.

**run for mac only**

``` r
system("./newick-utils-1.6/src/nw_ed PA_trees_combined.newick 'i & b<=10' o > PA_trees_collapse10.newick")
```

**on PC copy and paste the *PA\_trees\_collapse10.newick* from the
folder “intermediate\_files\_jic” to the working directory**

Then run astral on your set of collapsed trees. Double check that you
have extracted the zipped Astral file and that the path is correct.

**run for mac only**

``` r
system("java -jar ./Astral/astral.5.7.5.jar -t 2 -i PA_trees_collapse10.newick -o PA_astral2.tre")
```

**run for PC only** **not tested might not work**

``` r
system("java -jar .\Astral\astral.5.7.5.jar -t 2 -i PA_trees_collapse10.newick -o PA_astral2.tre")
```

Now we have a TREE\!\! Let’s look at it\! Here we can use a handy
program called [AstralPlane](https://github.com/chutter/AstralPlane) to
display our tree

``` r
astral.data = createAstralPlane(astral.tree = "PA_astral2.tre",
                                outgroups = "Bd_RFM380",
                                tip.length = 1)


astralProjection(astral.plane = astral.data,
                 local.posterior = TRUE,
                 pie.plot = TRUE,
                 pie.data = "genetree",
                 save.file = "astral_plot.pdf",
                 pie.colors = c("purple", "blue", "green"),
                 node.color.text = c("white"),
                 node.color.bg = c("black"),
                 node.label.size = 0.3,
                 tip.label.size = 0.5,
                 pie.chart.size = 0.75)
```

    ## quartz_off_screen 
    ##                 2

Here we have a lot of low support nodes showing up. Let’s go ahead and
collapse some of the low support nodes and then re-display our tree with
support values. Here we will run Astral one more time to use a different
notation (-t 3) use newick utils again to collapse low support nodes.

**run on mac only**

``` r
system("java -jar ./Astral/astral.5.7.5.jar -t 3 -i PA_trees_collapse10.newick -o PA_astral3.tre")
```

**run on PC only** **not tested, might not work**

``` r
system("java -jar .\Astral\astral.5.7.5.jar -t 3 -i PA_trees_collapse10.newick -o PA_astral3.tre")
```

Now let’s collapse nodes with low support

**run on mac only**

``` r
#collapse nodes with <0.5 posterior support
system("./newick-utils-1.6/src/nw_ed PA_astral3.tre 'i & b<=0.5' o > PA_astral_collapse5.tre")
```

**on PC copy and paste the *PA\_astral\_collapse5.tre* from the folder
“intermediate\_files\_jic” to the working directory**

Read in the new astral tree and plot it.

``` r
astral_collapse <- read.tree("PA_astral_collapse5.tre")
astral_collapse <- root(astral_collapse, "Bd_RFM380")

plot(astral_collapse, show.node.label=T)
```

    ## Warning in plot.phylo(astral_collapse, show.node.label = T): 43 branch length(s)
    ## NA(s): branch lengths ignored in the plot

![](PA_bd_tree_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Now it is probably helpful to read in a metadata file that has more
information about each sample so we can be more specific with our
labels.

``` r
#read in raw metadata file
PA_meta <- read.csv("PA_samples_metadata.csv")

#pull out names from the tree
samples <- as.data.frame(astral_collapse$tip.label)
colnames(samples) <- "Sample_ID"
#now join the info from the metatable to the sample names of the tree
sample_info <- left_join(samples, PA_meta, by = "Sample_ID")

#now make a copy of the tree so we can change the tip labels
astral_collapse_species <- astral_collapse
astral_collapse_species$tip.label <- sample_info$Species_code
refs <- which(astral_collapse_species$tip.label == "")
astral_collapse_species$tip.label[refs] <- sample_info$GENOASSIGN[refs]

#now plot
plot(astral_collapse_species, show.node.label=T, tip.color=as.numeric(as.factor(sample_info$GENOASSIGN)), cex=0.6, main="Astral Tree for PA Bd Samples")
```

    ## Warning in plot.phylo(astral_collapse_species, show.node.label = T, tip.color
    ## = as.numeric(as.factor(sample_info$GENOASSIGN)), : 43 branch length(s) NA(s):
    ## branch lengths ignored in the plot

![](PA_bd_tree_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->
