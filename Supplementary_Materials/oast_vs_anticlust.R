# oast_vs_anticlust.R

library(anticlust) # I used 0.8.6, but 0.8.5 (currently on CRAN) has the same implementation of fast_anticlustering() which is used below
library(OSAT) # I got version OSAT_1.52.0 from Bioconductor

sessionInfo()

# Reproduce initial example in OSAT vignette
inPath <- system.file('extdata', package='OSAT')
pheno <- read.table(
  file.path(inPath, 'samples.txt'),
  header=TRUE, sep="\t", colClasses="factor"
)


gs <- setup.sample(pheno, optimal=c("SampleType", "Race", "AgeGrp"))
gc <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates') # unfortunately I do not really understand these containers

# Perform OSAT assignment. nSim=5000 are commonly used according to the vignette
start <- Sys.time()
gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=5000)
Sys.time() - start 


# Merge OSAT groupings to initial data frame
df <- merge(pheno, gSetup@data$link[, c("OrigRow", "plates")], by.x = "ID", by.y = "OrigRow")
df <- df[order(df$ID), ]

# This reproduces the p values above!
chisq.test(table(df$SampleType, df$plates))$p.value
chisq.test(table(df$Race, df$plates))$p.value
chisq.test(table(df$AgeGrp, df$plates))$p.value

# Now do anticlustering

categories <- df[, c("SampleType", "Race", "AgeGrp")]
binary_categories <- categories_to_binary(categories)

# Perform k-means anticlustering using default exchange method
start <- Sys.time()
ac <- anticlustering(
  binary_categories,
  K = 6
)
Sys.time() - start 

# This is the k-means anticlustering objective function, which represents discrepancy
# in proportion of each category in each group (higher value = better)
variance_objective(binary_categories, ac)
variance_objective(binary_categories, df$plates)


# I get somewhat slightly different p values than the authors, but they are 
# pretty much the same. They are only numerically different, not regarding the
# order of magnitude; in both cases, the anticlust p values are higher. In the 
# paper, I include the p values from the PDF in their package.

# OSAT:
chisq.test(table(df$SampleType, df$plates))$p.value
chisq.test(table(df$Race, df$plates))$p.value
chisq.test(table(df$AgeGrp, df$plates))$p.value

# ANTICLUST
chisq.test(table(df$SampleType, ac))$p.value
chisq.test(table(df$Race, ac))$p.value
chisq.test(table(df$AgeGrp, ac))$p.value
