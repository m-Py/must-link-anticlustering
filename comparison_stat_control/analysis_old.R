

## old analysis stuff for only 1 run
pvalues <- results[, grepl("p_", colnames(results))]
cors <- results[, grepl("cor_", colnames(results))]
biases <- results[, grepl("balance_", colnames(results))]

# overall power per method:
sort(colMeans(pvalues < .05)) |> round(2)
t(apply(cors, 2, range)) # this one is more relevant!

# confound can also lower power if analysis controls for batch effects


## If we have a bunch of covariates that are not strongly related to the outcome, 
# stat. control is not necessarily good (i.e., controlling for the covariates)
## (e.g., in simulation set covariate effect = 0 and batch effect = 0; use 
# several covariates of them). Then anticlustering is still useful to prevent
# a reduction in power.

cat("For adjusted analysis: p value of anticlustering was lower (=better) in ", 
    mean(results$p_anticlust_control < results$p_rnd_control) * 100,
    "% of all cases,", 
    format_p(prop.test(sum(results$p_anticlust_control < results$p_rnd_control), nrow(results), p = .50)$p.value),
    "\n"
)

cat("For unadjusted analysis: p value of anticlustering was lower (=better) in ", 
    mean(results$p_anticlust_no_control < results$p_rnd_no_control) * 100,
    "% of all cases,", 
    format_p(prop.test(sum(results$p_anticlust_no_control < results$p_rnd_no_control), nrow(results), p = .50)$p.value),
    "\n"
)

# if there is no batch effect, we have the same p values for anticlustering and random assignment
# because the treatment variable is the same.

## regression to predict power (works better for fewer batches)

## maybe do implement confounding condition

plot_stuff <- function(cor_rnd, cor_anticlust, cor_confound, p_rnd, p_anticlust, p_confounded,
                       xlab = "correlation batch effect / treatment effect", ylab = "p value") {
  xlim <- range(cor_confound, cor_anticlust, cor_confound)
  plot(cor_rnd, p_rnd, xlim = xlim, pch = 2, cex = .5, col = "#ABCDEF",
       xlab = xlab, ylab = ylab)
  points(cor_anticlust, p_anticlust, col = "#333333", cex = .6, pch = 4)
  points(cor_confound, p_confounded, col = "#f0f", cex = .6, pch = 4)
  legend("topright", legend = c("Random Assignment", "Anticlustering", "Confounded Assignment"), 
         pch = c(2, 4), col = c("#ABCDEF", "#333333", "#f0f"), cex = .7)
}


plot_stuff(
  cors$cor_rnd, cors$cor_anticlust, cors$cor_confound,
  pvalues$p_rnd_no_control, pvalues$p_anticlust_no_control, pvalues$p_confound_no_control
)

plot_stuff(
  cors$cor_rnd, cors$cor_anticlust, cors$cor_confound,
  pvalues$p_rnd_control, pvalues$p_anticlust_control, pvalues$p_confound_control
)
## using bias values
plot_stuff(
  biases$balance_rnd, biases$balance_anticlust, biases$balance_confound,
  pvalues$p_rnd_no_control, pvalues$p_anticlust_no_control, pvalues$p_confound_no_control,
  xlab = "Batch balance"
)

plot_stuff(
  biases$balance_rnd, biases$balance_anticlust, biases$balance_confound,
  pvalues$p_rnd_control, pvalues$p_anticlust_control, pvalues$p_confound_control,
  xlab = "Batch balance"
)

significant <- pvalues$p_rnd_no_control < .05
summary(glm(significant ~ cors$cor_rnd, family = binomial()))

cor.test(pvalues$p_rnd_no_control, cors$cor_rnd)
cor.test(cors$cor_anticlust, pvalues$p_anticlust_no_control)
