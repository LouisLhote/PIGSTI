#!/usr/bin/env Rscript
# Chromosome residual-method sexing (parameterized).
#
# Usage:
#   Rscript sexing_residual_method.R input.idx output.pdf [output.tsv]
#
# The .idx file is built by run_pigsti_sexing.py from samtools idxstats:
#   - rows 1..N  = autosomes (chr column "1", "2", ... in numeric order)
#   - row N+1    = X chromosome (chr column "X")
#   - column 1: chr, column 2: length (bp), column 3+: read counts per sample
#
# Autosome count N is inferred from nrow(d) - 1 (no hard-coded karyotype).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript sexing_residual_method.R input.idx output.pdf [output.tsv]")
}

d <- read.table(args[1], header = TRUE)
pdf(args[2])

n_auto <- nrow(d) - 1L
x_row <- nrow(d)
if (n_auto < 3L) {
  stop("Need at least 3 autosomes plus X in .idx file; got nrow=", nrow(d))
}

# Regression residual df ≈ n_autosomes - 2 (matches legacy bos_capra: 29 autosomes, df=27)
df_t <- max(1L, n_auto - 2L)
auto_idx <- seq_len(n_auto)

probs <- numeric(max(0L, ncol(d) - 2L))
ratios <- probs
calls <- character(length(probs))

for (i in 3:ncol(d)) {
  y_reads <- d[x_row, i]
  x_len <- d$length[x_row]

  plot(
    d[auto_idx, i] ~ d$length[auto_idx],
    xlab = "Chromosome length (Mb)",
    ylab = "Number of reads",
    main = colnames(d)[i],
    pch = 16,
    col = "slategrey",
    xaxt = "n"
  )
  points(x_len, y_reads, col = "red", pch = 16)
  text(x_len, y_reads, labels = d$chr[x_row], col = "red", pos = 1)
  for (j in auto_idx) {
    text(d$length[j], d[j, i], labels = d$chr[j], col = "slategrey", pos = 1, cex = 0.5)
  }
  points(x_len, 2 * y_reads, col = "lightgrey", pch = 16)
  text(x_len, 2 * y_reads, labels = "2X", col = "lightgrey", pos = 3)
  axis(side = 1, at = axTicks(1), labels = axTicks(1) / 1e6)

  model <- lm(d[auto_idx, i] ~ d$length[auto_idx])
  abline(model)

  prd <- predict(model, data.frame(d$length[auto_idx]), interval = "confidence", level = 0.95)
  confidence <- data.frame(cbind(d$length[auto_idx], prd[, 2], prd[, 3]))
  confordered <- confidence[order(confidence[, 1]), ]
  lines(confordered$X1, confordered$X2, lty = 3)
  lines(confordered$X1, confordered$X3, lty = 3)

  modelwithX <- lm(d[c(auto_idx, x_row), i] ~ d$length[c(auto_idx, x_row)])
  studreswithX <- rstudent(modelwithX)
  femaleprob <- pt(studreswithX[x_row], df = df_t)

  modelwithtwiceX <- lm(
    c(d[auto_idx, i], d[x_row, i] * 2) ~ d$length[c(auto_idx, x_row)]
  )
  studreswithtwiceX <- rstudent(modelwithtwiceX)
  maleprob <- pt(studreswithtwiceX[x_row], df = df_t, lower.tail = FALSE)

  ratio <- femaleprob / maleprob
  probs[i - 2] <- femaleprob

  if (ratio >= 1) {
    lrt <- -2 * log(maleprob / femaleprob)
    lrtp <- 1 - pchisq(lrt, df = 1)
    legend(
      "topleft",
      c(
        paste0("Female (", signif(ratio, digits = 3), "X more likely)"),
        paste0("p = ", signif(lrtp, digits = 3)),
        paste0("autosomes=", n_auto, ", X=", d$chr[x_row])
      ),
      cex = 0.6
    )
    ratios[i - 2] <- ratio
    calls[i - 2] <- "Female"
  } else {
    lrt <- -2 * log(femaleprob / maleprob)
    lrtp <- 1 - pchisq(lrt, df = 1)
    legend(
      "topleft",
      c(
        paste0("Male (", signif(1 / ratio, digits = 3), "X more likely)"),
        paste0("p = ", signif(lrtp, digits = 3)),
        paste0("autosomes=", n_auto, ", X=", d$chr[x_row])
      ),
      cex = 0.6
    )
    ratios[i - 2] <- 1 / ratio
    calls[i - 2] <- "Male"
  }
}

dev.off()

if (length(args) >= 3) {
  samples <- colnames(d)[3:ncol(d)]
  outdf <- data.frame(
    sample = samples,
    call = calls,
    female_prob = probs,
    likelihood_ratio = ratios,
    n_autosomes = n_auto,
    x_chromosome = d$chr[x_row],
    stringsAsFactors = FALSE
  )
  write.table(outdf, args[3], sep = "\t", row.names = FALSE, quote = FALSE)
}
