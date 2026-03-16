args <- commandArgs(trailingOnly = TRUE)
raw_args <- commandArgs(trailingOnly = FALSE)
file_args <- raw_args[grep("^--file=", raw_args)]
script_path <- if (length(file_args) >= 1) sub("^--file=", "", file_args[[1]]) else "."
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

parse_args <- function(args) {
  out <- list(
    pareto_csv = file.path(repo_root, "generated", "figures", "optimization", "lhs200_reference_pareto_blackbox_full", "pareto_candidates.csv"),
    out = file.path(repo_root, "generated", "figures", "optimization", "lhs200_reference_pareto_blackbox_full", "pareto_front.png")
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--pareto-csv") {
      out$pareto_csv <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--out") {
      out$out <- args[[i + 1]]
      i <- i + 2
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
  }
  out
}

cfg <- parse_args(args)
pareto_csv <- normalizePath(cfg$pareto_csv, mustWork = TRUE)
out_path <- normalizePath(cfg$out, mustWork = FALSE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

df <- read.csv(pareto_csv, stringsAsFactors = FALSE)
df <- df[order(df$weight_tumor), , drop = FALSE]
pal <- grDevices::colorRampPalette(c("#1f77b4", "#d62728"))(nrow(df))
days <- c(0, 7, 14, 21, 42, 63)
ref_doses <- c(1.6, 10, 10, 20, 20, 20)

png(filename = out_path, width = 2200, height = 1200, res = 220)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.2), oma = c(0, 0, 2.5, 0))

nd <- tolower(as.character(df$nondominated)) == "true"
x <- df$tumor_loss
y <- df$il6_loss
xr <- range(c(x, 0))
yr <- range(c(y, 0))
plot(0, 0, type = "n", xlim = xr + c(-0.01, 0.01), ylim = yr + c(-0.01, 0.01),
     xlab = "Tumor loss vs reference", ylab = "IL6 loss vs reference",
     main = "Pareto objective space")
abline(h = 0, v = 0, lty = 2, col = "gray60")
points(x[!nd], y[!nd], pch = 21, bg = "gray75", col = "gray35", cex = 1.5)
points(x[nd], y[nd], pch = 21, bg = pal[nd], col = "black", cex = 1.8)
if (sum(nd) >= 2) {
  ord <- order(x[nd], y[nd])
  lines(x[nd][ord], y[nd][ord], lwd = 2.0, col = "#444444")
}
points(0, 0, pch = 4, cex = 1.7, lwd = 2, col = "black")
text(0, 0, labels = "ref", pos = 4, cex = 0.8, xpd = TRUE)
for (i in seq_len(nrow(df))) {
  lab <- sprintf("w=%.2g", df$weight_tumor[i])
  text(x[i], y[i], labels = lab, pos = 4, cex = 0.8, col = if (nd[i]) "black" else "gray35", xpd = TRUE)
}
legend("bottomleft",
       legend = c("Nondominated", "Dominated", "Reference"),
       pch = c(21, 21, 4),
       pt.bg = c("#1f77b4", "gray75", NA),
       pt.cex = c(1.4, 1.4, 1.2),
       col = c("black", "gray35", "black"),
       bty = "n", cex = 0.9)

plot(0, 0, type = "n", xlim = range(days), ylim = c(0, max(c(ref_doses, as.matrix(df[, grep("^dose[0-9]+_mg$", names(df))]))) * 1.08),
     xlab = "Dose day", ylab = "Dose (mg)", main = "Dose schedules for Pareto candidates")
lines(days, ref_doses, lwd = 2, lty = 2, col = "black")
points(days, ref_doses, pch = 16, cex = 0.8, col = "black")
for (i in seq_len(nrow(df))) {
  doses <- as.numeric(df[i, grep("^dose[0-9]+_mg$", names(df))])
  col_i <- if (nd[i]) pal[i] else "gray70"
  lwd_i <- if (nd[i]) 2.0 else 1.2
  lines(days, doses, col = col_i, lwd = lwd_i)
  points(days, doses, pch = 16, cex = 0.75, col = col_i)
}
legend("topright",
       legend = c("Reference 1.6/10/10/20/20/20", "Optimized candidates"),
       lty = c(2, 1), lwd = c(2, 2), col = c("black", "#1f77b4"),
       bty = "n", cex = 0.9)

mtext("Cohort-level Pareto optimization against the 2020 paper reference regimen", outer = TRUE, cex = 1.15)
dev.off()

cat(out_path, "\n")
