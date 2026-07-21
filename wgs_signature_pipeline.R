
# EDIT ONLY THESE PATHS
CFG <- list(
  isomut_snv    = "input/all_SNVs.isomut",
  blacklist_bed = "reference/hg38-blacklist.v2.bed.gz",
  clone_key_csv = "Table_S1_clone_key.csv",
  ref_genome    = "BSgenome.Hsapiens.UCSC.hg38",
  outdir        = "results"
)

# five references shown in Fig 2C
prespec_ref <- list(
  "Control"             = c("control"),
  "SSR (1.25 J)"        = c("SSR", "1.25"),
  "Cisplatin (12.5 µM)" = c("isplatin", "12.5"),
  "BPDE (0.125 µM)"     = c("BPDE", "0.125"),
  "ENU (400 µM)"        = c("ENU", "400")
)

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Biostrings)
  library(GenomeInfoDb); library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38); library(ggplot2)
  library(ggdendro); library(RColorBrewer)
})

dir.create(CFG$outdir, showWarnings = FALSE)
std_chr <- paste0("chr", c(1:22, "X", "Y"))

# Table S1 is the authoritative raw-ID -> protocol-name mapping
clone_key <- read.csv(CFG$clone_key_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_key_cols <- c("Protocol_name", "Original_clone_ID")
stopifnot(all(required_key_cols %in% colnames(clone_key)))
clone_key <- clone_key[!grepl("^\\(?omitted\\)?$", clone_key$Protocol_name, ignore.case = TRUE), ]
stopifnot(nrow(clone_key) == 42L, !anyDuplicated(clone_key$Protocol_name),
          !anyDuplicated(clone_key$Original_clone_ID))
key <- setNames(clone_key$Protocol_name, clone_key$Original_clone_ID)

# read the IsoMut SNV table and validate the complete joint-calling panel
isomut_cols <- c("sample_name","chr","pos","type","score","ref","mut","cov","mut_freq","cleanliness")
snv <- read.delim(CFG$isomut_snv, header = FALSE, comment.char = "#",
                  col.names = isomut_cols, stringsAsFactors = FALSE)
snv$raw <- sub("\\.final\\.bam$", "", snv$sample_name)
raw_samples <- sort(unique(snv$raw))
if (!all(grepl("^chr", unique(snv$chr)))) {
  stop("IsoMut chromosome names are not consistently UCSC-style (expected chr-prefixed contigs).")
}

if ("BAP2KO3" %in% raw_samples) {
  stop("BAP2KO3 is present in all_SNVs.isomut, but it should have been omitted from the IsoMut input panel before joint calling.")
}
unexpected <- setdiff(raw_samples, names(key))
missing <- setdiff(names(key), raw_samples)
if (length(unexpected) > 0L || length(missing) > 0L) {
  stop(paste0(
    "IsoMut sample panel mismatch.",
    if (length(unexpected) > 0L) paste0(" Unexpected: ", paste(unexpected, collapse = ", "), ".") else "",
    if (length(missing) > 0L) paste0(" Missing: ", paste(missing, collapse = ", "), ".") else ""
  ))
}

snv$sample <- unname(key[snv$raw])
stopifnot(all(nchar(snv$ref) == 1L), all(nchar(snv$mut) == 1L),
          all(snv$ref %in% c("A","C","G","T")), all(snv$mut %in% c("A","C","G","T")))

# GRanges per clone on the standard chromosomes only
make_gr <- function(d) {
  gr <- GRanges(seqnames = d$chr, ranges = IRanges(start = d$pos, width = 1L),
                REF = DNAStringSet(d$ref), ALT = DNAStringSetList(as.list(d$mut)))
  genome(gr) <- "hg38"; gr
}
grl_raw <- GRangesList(lapply(split(snv, snv$sample), make_gr))
grl_raw <- GRangesList(lapply(grl_raw, function(gr) {
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), std_chr), pruning.mode = "coarse")
  seqlevels(gr) <- std_chr
  seqinfo(gr)   <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)[std_chr]
  sort(gr)
}))

# drop ENCODE blacklist positions
blacklist <- rtracklayer::import(CFG$blacklist_bed)
seqlevelsStyle(blacklist) <- "UCSC"
blacklist <- keepSeqlevels(blacklist, intersect(seqlevels(blacklist), std_chr), pruning.mode = "coarse")
grl <- GRangesList(lapply(grl_raw, function(gr) subsetByOverlaps(gr, blacklist, invert = TRUE)))
stopifnot(length(grl) == 42L)

group_of <- function(s) sub("[0-9]+$", "", s)
burden <- data.frame(sample = names(grl), n_filtered = lengths(grl),
                     group = factor(group_of(names(grl)), levels = c("NT","UV","CP","BAP","ENU")),
                     row.names = NULL)
expected_group_counts <- c(NT = 9L, UV = 9L, CP = 7L, BAP = 8L, ENU = 9L)
observed_group_counts <- table(burden$group)
stopifnot(identical(as.integer(observed_group_counts), as.integer(expected_group_counts)),
          identical(names(observed_group_counts), names(expected_group_counts)))
message(sprintf("Analysis clones: %d", length(grl)))
cp_burden <- burden$n_filtered[burden$group == "CP"]
fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
message(sprintf(">>> CP mutation burden (median [range]): %s [%s - %s]",
                fmt_int(median(cp_burden)), fmt_int(min(cp_burden)), fmt_int(max(cp_burden))))

# supplementary table, SNV count per clone before and after the blacklist filter
rev_key <- setNames(clone_key$Original_clone_ID, clone_key$Protocol_name)
counts_tbl <- data.frame(
  Protocol_name         = burden$sample,
  Original_clone_ID     = rev_key[burden$sample],
  SNVs_before_blacklist = as.integer(lengths(grl_raw)[burden$sample]),
  SNVs_after_blacklist  = as.integer(burden$n_filtered), row.names = NULL)
counts_tbl$Removed_pct <- round(100 * (1 - counts_tbl$SNVs_after_blacklist /
                                           counts_tbl$SNVs_before_blacklist), 2)
counts_tbl <- counts_tbl[order(as.integer(factor(group_of(counts_tbl$Protocol_name),
                          levels = c("NT","UV","CP","BAP","ENU"))),
                          as.integer(sub("^[A-Z]+", "", counts_tbl$Protocol_name))), ]
write.csv(counts_tbl, file.path(CFG$outdir, "Table_S3_mutation_counts.csv"), row.names = FALSE)

# mut_matrix
mut_mat  <- mut_matrix(vcf_list = grl, ref_genome = CFG$ref_genome)
stopifnot(nrow(mut_mat) == 96L)
mut_prop <- sweep(mut_mat, 2, colSums(mut_mat), "/")
grp <- factor(group_of(colnames(mut_prop)), levels = c("NT","UV","CP","BAP","ENU"))
mut_prop_grp <- sapply(levels(grp), function(g) rowMeans(mut_prop[, grp == g, drop = FALSE]))

# reference catalogues
signal_exp <- get_known_signatures(muttype = "snv", source = "SIGNAL", sig_type = "exposure")
cosmic_sbs <- get_known_signatures(muttype = "snv", source = "COSMIC", genome = "GRCh38")
stopifnot(nrow(signal_exp) == 96L, nrow(cosmic_sbs) == 96L)
if (is.null(rownames(signal_exp))) rownames(signal_exp) <- rownames(mut_mat)
if (is.null(rownames(cosmic_sbs))) rownames(cosmic_sbs) <- rownames(mut_mat)
stopifnot(identical(rownames(mut_mat), rownames(signal_exp)),
          identical(rownames(mut_mat), rownames(cosmic_sbs)))
stopifnot(ncol(cosmic_sbs) == 60L)
message(sprintf(">>> COSMIC signatures returned by get_known_signatures(): %d  (use in caption)",
                ncol(cosmic_sbs)))

resolve_ref <- function(tokens, cols) { for (t in tokens) cols <- cols[grepl(t, cols, ignore.case = TRUE)]; cols }
prespec_native <- vapply(names(prespec_ref), function(lbl) {
  hit <- resolve_ref(prespec_ref[[lbl]], colnames(signal_exp))
  if (length(hit) != 1L)
    stop(sprintf("Reference '%s' matched %d SIGNAL columns (%s). Refine its tokens.",
                 lbl, length(hit), paste(hit, collapse = " | ")))
  hit
}, character(1))
message(">>> targeted references resolved: ",
        paste(sprintf("%s <- %s", names(prespec_native), prespec_native), collapse = " ; "))


cos_signal   <- cos_sim_matrix(mut_mat, signal_exp)
cos_targeted <- cos_signal[, prespec_native, drop = FALSE]
colnames(cos_targeted) <- names(prespec_native)

best_native <- colnames(signal_exp)[max.col(cos_signal, ties.method = "first")]
names(best_native) <- rownames(cos_signal)
all_expected <- c(
  UV = all(best_native[burden$sample[burden$group == "UV"]] == prespec_native[["SSR (1.25 J)"]]),
  ENU = all(best_native[burden$sample[burden$group == "ENU"]] == prespec_native[["ENU (400 µM)"]]),
  BAP_BPDE = all(best_native[burden$sample[burden$group == "BAP"]] == prespec_native[["BPDE (0.125 µM)"]])
)
message(">>> all-clone best matches: ", paste(names(all_expected), all_expected, sep = "=", collapse = " ; "))

cp_samples <- burden$sample[burden$group == "CP"]
message(sprintf(">>> CP on-target cosine (median [range]): %.3f [%.3f - %.3f]",
                median(cos_targeted[cp_samples, "Cisplatin (12.5 µM)"]),
                min(cos_targeted[cp_samples, "Cisplatin (12.5 µM)"]),
                max(cos_targeted[cp_samples, "Cisplatin (12.5 µM)"])))

# save plots
save3 <- function(p, name, w, h) {
  f <- file.path(CFG$outdir, name)
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(f, ".png"), p, width = w, height = h, dpi = 300, device = "png", type = "cairo")
  svg(paste0(f, ".svg"), width = w, height = h); print(p); invisible(dev.off())
}

# Figure 2
ord <- c("NT","UV","CP","BAP","ENU")
clone_number <- function(x) as.integer(sub("^[A-Z]+", "", x))
sample_order <- unlist(lapply(ord, function(g) {
  x <- burden$sample[burden$group == g]
  x[order(clone_number(x))]
}), use.names = FALSE)
stopifnot(setequal(sample_order, rownames(cos_targeted)))
cos_targeted_plot <- cos_targeted[sample_order, , drop = FALSE]

pA <- ggplot(burden, aes(x = group, y = n_filtered)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85", colour = "grey30") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, colour = "grey20") +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  labs(x = "Treatment", y = "SNVs per clone (log scale)") + theme_minimal()
pB <- plot_96_profile(as.data.frame(mut_prop_grp[, ord]), condensed = FALSE)
pC <- plot_cosine_heatmap(cos_targeted_plot, cluster_rows = FALSE, plot_values = TRUE)
save3(pA, "fig2A_burden",  5, 4)
save3(pB, "fig2B_spectra", 11, 7)
save3(pC, "fig2C_cosine",  6, 8)

# Figure S1
THR <- 0.6
grp_mat <- as.matrix(mut_prop_grp[, ord])

# Display-only lookup for Figure S1. All cosine calculations above use native SIGNAL column names.
# The lookup is validated against all 54 columns and is not used for reference resolution.
SIGNAL_NAMES <- c(
  "Potassium.bromate..875.uM."="Potassium bromate (875 µM)","DBADE..0.109.uM."="DBADE (0.109 µM)",
  "Formaldehyde..120.uM."="Formaldehyde (120 µM)","Semustine..150.uM."="Semustine (150 µM)",
  "Temozolomide..200.uM."="Temozolomide (200 µM)","DMH..11.6.mM....S9"="DMH (11.6 mM) +S9",
  "Benzidine..200.uM."="Benzidine (200 µM)","DBP..0.0039.uM."="DBP (0.0039 µM)",
  "MX..7.uM....S9"="MX (7 µM) +S9","Methyleugenol..1.25.mM."="Methyleugenol (1.25 mM)",
  "X4.ABP..300.uM....S9"="4-ABP (300 µM) +S9","DBPDE..0.000156.uM."="DBPDE (0.000156 µM)",
  "DBP..0.0313.uM....S9"="DBP (0.0313 µM) +S9","DBADE..0.0313.uM."="DBADE (0.0313 µM)",
  "X1.8.DNP..0.125.uM."="1,8-DNP (0.125 µM)","BPDE..0.125.uM."="BPDE (0.125 µM)",
  "MNU..350.uM."="MNU (350 µM)","ENU..400.uM."="ENU (400 µM)",
  "Cyclophosphamide..18.75.uM....S9"="Cyclophosphamide (18.75 µM) +S9","BaP..0.39.uM....S9"="BaP (0.39 µM) +S9",
  "X6.Nitrochrysene..12.5.uM....S9"="6-Nitrochrysene (12.5 µM) +S9","AAI..1.25.uM."="AAI (1.25 µM)",
  "Potassium.bromate..260.uM."="Potassium bromate (260 µM)","X6.Nitrochrysene..0.78.uM."="6-Nitrochrysene (0.78 µM)",
  "Ellipticine..0.375.uM....S9"="Ellipticine (0.375 µM) +S9","DBA..75.uM....S9"="DBA (75 µM) +S9",
  "PhIP..3.uM....S9"="PhIP (3 µM) +S9","AFB1..0.25.uM....S9"="AFB1 (0.25 µM) +S9",
  "X3.NBA..0.025.uM."="3-NBA (0.025 µM)","X1.6.DNP..0.09.uM."="1,6-DNP (0.09 µM)",
  "X5.Methylchrysene..1.6.uM....S9"="5-Methylchrysene (1.6 µM) +S9","Furan..100.mM....S9"="Furan (100 mM) +S9",
  "SSR..1.25.J."="SSR (1.25 J)","AAII..37.5.uM."="AAII (37.5 µM)",
  "Propylene.oxide..10.mM."="Propylene oxide (10 mM)","N.Nitrosopyrrolidine..50.mM."="N-Nitrosopyrrolidine (50 mM)",
  "Mechlorethamine..0.3.uM."="Mechlorethamine (0.3 µM)","DES..0.938.mM."="DES (0.938 mM)",
  "DMS..0.078.mM."="DMS (0.078 mM)","Cisplatin..3.125.uM."="Cisplatin (3.125 µM)",
  "OTA..0.08.uM....S9"="OTA (0.08 µM) +S9","Carboplatin..5.uM."="Carboplatin (5 µM)",
  "DBAC..5.uM....S9"="DBAC (5 µM) +S9","Temozolomide..200.uM..1"="Temozolomide (200 µM) [2]",
  "Cisplatin..12.5.uM."="Cisplatin (12.5 µM)","AZD7762..1.625.uM."="AZD7762 (1.625 µM)",
  "X3.NBA..0.1.uM."="3-NBA (0.1 µM)","PhIP..4.uM....S9"="PhIP (4 µM) +S9",
  "BaP..2.uM....S9"="BaP (2 µM) +S9","X6.Nitrochrysene..50.uM....S9"="6-Nitrochrysene (50 µM) +S9",
  "X6.Nitrochrysene..50.uM."="6-Nitrochrysene (50 µM)","X1.8.DNP..8.uM."="1,8-DNP (8 µM)",
  "DBPDE..0.000625.uM."="DBPDE (0.000625 µM)","Control"="Control")
signal_disp <- signal_exp
stopifnot(ncol(signal_disp) == 54L,
          setequal(colnames(signal_disp), names(SIGNAL_NAMES)))
colnames(signal_disp) <- unname(SIGNAL_NAMES[colnames(signal_disp)])

cosA <- cos_sim_matrix(cosmic_sbs,  grp_mat)   # COSMIC, 60 signatures
cosB <- cos_sim_matrix(signal_disp, grp_mat)   # SIGNAL, 54 signatures
pal <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
nt <- length(ord)
Axh <- 1:nt; Adnd <- c(nt + 0.6, nt + 3.6); Bdnd <- c(nt + 8.5, nt + 11.5)
Bxh <- Bdnd[2] + 0.5 + (1:nt); Alab_x <- 0.4; Blab_x <- max(Bxh) + 0.6
TOP <- nrow(cosA)  

build <- function(cm, xcols, dndband, dnd_dir) {
  hc <- hclust(dist(cm), method = "complete")
  pos <- setNames(seq_along(hc$order), rownames(cm)[hc$order])
  yof <- function(sig) TOP - pos[sig] + 1
  tl <- expand.grid(sig = rownames(cm), grp = colnames(cm), stringsAsFactors = FALSE)
  tl$x <- xcols[match(tl$grp, colnames(cm))]; tl$y <- yof(tl$sig)
  tl$v <- mapply(function(s, g) cm[s, g], tl$sig, tl$grp)
  tl$lab <- ifelse(tl$v >= THR, sprintf("%.2f", tl$v), "")
  tl$tcol <- ifelse(tl$v >= 0.78, "white", "grey15")
  rl <- data.frame(sig = rownames(cm), y = yof(rownames(cm)))
  seg <- ggdendro::dendro_data(as.dendrogram(hc))$segments; mh <- max(seg$y)
  fx <- if (dnd_dir == "right") function(h) dndband[1] + h/mh * diff(dndband)
        else                    function(h) dndband[2] - h/mh * diff(dndband)
  seg <- transform(seg, X = fx(y), Xend = fx(yend), Y = TOP - x + 1, Yend = TOP - xend + 1)
  list(tiles = tl, rows = rl, seg = seg)
}
A <- build(cosA, Axh, Adnd, "right"); B <- build(cosB, Bxh, Bdnd, "left")
colA <- data.frame(x = Axh, y = 0.0, lab = ord)
colB <- data.frame(x = Bxh, y = TOP - nrow(cosB) - 1, lab = ord)
pS1 <- ggplot() +
  geom_tile(data = A$tiles, aes(x, y, fill = v), width = 0.96, height = 0.9) +
  geom_tile(data = B$tiles, aes(x, y, fill = v), width = 0.96, height = 0.9) +
  scale_fill_gradientn(colours = pal, limits = c(0, 1), name = "Cosine\nsimilarity",
                       guide = guide_colorbar(barheight = 9, barwidth = 1)) +
  geom_text(data = subset(A$tiles, lab != ""), aes(x, y, label = lab, colour = tcol), size = 2.0) +
  geom_text(data = subset(B$tiles, lab != ""), aes(x, y, label = lab, colour = tcol), size = 2.0) +
  scale_colour_identity() +
  geom_segment(data = A$seg, aes(x = X, y = Y, xend = Xend, yend = Yend), linewidth = 0.25, colour = "grey35") +
  geom_segment(data = B$seg, aes(x = X, y = Y, xend = Xend, yend = Yend), linewidth = 0.25, colour = "grey35") +
  geom_text(data = A$rows, aes(Alab_x, y, label = sig), hjust = 1, size = 2.1) +
  geom_text(data = B$rows, aes(Blab_x, y, label = sig), hjust = 0, size = 2.1) +
  geom_text(data = colA, aes(x, y, label = lab), size = 3.1, fontface = "bold", vjust = 1) +
  geom_text(data = colB, aes(x, y, label = lab), size = 3.1, fontface = "bold", vjust = 1) +
  annotate("text", x = mean(Axh), y = TOP + 3, label = "COSMIC (reference)", fontface = "bold", size = 4) +
  annotate("text", x = mean(Bxh), y = TOP + 3, label = "SIGNAL (exposure)", fontface = "bold", size = 4) +
  annotate("text", x = Alab_x - 3.5, y = TOP + 3, label = "A", fontface = "bold", size = 6) +
  annotate("text", x = Bdnd[1] - 0.3, y = TOP + 3, label = "B", fontface = "bold", size = 6) +
  coord_cartesian(xlim = c(Alab_x - 4, Blab_x + 10), ylim = c(-2, TOP + 4), expand = FALSE, clip = "off") +
  theme_void() +
  theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        plot.margin = margin(6, 6, 6, 6))
xlo <- Alab_x - 4; xhi <- Blab_x + 10
pS1 <- pS1 + theme(legend.position = c((mean(c(Adnd[2], Bdnd[1])) - xlo) / (xhi - xlo), 0.55))
ggsave(file.path(CFG$outdir, "figS1_cosmic_signal.pdf"), pS1,
       width = 12, height = 11, bg = "white", device = cairo_pdf, limitsize = FALSE)

message("Done. Outputs in ", normalizePath(CFG$outdir))
