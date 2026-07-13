## 1) Define low/high sample sets
## -----------------------------
ch4net_low_samples <- c(
  "SRR873595","SRR873596","SRR873597","SRR873598",
  "SRR873599","SRR873600","SRR873601","SRR873602"
)

ch4net_high_samples <- c(
  "SRR1206671","SRR873604","SRR873605","SRR873606",
  "SRR873607","SRR873608","SRR873609","SRR873610"
)

keep_srr <- c(ch4net_low_samples, ch4net_high_samples)

## -----------------------------
## 2) Clean + filter to low/high + remove NA pN/pS
## -----------------------------
dt <- copy(ko_pnps_df)

# Ensure numeric pN/pS and remove NAs/infinite
dt[, pN_pS := as.numeric(pN_pS)]
dt <- dt[!is.na(pN_pS) & is.finite(pN_pS)]

# Normalize SRR IDs (handles SRRxxxx.1)
dt[, SRR_base := sub("\\.\\d+$", "", SRR)]

# Keep only the 16 samples of interest
dt <- dt[SRR_base %in% keep_srr]

# Label condition
dt[, condition := fifelse(SRR_base %in% ch4net_low_samples, "low",
                          fifelse(SRR_base %in% ch4net_high_samples, "high", NA_character_))]
dt <- dt[!is.na(condition)]

## -----------------------------
## 3) Enforce exactly 8 low + 8 high per genome×KO
## -----------------------------
counts <- dt[, .(
  n_low  = uniqueN(SRR_base[condition == "low"]),
  n_high = uniqueN(SRR_base[condition == "high"])
), by = .(genome, KO_number)]

eligible_pairs <- counts[n_low == 8 & n_high == 8, .(genome, KO_number)]

# Subset to eligible genome×KO pairs
dt_elig <- dt[eligible_pairs, on = .(genome, KO_number)]

## Optional sanity checks
# How many eligible pairs?
cat("Eligible genome×KO pairs:", nrow(eligible_pairs), "\n")
# How many genomes represented?
cat("Eligible genomes:", uniqueN(dt_elig$genome), "\n")

## -----------------------------
## 4) Parallel Wilcoxon by genome (compute p-values per KO)
## -----------------------------
plan(multisession, workers = 10)

# Split into a list by genome for parallel processing
split_list <- split(dt_elig, by = "genome", keep.by = TRUE)

wilcox_raw <- rbindlist(
  future_lapply(split_list, function(dg) {
    
    # dg = one genome
    out <- dg[, {
      x <- pN_pS[condition == "low"]
      y <- pN_pS[condition == "high"]
      
      p <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
      
      .(
        n_low = length(x),
        n_high = length(y),
        median_low = median(x),
        median_high = median(y),
        p_value = p
      )
    }, by = KO_number]
    
    out[, genome := unique(dg$genome)]
    out
  }),
  fill = TRUE
)

## -----------------------------
## 5) BH correction WITHIN each genome (not global)
## -----------------------------
wilcox_raw[, q_value_genome := p.adjust(p_value, method = "BH"), by = genome]

## Direction label
wilcox_raw[, direction := fifelse(median_high > median_low, "High", "Low")]

## Order results
setorder(wilcox_raw, genome, q_value_genome)

## -----------------------------
## 6) Results
## -----------------------------
# View top hits overall
print(wilcox_raw[order(q_value_genome)][1:20])

# Example: top hits for a single genome
# print(wilcox_raw[genome == "Fibrobacter succinogenes subsp. succinogenes S85"][order(q_value_genome)][1:20])

## Optional: save output
# fwrite(wilcox_raw, "/Users/stephencourtney/Desktop/wilcox_pnps_byGenome_byKO.tsv", sep = "\t")

wilcox_raw


# keep genome-level significant KO tests
sig_wilcox <- wilcox_raw[q_value_genome < 0.05]

# how many?
nrow(sig_wilcox)

# look at top hits
sig_wilcox[order(q_value_genome)][1:20]


sig_wilcox2 <- wilcox_raw[!is.na(q_value_genome) & q_value_genome < 0.05]

fwrite(sig_wilcox,
       "/Users/stephencourtney/Desktop/wilcox_pnps_byGenome_byKO_significant_q0.05.tsv",
       sep = "\t")


library(data.table)
setDT(sig_wilcox)

# 1) Number of unique genomes with ≥1 significant KO
n_sig_genomes <- uniqueN(sig_wilcox$genome)
n_sig_genomes

# 2) Number of significant KOs per genome
sig_kos_per_genome <- sig_wilcox[, .(n_sig_KOs = uniqueN(KO_number)), by = genome][
  order(-n_sig_KOs)
]

sig_kos_per_genome




library(data.table)
setDT(module_dfmod_clean)

# number of unique genomes with at least one significant module
n_genomes_with_sig_module <- module_dfmod_clean[significant == 1, uniqueN(genome)]
n_genomes_with_sig_module




library(data.table)

setDT(module_dfmod_clean)
setDT(sig_wilcox)

# 1) Genomes significant for modules
genomes_sig_module <- unique(module_dfmod_clean[significant == 1, genome])

# 2) Genomes significant for KO pN/pS
genomes_sig_ko <- unique(sig_wilcox[q_value_genome < 0.05, genome])

# Overlap / intersection
genomes_sig_both <- intersect(genomes_sig_module, genomes_sig_ko)

# Summaries
cat("Significant module genomes:", length(genomes_sig_module), "\n")
cat("Significant KO pN/pS genomes:", length(genomes_sig_ko), "\n")
cat("Significant in BOTH:", length(genomes_sig_both), "\n")

# Are all module-significant genomes also KO-significant?
all_module_in_ko <- all(genomes_sig_module %in% genomes_sig_ko)
cat("All sig-module genomes also sig-KO?", all_module_in_ko, "\n")

# Who is missing from which set?
module_only <- setdiff(genomes_sig_module, genomes_sig_ko)
ko_only <- setdiff(genomes_sig_ko, genomes_sig_module)

cat("Module-only genomes:", length(module_only), "\n")
cat("KO-only genomes:", length(ko_only), "\n")

# View them if you want
module_only
ko_only
genomes_sig_both




sig_modules_per_genome <- module_dfmod_clean[significant == 1,
                                             .(n_sig_modules = uniqueN(Modules)), by = genome
]

sig_kos_per_genome <- sig_wilcox[q_value_genome < 0.05,
                                 .(n_sig_KOs = uniqueN(KO_number)), by = genome
]

both_summary <- merge(sig_modules_per_genome, sig_kos_per_genome, by = "genome", all = FALSE)
both_summary <- both_summary[genome %in% genomes_sig_both][order(-n_sig_modules, -n_sig_KOs)]

both_summary





suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDT(sig_wilcox)
setDT(module_dfmod_clean)

# --- KO-level: count significant KOs by genome + direction
ko_sig <- sig_wilcox[q_value_genome < 0.05,
                     .(n_sig = uniqueN(KO_number)),
                     by = .(genome, direction)
]

ko_sig[, test_type := "KO"]

# Ensure direction labels match your caption (HME vs LME)
# In your wilcox code: direction=="High" means median_high > median_low
# Map that to HME/LME wording if needed:
ko_sig[, direction := fifelse(direction == "High", "HME", "LME")]

# --- Module-level: count significant modules by genome + direction
mod_sig <- module_dfmod_clean[significant == 1,
                              .(n_sig = uniqueN(Modules)),
                              by = .(genome, median_comparison)
]
mod_sig[, test_type := "Module"]
setnames(mod_sig, "median_comparison", "direction")

# Map direction to HME/LME if your column uses "High"/"Low"
mod_sig[, direction := fifelse(direction == "High", "HME",
                               fifelse(direction == "Low",  "LME", direction))]

# --- Combine
all_sig <- rbindlist(list(ko_sig, mod_sig), fill = TRUE)

# Fill missing genome×type×direction combinations with 0
all_sig <- dcast(
  all_sig,
  genome + test_type ~ direction,
  value.var = "n_sig",
  fill = 0
)

# Long format for plotting (two columns per type: HME and LME)
plot_dt <- melt(
  all_sig,
  id.vars = c("genome","test_type"),
  variable.name = "direction",
  value.name = "n_sig"
)

# Order genomes by total significant results (KO + module, both directions)
genome_order <- plot_dt[, .(total = sum(n_sig)), by = genome][order(-total), genome]
plot_dt[, genome := factor(genome, levels = genome_order)]

# Make a signed value for red/blue direction:
# positive (HME) -> red, negative (LME) -> blue
plot_dt[, signed_n := ifelse(direction == "HME", n_sig, -n_sig)]

# --- Plot
p <- ggplot(plot_dt, aes(x = test_type, y = genome, fill = signed_n)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Significant tests\n(+HME / −LME)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7)
  )

print(p)

# Optional save
# ggsave("Figure1_sig_KO_module_heatmap.png", p, width = 7, height = 10, dpi = 300)



suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setDT(sig_wilcox)
setDT(module_dfmod_clean)

alpha <- 0.05

# -----------------------------
# 1) KO-level significant counts by genome + direction
# -----------------------------
ko_sig <- sig_wilcox[q_value_genome < alpha,
                     .(n_sig = uniqueN(KO_number)),
                     by = .(genome, direction)
]
ko_sig[, test_type := "KO"]

# Map direction to HME/LME
# direction == "High" means median_high > median_low => higher pN/pS in HME
ko_sig[, direction := fifelse(direction == "High", "HME", "LME")]

# -----------------------------
# 2) Module-level significant counts by genome + direction
# -----------------------------
mod_sig <- module_dfmod_clean[significant == 1,
                              .(n_sig = uniqueN(Modules)),
                              by = .(genome, median_comparison)
]
mod_sig[, test_type := "Module"]
setnames(mod_sig, "median_comparison", "direction")

# Map direction to HME/LME if needed
mod_sig[, direction := fifelse(direction == "High", "HME",
                               fifelse(direction == "Low", "LME", direction))]

# -----------------------------
# 3) Combine and fill missing combos with 0
# -----------------------------
all_sig <- rbindlist(list(ko_sig, mod_sig), fill = TRUE)

# Ensure every genome has both directions for both test types (fill 0)
all_sig <- dcast(all_sig, genome + test_type ~ direction, value.var = "n_sig", fill = 0)
all_sig <- melt(all_sig,
                id.vars = c("genome", "test_type"),
                variable.name = "direction",
                value.name = "n_sig"
)

# -----------------------------
# 4) Order genomes by total significant results
# -----------------------------
genome_order <- all_sig[, .(total = sum(n_sig)), by = genome][order(-total), genome]
all_sig[, genome := factor(genome, levels = genome_order)]

# (Optional) order bars KO then Module within each genome
all_sig[, test_type := factor(test_type, levels = c("KO", "Module"))]

# -----------------------------
# 5) Plot: stacked bars, two bars per genome (KO vs Module)
# -----------------------------
p <- ggplot(all_sig, aes(x = genome, y = n_sig, fill = direction)) +
  geom_col(width = 0.8) +
  facet_wrap(~ test_type, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Number of significant tests", fill = "Direction") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )

print(p)

# Optional save
# ggsave("Figure1_stackedbars_KO_Module_byGenome.png", p, width = 12, height = 5, dpi = 300)






sig_KOs_direction <- sig_wilcox[q_value_genome < 0.05,
                                .(n_sig_KOs = uniqueN(KO_number)),
                                by = .(genome, direction)
][order(genome, direction)]

sig_KOs_direction





'/Users/stephencourtney/Library/Mobile Documents/com~apple~CloudDocs/analysis_folder/MG_organism_DESeq2_High_vs_Low_filter10.tsv'
'/Users/stephencourtney/Library/Mobile Documents/com~apple~CloudDocs/analysis_folder/MT_organism_DESeq2_High_vs_Low_filter10.tsv'