############################
# G4Hunter Conservation Pipeline
# Robust final version
############################

setwd('/home/dimitri/cedarHazaraToscana/others/ebola/testAlternativeApproach/')
getwd()

############################
# 0. LIBRARIES
############################

library(Biostrings)
library(dplyr)
library(stringr)
library(seqinr)

############################
# 1. LOAD DATA
############################

message("=== Loading data ===")


# align all genomes first using mafft
# Aligned genomes (MAFFT output)
aln <- readDNAStringSet("all_genomes.aln")

# this is an original .fasta file that contains all the COMPLETE GENOMES
# Original (unaligned) genomes — used for sequence extraction
genomes <- readDNAStringSet("zaire_ebola_completeGnms.fa")

# this is a table that is created using the G4hunter per genome
# G4Hunter output (skip header comment line)
g4 <- read.csv("zaire_genome_G4Hunter.csv", skip = 1)

message("Loaded: ", length(aln), " aligned sequences")
message("Loaded: ", length(genomes), " original genomes")
message("Loaded: ", nrow(g4), " G4Hunter hits")

############################
# 2. CLEAN & VALIDATE NAMES
############################

message("\n=== Cleaning sequence names ===")

# Strip path-like prefix added by G4Hunter (e.g. "completeGnms.fa_ACC")
g4$seqnames <- sub("^zaire_ebola_completeGnms\\.fa_", "", g4$seqnames)

# Keep only the accession (everything before the first space) in all three objects
names(aln)     <- sub(" .*", "", names(aln))
names(genomes) <- sub(" .*", "", names(genomes))
g4$seqnames    <- sub(" .*", "", g4$seqnames)

# --- Validation checks ---
missing_in_aln <- setdiff(unique(g4$seqnames), names(aln))
if (length(missing_in_aln) > 0) {
  warning(length(missing_in_aln), " G4 seqname(s) not found in alignment — they will be dropped:\n  ",
          paste(missing_in_aln, collapse = ", "))
}

missing_in_genomes <- setdiff(unique(g4$seqnames), names(genomes))
if (length(missing_in_genomes) > 0) {
  warning(length(missing_in_genomes), " G4 seqname(s) not found in original genomes — flanking extraction will fail:\n  ",
          paste(missing_in_genomes, collapse = ", "))
}

aln_not_in_g4 <- setdiff(names(aln), unique(g4$seqnames))
if (length(aln_not_in_g4) > 0) {
  message("Note: ", length(aln_not_in_g4), " aligned genome(s) have no G4 hits.")
}

############################
# 3. FILTER BY G4HUNTER SCORE
############################

message("\n=== Filtering by G4Hunter score ===")

# Retain only high-confidence hits (|score| >= 1.5 is a common threshold)
G4_SCORE_THRESHOLD <- 1.5

n_before <- nrow(g4)
g4 <- g4 %>% filter(abs(score) >= G4_SCORE_THRESHOLD)
message("Retained ", nrow(g4), " / ", n_before, " hits with |score| >= ", G4_SCORE_THRESHOLD)

if (nrow(g4) == 0) stop("No G4 hits remain after score filtering. Check your threshold or input file.")

############################
# 4. BUILD GENOME → ALIGNMENT COORDINATE MAPS
############################

message("\n=== Building coordinate maps ===")

build_coord_map <- function(aligned_seq) {
  chars <- strsplit(as.character(aligned_seq), "")[[1]]
  genome_pos <- 0L
  map <- integer(length(chars))
  for (i in seq_along(chars)) {
    if (chars[i] != "-") {
      genome_pos <- genome_pos + 1L
      map[i] <- genome_pos
    } else {
      map[i] <- NA_integer_
    }
  }
  return(map)
}

# Only build maps for sequences actually present in g4
seqs_needed <- intersect(unique(g4$seqnames), names(aln))
coord_maps  <- lapply(aln[seqs_needed], build_coord_map)
message("Built coordinate maps for ", length(coord_maps), " sequences")

############################
# 5. MAP G4 POSITIONS → ALIGNMENT COORDINATES
############################

message("\n=== Mapping G4 positions to alignment coordinates ===")

g4$aln_start <- NA_integer_
g4$aln_end   <- NA_integer_

for (i in seq_len(nrow(g4))) {
  seqname <- g4$seqnames[i]
  map     <- coord_maps[[seqname]]
  if (is.null(map)) next
  
  hits_start <- which(map == g4$start[i])
  hits_end   <- which(map == g4$end[i])
  
  if (length(hits_start) == 0 || length(hits_end) == 0) next
  
  g4$aln_start[i] <- hits_start[1]
  g4$aln_end[i]   <- hits_end[1]
}

n_failed <- sum(is.na(g4$aln_start) | is.na(g4$aln_end))
if (n_failed > 0) warning(n_failed, " G4 hit(s) could not be mapped to alignment coordinates and will be dropped.")

g4 <- g4 %>% filter(!is.na(aln_start), !is.na(aln_end))
message(nrow(g4), " hits successfully mapped")

############################
# 6. CLUSTER IN ALIGNMENT SPACE
############################

message("\n=== Clustering G4s in alignment space ===")

CLUSTER_GAP <- 10L   # merge hits within this many alignment columns of each other

g4 <- g4 %>% arrange(aln_start)

cluster_id  <- 1L
max_end_seen <- g4$aln_end[1]
clusters    <- integer(nrow(g4))
clusters[1] <- cluster_id

for (i in 2:nrow(g4)) {
  if (g4$aln_start[i] <= max_end_seen + CLUSTER_GAP) {
    # extend current cluster
    clusters[i]  <- cluster_id
    max_end_seen <- max(max_end_seen, g4$aln_end[i])
  } else {
    # start a new cluster
    cluster_id   <- cluster_id + 1L
    clusters[i]  <- cluster_id
    max_end_seen <- g4$aln_end[i]
  }
}

g4$cluster <- clusters
message("Identified ", max(clusters), " clusters")

############################
# 7. CONSERVATION PER CLUSTER
############################

message("\n=== Computing conservation ===")

n_genomes <- length(aln)

cluster_summary <- g4 %>%
  group_by(cluster) %>%
  summarise(
    aln_start_median  = median(aln_start),
    aln_end_median    = median(aln_end),
    genomes_with_g4   = n_distinct(seqnames),
    conservation      = genomes_with_g4 / n_genomes,
    mean_g4_length    = mean(end - start + 1),
    mean_score        = mean(abs(score)),
    max_score         = max(abs(score)),
    .groups = "drop"
  ) %>%
  arrange(desc(conservation), desc(mean_score))

message("Conservation summary (top 10 clusters):")
print(head(cluster_summary, 10))

############################
# 8. EXTRACT FLANKING SEQUENCES
############################

message("\n=== Extracting flanking sequences ===")

FLANK_SIZE <- 30L

g4$flank_seq <- NA_character_

for (i in seq_len(nrow(g4))) {
  seqname <- g4$seqnames[i]
  seq     <- genomes[[seqname]]
  if (is.null(seq)) next
  
  seq_len <- length(seq)
  left    <- max(1L, g4$start[i] - FLANK_SIZE)
  right   <- min(seq_len, g4$end[i] + FLANK_SIZE)
  
  g4$flank_seq[i] <- as.character(subseq(seq, left, right))
}

n_missing_flank <- sum(is.na(g4$flank_seq))
if (n_missing_flank > 0) warning(n_missing_flank, " hit(s) have no flanking sequence (seqname not found in genomes).")

############################
# 9. EXPORT CLUSTER FASTA FILES
############################

message("\n=== Exporting per-cluster FASTA files ===")

dir.create("G4_clusters_flanked", showWarnings = FALSE)

for (cl in sort(unique(g4$cluster))) {
  sub  <- g4 %>% filter(cluster == cl, !is.na(flank_seq))
  if (nrow(sub) == 0) next
  
  seqs        <- as.list(sub$flank_seq)
  seq_names   <- paste0(sub$seqnames, "_", sub$start, "_", sub$end)
  
  write.fasta(
    sequences = seqs,
    names     = seq_names,
    file.out  = paste0("G4_clusters_flanked/G4_cluster_", cl, ".fa")
  )
}

message("FASTA files written to G4_clusters_flanked/")

############################
# 10. FILTER HIGHLY CONSERVED G4s
############################

message("\n=== Filtering highly conserved G4s ===")

CONSERVATION_THRESHOLD <- 0.8

high_conf <- g4 %>%
  inner_join(cluster_summary, by = "cluster") %>%
  filter(conservation >= CONSERVATION_THRESHOLD) %>%
  arrange(desc(conservation), desc(mean_score))

message(nrow(high_conf), " hits in clusters with conservation >= ", CONSERVATION_THRESHOLD)

############################
# 11. G4 DENSITY PER GENOME
############################

message("\n=== Computing G4 density per genome ===")

genome_lengths        <- width(genomes)
names(genome_lengths) <- names(genomes)

g4_density <- g4 %>%
  group_by(seqnames) %>%
  summarise(
    genome_length    = genome_lengths[seqnames[1]],
    g4_count         = n(),
    density_per_kb   = g4_count / genome_length * 1000,
    mean_score       = mean(abs(score)),
    .groups = "drop"
  ) %>%
  arrange(desc(density_per_kb))

if (any(is.na(g4_density$genome_length))) {
  warning("Some genomes have NA length — check that seqnames match between g4 and genomes.")
}

message("Density summary (top 10 genomes):")
print(head(g4_density, 10))

############################
# 12. EXPORT RESULTS
############################

message("\n=== Exporting results ===")

write.csv(cluster_summary, "G4_cluster_summary.csv",    row.names = FALSE)
write.csv(high_conf,       "G4_highly_conserved.csv",   row.names = FALSE)
write.csv(g4_density,      "G4_density.csv",            row.names = FALSE)
write.csv(g4,              "G4_all_hits_annotated.csv", row.names = FALSE)

message("\n=== Pipeline complete ===")
message("Output files:")
message("  G4_cluster_summary.csv      — per-cluster conservation stats")
message("  G4_highly_conserved.csv     — hits in clusters with conservation >= ", CONSERVATION_THRESHOLD)
message("  G4_density.csv              — G4 density per genome")
message("  G4_all_hits_annotated.csv   — all hits with cluster & alignment coordinates")
message("  G4_clusters_flanked/        — per-cluster FASTA files with flanking sequences")

