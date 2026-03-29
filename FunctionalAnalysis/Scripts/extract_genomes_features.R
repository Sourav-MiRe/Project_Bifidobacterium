#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# ========================= Options / flags =========================
VERBOSE <- TRUE
base_dir <- "/storage/alisha/molecular_mimicry"
fasta_dir <- file.path(base_dir, "protein_fasta_66_species_dump")
anno_dir  <- file.path(base_dir, "protein_annotation_66_species_dump")
out_dir   <- file.path(base_dir, "Feature_fasta")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========================= Logger helper =========================
log_msg <- function(fmt, ...) {
  if (!VERBOSE) return(invisible(NULL))
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, sprintf(fmt, ...)))
}

# ===================== Normalization =====================
normalize_token <- function(x) {
  if (is.null(x)) return(character(0))
  x <- unlist(strsplit(x, ",", fixed = TRUE))   # split comma separated
  x <- gsub("^ko:", "", x, ignore.case = TRUE)  # drop ko: prefix
  x <- sub("@.*", "", x)                        # drop everything after @
  toupper(trimws(x[nzchar(x)]))                 # clean + uppercase
}

# ================= Genome ↔ Species map =================
map_path <- file.path(base_dir, "genome_species_66.csv")
if (!file.exists(map_path)) stop("Missing genome map: ", map_path)
genome_map <- fread(map_path, header = FALSE)
setnames(genome_map, c("genome_id", "species"))
genome_to_species <- setNames(genome_map$species, genome_map$genome_id)
log_msg("Loaded genome map: %d rows.", nrow(genome_map))

# ============== Features + species list =================
feat_path <- file.path(base_dir, "annotations_196_featureid.csv")
if (!file.exists(feat_path)) stop("Missing features file: ", feat_path)
features_df <- fread(feat_path, header = TRUE, showProgress = FALSE)
features_df <- features_df[
  !is.na(Feature_ids) & Feature_ids != "" &
    !is.na(Detection_NonBif_list) & Detection_NonBif_list != ""
]
feature_list <- unique(trimws(features_df$Feature_ids))

log_msg("Loaded features table: %d valid rows; %d unique Feature IDs.",
        nrow(features_df), length(feature_list))

feature_to_species <- lapply(feature_list, function(fid) {
  x <- trimws(unlist(strsplit(
    features_df[Feature_ids == fid, Detection_NonBif_list],
    ";", fixed = TRUE
  )))
  unique(x[nzchar(x)])
})
names(feature_to_species) <- feature_list

species_to_features <- list()
for (fid in names(feature_to_species)) {
  for (sp in feature_to_species[[fid]]) {
    if (!nzchar(sp)) next
    species_to_features[[sp]] <- unique(c(species_to_features[[sp]], fid))
  }
}

needed_species <- names(species_to_features)
needed_genomes <- genome_map[species %in% needed_species, unique(genome_id)]
log_msg("Species with features: %d | Genomes to process: %d",
        length(needed_species), length(needed_genomes))

# Prepare empty files
for (fid in feature_list) {
  fp <- file.path(out_dir, paste0(fid, ".fasta"))
  if (!file.exists(fp)) file.create(fp)
}

# ================= Matching settings =================
search_cols <- c("eggNOG_OGs", "EC", "KEGG_ko", "KEGG_Module",
                 "KEGG_Reaction", "KEGG_rclass", "CAZy", "BiGG_Reaction", "PFAMs")

# Exact match check
matches_feature <- function(tokens_per_row, fid) {
  norm_fid <- toupper(fid)
  for (tokens in tokens_per_row) {
    if (length(tokens) && norm_fid %in% tokens) return(TRUE)
  }
  return(FALSE)
}

# Progress tracking
total_genomes <- length(needed_genomes)
pb <- txtProgressBar(min = 0, max = total_genomes, style = 3)
processed <- 0L
total_written_overall <- 0L
feature_write_count <- integer(length(feature_list))
names(feature_write_count) <- feature_list

log_msg("👉 Total genomes to process (once each): %d", total_genomes)

# ===================== Main loop ======================
for (gid in needed_genomes) {
  t0 <- Sys.time()
  processed <- processed + 1L
  setTxtProgressBar(pb, processed)

  sp <- genome_to_species[[gid]]
  if (is.null(sp) || !nzchar(sp)) next
  
  features_for_gid <- species_to_features[[sp]]
  if (!length(features_for_gid)) next
  
  anno_file  <- file.path(anno_dir,  paste0(gid, ".annotations.tsv"))
  fasta_file <- file.path(fasta_dir, paste0(gid, ".emapper.genepred.fasta"))
  if (!file.exists(anno_file) || !file.exists(fasta_file)) next
  
  # Read annotation
  header_probe <- readLines(anno_file, n = 200, warn = FALSE)
  skip_line <- which(grepl("^#query\\b", header_probe))[1]
  if (is.na(skip_line)) {
    anno_data <- fread(anno_file, sep = "\t", header = TRUE, skip = "#",
                       quote = "", fill = TRUE, showProgress = FALSE)
  } else {
    anno_data <- fread(anno_file, sep = "\t", header = TRUE, skip = skip_line - 1L,
                       quote = "", fill = TRUE, showProgress = FALSE)
  }
  if (!nrow(anno_data)) next
  
  for (col in search_cols) if (!col %in% names(anno_data)) anno_data[[col]] <- NA_character_
  query_col <- if ("#query" %in% names(anno_data)) "#query" else if ("query" %in% names(anno_data)) "query" else NA_character_
  if (is.na(query_col)) next
  
  # Feature → protein IDs mapping
  feature_to_pids <- list()
  for (fid in features_for_gid) {
    hit_idx <- which(apply(anno_data[, ..search_cols], 1, function(row) {
      row_tokens <- lapply(row, normalize_token)   # normalize each cell
      matches_feature(row_tokens, fid)
    }))
    if (length(hit_idx)) {
      feature_to_pids[[fid]] <- trimws(anno_data[[query_col]][hit_idx])
    }
  }
  
  feature_to_pids <- feature_to_pids[lengths(feature_to_pids) > 0]
  if (!length(feature_to_pids)) next
  
  # Stream FASTA
  needed_pids <- unique(unlist(feature_to_pids))
  needed_set <- new.env(parent = emptyenv())
  for (p in needed_pids) assign(p, TRUE, envir = needed_set)
  
  pid_to_block <- new.env(parent = emptyenv())
  con <- file(fasta_file, open = "r")
  current_pid <- NULL
  current_block <- NULL
  flush_block <- function() {
    if (!is.null(current_pid) && !is.null(current_block)) {
      assign(current_pid, current_block, envir = pid_to_block)
    }
  }
  while (length(ln <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(ln, ">")) {
      flush_block()
      main_token <- strsplit(sub("^>", "", ln), " ", fixed = TRUE)[[1]][1]
      pid <- strsplit(main_token, "___", fixed = TRUE)[[1]][1]
      if (exists(pid, envir = needed_set, inherits = FALSE)) {
        current_pid <- pid
        current_block <- ln
      } else {
        current_pid <- NULL
        current_block <- NULL
      }
    } else if (!is.null(current_block)) {
      current_block <- c(current_block, ln)
    }
  }
  flush_block()
  close(con)
  
  # Write output
  written_this_genome <- 0L
  for (fid in names(feature_to_pids)) {
    fasta_path <- file.path(out_dir, paste0(fid, ".fasta"))
    pids <- feature_to_pids[[fid]]
    out_lines <- character(0)
    for (pid in pids) {
      blk <- mget(pid, envir = pid_to_block, ifnotfound = list(NULL))[[1]]
      if (is.null(blk)) next
      header <- blk[1]
      pieces <- strsplit(sub("^>", "", header), " ", fixed = TRUE)[[1]]
      main_id <- pieces[1]
      meta_info <- if (length(pieces) > 1) paste(pieces[-1], collapse = " ") else ""
      main_parts <- strsplit(main_id, "___", fixed = TRUE)[[1]]
      if (length(main_parts) >= 3) {
        main_id_mod <- paste0(main_parts[1], "___", main_parts[2], "___", main_parts[3], "___", fid)
      } else {
        main_id_mod <- paste0(main_id, "___", fid)
      }
      new_header <- paste0(">", main_id_mod, if (nzchar(meta_info)) paste0(" ", meta_info) else "")
      out_lines <- c(out_lines, new_header, blk[-1])
      written_this_genome <- written_this_genome + 1L
    }
    if (length(out_lines)) {
      cat(out_lines, file = fasta_path, append = TRUE, sep = "\n")
      feature_write_count[fid] <- feature_write_count[fid] + length(out_lines) / 2
    }
  }
  
  total_written_overall <- total_written_overall + written_this_genome
  log_msg("✓ %s | features: %d | proteins written: %d | time: %.2fs",
          gid, length(feature_to_pids), written_this_genome,
          as.numeric(difftime(Sys.time(), t0, units = "secs")))
  
  rm(anno_data, feature_to_pids, pid_to_block, needed_set)
  gc(verbose = FALSE)
}

close(pb)

# ================= Final reports =================
log_msg("✅ All genomes processed. Total proteins written: %d", total_written_overall)
message(sprintf("Output FASTAs are in: %s", out_dir))

# Empty FASTA reporting
empty_fastas <- vapply(file.path(out_dir, list.files(out_dir)), function(f) file.size(f) == 0, logical(1))
if (any(empty_fastas)) {
  log_msg("⚠ Empty FASTA files: %d", sum(empty_fastas))
  log_msg("List: %s", paste(basename(names(empty_fastas[empty_fastas])), collapse = ", "))
}

# Top features by sequence count
top_features <- sort(feature_write_count, decreasing = TRUE)[1:10]
log_msg("Top 10 features by protein count:\n%s",
        paste(names(top_features), top_features, sep = ": ", collapse = "\n"))
