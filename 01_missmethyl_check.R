# 01_missmethyl_check.R - Gene Ontology test of the age-related methylation sites with the Bioconductor package missMethyl.
# What it does: takes the sites the elastic net chose (from the methylation project) and asks which
# Gene Ontology terms are over-represented among their genes. The function gometh corrects for the fact that
# some genes have many more sites on the array than others (it uses the annotation of the EPIC array).
# This is the reference result that the Python enrichment in the notebook is compared with.
#
# Run from the repository folder:  Rscript 01_missmethyl_check.R
# Expected output: one line per list of sites, the files results/missmethyl_*.csv and data/site_to_gene_table.csv.
# A list whose result file already exists is skipped. To compute everything again, set RERUN_R=1 before running.

# ---------- Settings ----------
input_folder <- file.path("input", "methylation_repo")
results_folder <- "results"
seeds <- c(0, 1, 42)
ontology_to_keep <- "BP"               # BP = biological process, the same as the Python notebook
# GPL21145 is the Infinium MethylationEPIC array (version 1), see bridge_information.json
array_type <- "EPIC"

dir.create(results_folder, showWarnings = FALSE)

# ---------- Step 1: install the packages once ----------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (package_name in c("missMethyl", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        BiocManager::install(package_name, ask = FALSE, update = FALSE)
    }
}
library(missMethyl)
cat("missMethyl version:", as.character(packageVersion("missMethyl")), "\n")
cat("R version:", R.version.string, "\n")

# ---------- Step 2: read the inputs ----------
bridge_text <- paste(readLines(file.path(input_folder, "bridge_information.json"), warn = FALSE), collapse = " ")
if (!grepl("GPL21145", bridge_text)) {
    stop("bridge_information.json does not name the platform GPL21145. Check that the array type is still correct.")
}

all_sites_in_study <- readLines(file.path(input_folder, "all_sites_in_study.txt"), warn = FALSE)
all_sites_in_study <- all_sites_in_study[nchar(all_sites_in_study) > 0]
# Only sites whose name starts with cg have a gene annotation. The others (names starting with ch.) are removed.
background_sites <- all_sites_in_study[startsWith(all_sites_in_study, "cg")]
cat("Sites in the study:", length(all_sites_in_study), "| background sites (cg only):", length(background_sites), "\n")

sites_with_weight_by_seed <- list()
for (seed in seeds) {
    table_for_seed <- read.csv(file.path(input_folder, paste0("elastic_net_sites_seed_", seed, ".csv")), stringsAsFactors = FALSE)
    chosen <- table_for_seed$site[as.logical(table_for_seed$weight_is_not_zero)]
    sites_with_weight_by_seed[[as.character(seed)]] <- chosen[startsWith(chosen, "cg")]
    cat("Seed", seed, "| weights other than zero:", sum(as.logical(table_for_seed$weight_is_not_zero)),
        "| after removing ch. sites:", length(sites_with_weight_by_seed[[as.character(seed)]]), "\n")
}
sites_in_all_seeds <- Reduce(intersect, sites_with_weight_by_seed)
cat("Sites with a weight other than zero in all seeds:", length(sites_in_all_seeds), "\n")

lists_to_test <- list(
    seed_0 = sites_with_weight_by_seed[["0"]],
    seed_1 = sites_with_weight_by_seed[["1"]],
    seed_42 = sites_with_weight_by_seed[["42"]],
    all_seeds = sites_in_all_seeds
)

# ---------- Step 3: the test ----------
# The annotation of the array is loaded once and reused for every list. Without it, gometh builds it again for each list.
annotation_table <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
cat("Annotation loaded:", nrow(annotation_table), "rows at", format(Sys.time()), "\n")

# Table for the Python notebook: the gene names of every site on the array.
dir.create("data", showWarnings = FALSE)
site_to_gene_table <- data.frame(
    site = rownames(annotation_table),
    chromosome = annotation_table$chr,
    position = annotation_table$pos,
    gene_names = annotation_table$UCSC_RefGene_Name,
    gene_groups = annotation_table$UCSC_RefGene_Group,
    stringsAsFactors = FALSE
)
write.csv(site_to_gene_table, file.path("data", "site_to_gene_table.csv"), row.names = FALSE)
cat("Site-to-gene table written:", nrow(site_to_gene_table), "rows\n")

for (list_name in names(lists_to_test)) {
    significant_sites <- lists_to_test[[list_name]]
    if (length(significant_sites) < 5) {
        cat(list_name, ": fewer than 5 sites, skipped.\n")
        next
    }
    output_path <- file.path(results_folder, paste0("missmethyl_", list_name, ".csv"))
    if (file.exists(output_path) && Sys.getenv("RERUN_R") != "1") {
        cat(list_name, ": result file already exists, skipped.\n")
        next
    }

    cat("Starting", list_name, "at", format(Sys.time()), "\n")
    # all.cpg = the fair comparison group (every site of the study). prior.prob = TRUE corrects for the
    # unequal number of sites per gene. sig.genes is off: it only lists the genes behind each term, and it
    # made the call far slower (tested: 156 seconds without it, no result after 35 minutes with it).
    gene_ontology_result <- gometh(
        sig.cpg = significant_sites,
        all.cpg = background_sites,
        collection = "GO",
        array.type = array_type,
        anno = annotation_table,
        prior.prob = TRUE,
        plot.bias = FALSE,
        sig.genes = FALSE
    )
    gene_ontology_result <- gene_ontology_result[gene_ontology_result$ONTOLOGY == ontology_to_keep, ]
    gene_ontology_result <- gene_ontology_result[order(gene_ontology_result$P.DE), ]

    # Genes behind the 20 lowest-p terms: the genes of your sites that belong to each term.
    # (A plain lookup with the array annotation and org.Hs.eg.db; it does not change any p-value.)
    site_gene_names <- annotation_table[significant_sites, "UCSC_RefGene_Name"]
    genes_of_list <- unique(unlist(strsplit(site_gene_names[site_gene_names != ""], ";")))
    top_terms <- head(rownames(gene_ontology_result), 20)
    term_gene_table <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = top_terms,
                                             keytype = "GOALL", columns = "SYMBOL")
    gene_ontology_result$genes_in_list <- ""
    gene_ontology_result$genes_listed_count <- NA_integer_
    for (term_id in top_terms) {
        genes_here <- sort(intersect(genes_of_list, term_gene_table$SYMBOL[term_gene_table$GOALL == term_id]))
        gene_ontology_result[term_id, "genes_in_list"] <- paste(genes_here, collapse = ";")
        gene_ontology_result[term_id, "genes_listed_count"] <- length(genes_here)
    }

    write.csv(gene_ontology_result, output_path)
    cat(list_name, "| sites tested:", length(significant_sites), "| terms with FDR below 0.05:",
        sum(gene_ontology_result$FDR < 0.05), "| written to", output_path, "at", format(Sys.time()), "\n")
}
cat("missMethyl check finished.\n")