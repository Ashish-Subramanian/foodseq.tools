#' Quality Controls
#'
#' This function outputs quality control plots as relate to controls and possible contamination.
#'
#' @param marker The marker gene, either "trnL" or "12S"
#' @param proj The name of the project
#' @param ps The phyloseq object
#' @param control_species A character vector of positive control species names
#' @param out_dir The path to the output directory
#' @return A plot of control compositions, a plot of samples with positive controls detected, and a plate map
#' @export
qc_controls <- function(marker,
                        proj = project,
                        ps = NULL,
                        control_species = NULL,
                        out_dir = NULL) {

  if (is.null(ps)) {
    ps <- get(paste0("ps.", marker), envir = parent.frame())
  }

  if (is.null(control_species)) {
    control_species <- get(paste0("control.", marker), envir = parent.frame())
  }

  if (is.null(out_dir)) {
    out_dir <- get(paste0("qiime.dir.", marker), envir = parent.frame())
  }

  if (is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE)) ||
      ncol(phyloseq::sample_data(ps)) == 0) {
    message("[SKIP] phyloseq object has no sample metadata; nothing to run.")
    return(invisible(NULL))
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 1. Subset controls
  ps.subset <- phyloseq::subset_samples(ps, type != "sample")
  ps.controls <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps.subset) > 0, ps.subset)

  # 2. Build tax labels
  taxtab.controls <- phyloseq::tax_table(ps.controls)@.Data %>%
    data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(label = dplyr::coalesce(species, genus, family, order, phylum))
  phyloseq::tax_table(ps.controls) <- as.matrix(taxtab.controls)

  # 3. Plot controls composition
  p.controls <- ps.controls %>%
    phyloseq::psmelt() %>%
    ggplot2::ggplot(ggplot2::aes(x = Sample, y = Abundance, fill = label)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", "type")), scales = "free") +
    ggplot2::labs(
      x = 'Control',
      y = 'Number of reads',
      fill = 'ASV identity',
      title = paste0(proj, ": ", marker, " positive controls, negative controls, and blanks")
    )

  ggplot2::ggsave(filename = file.path(out_dir, "QC_controls.png"), plot = p.controls, width = 12, height = 6)

  # 4. Positive control species detected in non-positive-control samples
  df <- phyloseq::psmelt(ps)
  df$type <- as.character(df[["type"]])

  p.samples <- df %>%
    dplyr::filter(species %in% control_species,
           Abundance > 0,
           tolower(type) != "positive control") %>%
    ggplot2::ggplot(ggplot2::aes(x = Sample, y = Abundance, fill = species)) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", "type")), scales = "free") +
    ggplot2::labs(
      x = 'Sample',
      y = 'Number of reads',
      fill = 'Species',
      title = paste0(proj, ": Samples where ", marker, " positive controls were detected")
    )

  ggplot2::ggsave(filename = file.path(out_dir, "QC_Pos.Control-Detections.png"), plot = p.samples, width = 12, height = 6)

  # 5. Which samples show control species
  yesControl <- df %>%
    dplyr::filter(species %in% control_species) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(yesControl = any(Abundance > 0)) %>%
    dplyr::filter(yesControl) %>%
    dplyr::pull(Sample)

  sample_plot <- df %>%
    dplyr::filter(Sample %in% yesControl,
           tolower(type) %in% "sample") %>%
    ggplot2::ggplot(ggplot2::aes(x = Sample, y = Abundance, fill = species)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~Sample, scales = "free", nrow = 1) +
    ggplot2::labs(
      x = 'Sample',
      y = 'Number of reads',
      fill = 'Species',
      title = paste0(proj, ": Samples where ", marker, " positive controls were detected")
    ) +
    ggplot2::theme(legend.position = "bottom")

  # ggplot2::ggsave(filename = file.path(out_dir, "QC_Pos.Control-Detections_SamplesOnly.png"),
  #        plot = sample_plot, width = 14, height = 6)

  # 6. Plate map
  if (requireNamespace("ggplate", quietly = TRUE)) {
    contam_plot <- df %>%
      dplyr::filter(species %in% control_species) %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarise(reads = sum(as.numeric(Abundance)),
                yesControl = ifelse(any(Abundance > 0), "+", ""),
                type = dplyr::first(type)) %>%
      tidyr::separate(Sample, into = c("plate", "well"), sep = "-", remove = FALSE)

    for (Plate in unique(contam_plot$plate)) {
      dat <- dplyr::filter(contam_plot, plate == Plate)
      p.plate <- ggplate::plate_plot(
        data = dat,
        position = well,
        label = yesControl,
        value = type,
        plate_size = 96,
        plate_type = "round",
        title = paste0(proj, ": ", marker, " positive controls detected on barcode plate ", Plate),
        limits = c(1, 2),
        colour = c(
          "sample" = "#f0f0f0",
          "negative control" = "#ff8080",
          "positive control"    = "#80ff80",
          "blank"  = "#fff"
        )
      ) + ggplot2::labs(fill = "Type")
      ggplot2::ggsave(filename = file.path(out_dir, paste0("QC_Pos.ControlMap_Plate_", Plate, ".png")),
             plot = p.plate, width = 7, height = 7)
    }
  } else {
    message("ggplate not installed; skipping plate map.")
  }

  # Return results
  list(p.controls = p.controls,
       p.samples = p.samples,
       p.plate = p.plate)
}
