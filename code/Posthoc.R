# Run from the repository root after the four main experiments.
results <- lapply(c("AR2", "NL", "elec", "cafe"), function(dataset) {
  message("Running post-hoc comparison: ", dataset)
  e <- new.env(parent = globalenv())
  sys.source(paste0("code/Sim_", dataset, "_param.R"), envir = e)
  do.call(rbind, lapply(c("macp", "acmcp"), function(method) {
    cov <- e$cov_list[[method]]
    wid <- e$wid_list[[method]]
    data.frame(
      dataset = dataset, method = if (method == "macp") "MACP" else "AcMCP",
      horizon = seq_len(e$horizon),
      gamma = if (method == "macp") e$gamma else NA_real_,
      coverage = 100 * as.vector(cov$mean),
      coverage_min = 100 * apply(cov$rollmean, 2, min, na.rm = TRUE),
      coverage_max = 100 * apply(cov$rollmean, 2, max, na.rm = TRUE),
      width_mean = as.vector(wid$mean), width_median = as.vector(wid$median),
      row.names = NULL
    )
  }))
})
saveRDS(do.call(rbind, results), "result/posthoc_results.rds")
