# Main function
ABarrayTest <- function(rawData, Ref, Sample, outDir, B = 1e4, today=format(Sys.Date(), format="%Y-%m-%d")){
  idx <- grep(sprintf("%s$|%s$", Ref, Sample), colnames(rawData))
  cat("Running tests...\n")
  ratio <- computeLRatios(rawData, Ref, Sample, B)
  cat("Saving plots...\n")
  png(file.path(outDir, sprintf("%s_%s_Signal_%s.png", Sample, Ref, today)), width = 1200, height = 600)
  plotSignal(rawData, idx)
  dev.off()
  png(file.path(outDir, sprintf("%s_%s_Profile_%s.png", Sample, Ref, today)), width = 1200, height = 600)
  plotLRatios(ratio, main = paste(Sample, 'Vs.', Ref), useP = "adjp",thresh=.05)
  dev.off()
  png(file.path(outDir, sprintf("%s_%s_vPlot_%s.png", Sample, Ref, today)), width = 700, height = 600)
  vPlot(ratio, useP = "adjp", thresh=.05)
  dev.off()
  cat("Saving results...\n")
  write.table(ratio, file.path(outDir, sprintf('Results_%s_vs_%s_%s.xls', Sample, Ref, today)),
              sep = '\t', row.names = FALSE)
  cat("Done.\n\n")
  return(ratio)
}
