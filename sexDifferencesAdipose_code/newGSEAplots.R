
GSEAplots = function (input.ds.name = "", input.cls.name = "", gene.set.input = "", 
          doc.string = "", nperm = 1000, fdr.q.val.threshold = 0.25, 
          bar_percent = 0.1, gap_percent = 0.1, under_percent = 0.1, 
          upper_percent = 0.1, color_line = "black", color_tick = "black", 
          abs.val = F, gs.size.threshold.max = 1000) 
{
  wd_new = getwd()
  if (file.exists(paste0(wd_new, "/", doc.string)) == FALSE) {
    dir.create(doc.string)
  }
  
  results_new = GSEA(input.ds = input.ds.name, input.cls = input.cls.name, 
                     gs.db = gene.set.input, output.directory = paste0(wd_new, 
                                                                       "/", doc.string, "/"), output.directory2 = paste0(wd_new, 
                                                                                                                         "/"), doc.string = doc.string, non.interactive.run = T, 
                     reshuffling.type = "sample.labels", nperm = nperm, weighted.score.type = 1, 
                     nom.p.val.threshold = -1, fwer.p.val.threshold = -1, 
                     fdr.q.val.threshold = 0.25, topgs = 20, adjust.FDR.q.val = F, 
                     gs.size.threshold.min = 15, gs.size.threshold.max = gs.size.threshold.max, 
                     reverse.sign = F, preproc.type = 0, random.seed = 3338, 
                     perm.type = 0, fraction = 1, replace = F, save.intermediate.results = F, 
                     OLD.GSEA = F, use.fast.enrichment.routine = T, abs.val = abs.val)
  
  
  if (length(which(sapply(results_new$out5, is.null))) == 0) {
    ES.tags.files <- results_new$out5
    ES.data.files <- results_new$out6
    ES.report.files <- results_new$out7
    gene.set.numbers <- results_new$out4
  } else {
    ES.tags.files <- results_new$out5[-which(sapply(results_new$out5, 
                                                    is.null))]
    ES.data.files <- results_new$out6[-which(sapply(results_new$out6, 
                                                    is.null))]
    ES.report.files <- results_new$out7[-which(sapply(results_new$out7, 
                                                      is.null))]
    gene.set.numbers <- results_new$out4[-which(sapply(results_new$out5, 
                                                       is.null))]
  }
  
  gene.set.reference.matrix <- results_new$gene.set.reference.matrix
  gene.set.leading <- rep(list("null"), length(gene.set.numbers))
  ES <- rep(list("null"), length(gene.set.numbers))
  enrichind <- rep(list("null"), length(gene.set.numbers))
  report1 <- results_new$report1
  report2 <- results_new$report2
  
  # if (regexpr(pattern = "HALLMARK_", gene.set.numbers[1]) == 
  #     -1) {
  # } else {
  #   for (i in 1:length(gene.set.numbers)) {
  #     g <- strsplit(gene.set.numbers[[i]], split = "HALLMARK_")
  #     h <- g[[1]][2]
  #     gene.set.numbers[[i]] <- h
  #   }
  # }
  
  plots <- vector(mode = "list", length = length(ES.data.files))
  for (i in 1:length(ES.tags.files)) {
    dat1_name = paste(wd_new, "/", doc.string, "/", doc.string, 
                      ".", ES.data.files[[i]], sep = "", collapse = NULL)
    dat2_name = paste(wd_new, "/", doc.string, "/", doc.string, 
                      ".", ES.tags.files[[i]], sep = "", collapse = NULL)
    report_name = paste(wd_new, "/", doc.string, "/", doc.string, 
                        ".", ES.report.files[[i]], sep = "", collapse = NULL)
    dat1 = read.table(dat1_name, header = T, sep = "\t")
    dat2 = read.table(dat2_name, header = T, sep = "\t")
    report = read.table(report_name, sep = "\t")
    datcomb <- cbind(dat1, dat2)
    datcomb = datcomb[, -5]
    ES[[i]] = datcomb
    enrich_ind = which(dat2$EStag == 1)
    height = max(dat1$RES) - min(dat1$RES)
    bar_height = bar_percent * height
    gap_height = gap_percent * height
    under_height = under_percent * height
    upper_height = upper_percent * height
    window_height = height + bar_height + gap_height + under_height + 
      upper_height
    y_lower = min(dat1$RES) - gap_height - bar_height
    window_low = y_lower - under_height
    window_high = max(dat1$RES) + upper_height
    d = data.frame(x = enrich_ind, y = matrix(y_lower, length(enrich_ind), 
                                              1), vx = matrix(0, length(enrich_ind), 1), vy = matrix(bar_height, 
                                                                                                     length(enrich_ind), 1))
    p <- ggplot(datcomb, aes(index, RES)) + geom_line(col = color_line)
    p <- p + geom_segment(data = d, mapping = aes(x = x, 
                                                  y = y, xend = x + vx, yend = y + vy), col = color_tick)
    p <- p + theme_classic()
    p <- p + ylim(window_low, window_high)
    p <- p + ggtitle(gene.set.numbers[[i]])
    p <- p + theme(plot.title = element_text(size = 12))
    plots[[i]] <- p
    leading_index <- which.max(abs(datcomb$RES))
    if (datcomb$RES[leading_index] > 0) {
      leading.edge.names <- report$V2[which(report$V5 < 
                                              leading_index)]
      leading.edge.RES <- report$V7[which(report$V5 < leading_index)]
    }
    else {
      leading.edge.names <- report$V2[which(report$V5 > 
                                              leading_index)]
      leading.edge.RES <- report$V7[which(report$V5 > leading_index)]
    }
    leading.edge.set <- data.frame(leading.edge.names, leading.edge.RES)
    gene.set.leading[[i]] <- leading.edge.names
  }
  names(gene.set.leading) <- gene.set.numbers
  names(ES) <- gene.set.numbers
  gene.set.leading[] <- lapply(gene.set.leading, as.character)
  
  return(list(plots = plots, gene.set.reference.matrix = gene.set.reference.matrix, 
              gene.set.leading = gene.set.leading, report1 = report1, 
              report2 = report2, ES = ES))
}

  