setwd("D:/Git/project/missing-data/result/miss-result")
require(ggplot2)
require(ggpubr)

for (p in c(200, 400, 800)) {
  for (snr in c(0.5)) {
    file = paste("./p=", p, "--snr=", snr, "/miss.RData", sep = "")
    plot.name <- paste("./p=", p, "--snr=", snr*100, ".eps", sep = "")
    load(file)
    table$miss <- as.character(table$miss)
    table$ErXC <- table$ErXC / 1000

    # ErC
    p.c <- ggplot(table) +
      aes(x = miss, y = ErC, fill = method) +
      geom_boxplot(notch = T) +
      xlab("missing rate") +
      ylab(expression("Er(C)" %*% 10 ^ 3)) + 
      scale_fill_manual(breaks = c("peer(l1)", "peer(l0)"), 
                        labels = c(expression("PEER-l"[1]), expression("PEER-l"[0])), 
                        values = c("gray80", "gray40"))

    # ErXC
    p.xc <- ggplot(table) +
      aes(x = miss, y = ErXC, fill = method) +
      geom_boxplot(notch = T) +
      xlab("missing rate") +
      ylab("Er(XC)")  + 
      scale_fill_manual(breaks = c("peer(l1)", "peer(l0)"), 
                        labels = c(expression("PEER-l"[1]), expression("PEER-l"[0])), 
                        values = c("gray80", "gray40"))
    # FPR
    p.fpr <- ggplot(table) +
      aes(x = miss, y = FPR, fill = method) +
      geom_boxplot(notch = T) +
      xlab("missing rate") +
      ylab("FPR(%)")  + 
      scale_fill_manual(breaks = c("peer(l1)", "peer(l0)"), 
                        labels = c(expression("PEER-l"[1]), expression("PEER-l"[0])), 
                        values = c("gray80", "gray40"))
    # FNR
    p.fnr <- ggplot(table) +
      aes(x = miss, y = FNR, fill = method) +
      geom_boxplot(notch = F) +
      xlab("missing rate") +
      ylab("FNR(%)")  + 
      scale_fill_manual(breaks = c("peer(l1)", "peer(l0)"), 
                        labels = c(expression("PEER-l"[1]), expression("PEER-l"[0])), 
                        values = c("gray80", "gray40"))
    
    ratio <- 978 / 289 # width / height
    width <- 10.1875
    
    p <- ggarrange(p.c, p.xc, p.fpr, p.fnr, ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom")
    
    setEPS()
    postscript(file = plot.name, width = width, height = width / ratio)
    plot(p)
    dev.off()
  }
}



print(p)
