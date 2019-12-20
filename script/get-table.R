#setwd("./project/missing-data/result/table-result/")

require(xlsx)

my.mean <- function(x){
  return(mean(x, na.rm = TRUE))
}

my.sd <- function(x){
  return(sd(x, na.rm = TRUE))
}

for (snr in c(0.25, 0.5, 1)) {
  for (p in c(100, 200, 400)) {
    for (r in c(3, 6)) {
      file = paste("./snr=", snr, "/p=", p, "--r=", r, ".Rdata", sep = "")
      load(file) 
      table$time <- as.double(table$time)
      table$ErXC <- table$ErXC / 1000
      method <- c('RRR', 'SRRR', 'SEED', 'PEER(l1)', 'PEER(l0)')
      performance <- colnames(table)[1:5]
      ave <- std <- data.frame()
      for (i in c(1:length(method))){
        temp <- table[table$method == method[i], 1:5]
        colnames(temp) <- NULL
        ave <- rbind(ave, apply(temp, 2, my.mean))
        std <- rbind(std, apply(temp, 2, my.sd))
      }
      rownames(ave) <- method
      rownames(std) <- method 
      colnames(ave) <- performance 
      colnames(std) <- performance
      
      # wirte table
      filename <- paste("p=", p, "--r=", r, ".xlsx", sep = "")
      path <- paste("./snr=", snr, "/", filename, sep = "")
      write.xlsx(ave, file = path, sheetName = "ave")
      write.xlsx(std,file = path,sheetName = "std",append = TRUE)
    }
  }
}
