# summarize simulation results

S <- 77
prefix <- "~/Dropbox/confounding/output"

dwpr <- matrix(NA, 0, 4)
ortho <- non.ortho <- matrix(NA, 0, 5)

for (s in 1:S) {
  file <- paste(paste(prefix, "dwpr", sep='/'), s, 'out', sep='.')
  data <- scan(file, what=character())
  dwpr <- rbind(dwpr, as.numeric(data[5:8]))
  
  file <- paste(paste(prefix, "ortho", sep='/'), s, 'out', sep='.')
  data <- scan(file, what=character())
  ortho <- rbind(ortho, as.numeric(data[1:5]))
  
  file <- paste(paste(prefix, "non.ortho", sep='/'), s, 'out', sep='.')
  data <- scan(file, what=character())
  non.ortho <- rbind(non.ortho, as.numeric(data[1:5]))
}