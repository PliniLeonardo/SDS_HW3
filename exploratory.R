b2 = read.csv('data/cell_signaling/b2camp.csv')
names(b2)
praf = as.numeric(b2$praf)
pmek = as.numeric(b2$pmek)
hist(praf)
hist(pmek)
boxplot(praf)
qqplot(praf, pmek)

par(mfrow=c(4,3))
for (n in names(b2)){
  hist(as.numeric(b2[[n]]), main=n)
}

df = read.csv('data/cell_signaling/cd3cd28+aktinhib.csv', colClasses = 'num.with.commas')
par(mfrow=c(4,3))
for (n in names(df)){
  hist(as.numeric(df[[n]]), main=n, breaks = 30)
}
plot(1:5, 5:9)


setClass("num.with.commas")
setAs("character", "num.with.commas", 
      function(from) as.numeric(gsub(',', '.', gsub("\\.", "",from))) )

base_dir = './data/cell_signaling'
for(fname in list.files(base_dir, recursive = T, full.names = F)){
  csv_path = paste(base_dir, fname, sep = '/')
  config = substr(fname, start=0, stop=(nchar(fname)-4))
  png_path = paste(base_dir, '/', config, '.png', sep = '')
  png(png_path, width = 1024, height = 700)
  par(mfrow=c(4,3))
  df= read.csv(csv_path, colClasses = 'num.with.commas')
  for (n in names(df)){
    hist(as.numeric(df[[n]]), main=n, breaks = 30)
  }
  dev.off()
}

