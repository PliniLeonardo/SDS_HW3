b2 = read.csv('data/cell_signaling/b2camp.csv', colClasses = 'num.with.commas')
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
  data = df[[n]]
  scaled = scale(data)
  normalized = normalize(data)
  
  data_hist = hist(scaled, breaks=30,plot=F)
  norm_data_hist = hist(normalized, breaks=30,plot=F)
  
  hist(scaled, main=n, breaks = 30, col = rgb(0,1,0, alpha = 0.6),
       xlim=c(min(min(data_hist$breaks),min(norm_data_hist$breaks)), max(max(data_hist$breaks),max(norm_data_hist$breaks))),
       ylim=c(0, max(max(data_hist$counts),max(norm_data_hist$counts)))
  )
  hist(normalized, add=T, col = rgb(0,0,1,alpha=0.3))
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

boxcox.transform= function(data){
  bx = boxcox(data~1, plotit = F)
  lambda = bx$x[which(bx$y==max(bx$y))[1]]
  if(lambda != 0)
    return((data^lambda-1)/lambda)
  else
    return(log(data))
}
normalize = function(data)
  return(scale(boxcox.transform(data)))
