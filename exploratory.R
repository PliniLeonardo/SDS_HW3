setClass("num.with.commas")
setAs("character", "num.with.commas", 
      function(from) as.numeric(gsub(',', '.', gsub("\\.", "",from))) )

base_dir = './data/cell_signaling/original'

df = read.csv('data/cell_signaling/cd3cd28+aktinhib.csv', colClasses = 'num.with.commas')

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

plot_normalized_distros =function(df){
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
    hist(normalized, add=T, col = rgb(0,0,1,alpha=0.3), breaks=30)
  }
  par(mfrow=c(1,1))  
}

store_normalized_distros = function(base_dir){
  for(fname in list.files(base_dir, recursive = T, full.names = F)){
    csv_path = paste(base_dir, fname, sep = '/')
    config = substr(fname, start=0, stop=(nchar(fname)-4))
    png_path = paste(base_dir, '/', config, '.png', sep = '')
    png(png_path, width = 1024, height = 700)
    
    df= read.csv(csv_path, colClasses = 'num.with.commas')
    plot_normalized_distros(df)
    dev.off()
  }
}

normalize_dataframe = function(df){
  out = list()
  for (n in names(df)){
    data = df[[n]]
    tmp=setNames(list(as.numeric(normalize(data))), n)
    out=append(out, tmp)
  }
  return(out)
}

store_normalized_dataframes = function(base_dir){
  for(fname in list.files(base_dir, recursive = T, full.names = F)){
    csv_path = paste(base_dir, fname, sep = '/')
    config = substr(fname, start=0, stop=(nchar(fname)-4))
    out_path = paste(base_dir, '/', config, '_normalized.csv', sep = '')
    
    df= read.csv(csv_path, colClasses = 'num.with.commas')
    out = data.frame(normalize_dataframe(df))
    write.csv(out, file = out_path, row.names = F)
  }
}

plot_df = function(df){
  par(mfrow=c(4,3))
  for (n in names(df)){
    data = df[[n]]
    hist(data, main=n, breaks = 30, col = rgb(0,1,0))
  }
  par(mfrow=c(1,1))
}

