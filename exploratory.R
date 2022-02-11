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

df = read.csv('data/cell_signaling/cd3cd28+aktinhib.csv')
par(mfrow=c(4,3))
for (n in names(df)){
  hist(as.numeric(df[[n]]), main=n)
}
#list.files('./data/cell_signaling', recursive = T, full.names = T)