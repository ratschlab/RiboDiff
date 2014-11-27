library('babel')
test.rna <- read.table('/Users/zhongyi/Yi/Work/RiboDiff/exp/Sim/Sim.Rna.Rep3.G10K.Diff0K.cnt.babel.txt', header=TRUE)
test.rp <- read.table('/Users/zhongyi/Yi/Work/RiboDiff/exp/Sim/Sim.Ribo.Rep3.G10K3.Diff1K.Sh1.5.Sc0.5.cnt.babel.txt', header=TRUE)
test.group <- c('A', 'B', 'C', 'A', 'B', 'C')
options(mc.cores = 2)
set.seed(12345)
test.babel <- babel(test.rna, test.rp, group = test.group, nreps = 1e+06, min.rna = 10)
within.babel <- test.babel$within
combined.babel <- test.babel$combined
between.babel <- test.babel$between
write.table(as.data.frame(between.babel[[1]]), sep = '\t', file='/Users/zhongyi/Yi/Work/RiboDiff/exp/Sim/Sim.Merged.Rep3.G10K3.Diff1KRb.Sh1.5.Sc0.5.res.babel.txt')