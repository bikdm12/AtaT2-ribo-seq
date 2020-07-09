library('compositions')
library('ICSNP')
library('energy')

ctrl <- read.csv('Control_paired_raw.csv', header = F)
atat2 <- read.csv('TacT2_paired_raw.csv', header = F)
ctrl_ilr <- read.csv('Control_paired_ilr.csv', header = F)
atat2_ilr <- read.csv('TacT2_paired_ilr.csv', header = F)

acompNormalGOF.test(acomp(ctrl), R = 1000)
acompNormalGOF.test(acomp(atat2), R = 1000)
# the next two lines are equivalent to the previous two
mvnorm.etest(ctrl_ilr, R = 1000)
mvnorm.etest(atat2_ilr, R = 1000)

d <- ctrl_ilr - atat2_ilr
o <- rep(0, times = 9)
HotellingsT2(d, mu = o)
