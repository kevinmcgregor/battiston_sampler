# Testing Battiston sampler
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

# Reading in parameters for simulation
args <- commandArgs(trailingOnly = TRUE)
sim_number <- as.numeric(args[1])
arg_file <- args[2]

f.args <- read.table(arg_file, header = FALSE, stringsAsFactors = FALSE)$V1
N <- as.numeric(f.args[1])
J <- as.numeric(f.args[2])
conc.top <- as.numeric(f.args[3])
dsct.top <- as.numeric(f.args[4])
conc.local <- as.numeric(f.args[5])
dsct.local <- as.numeric(f.args[6])
out.dir <- f.args[7]

# Sampler parameters
n.iter <- 2000
n.burn <- 500

# Sampling from HPY
s.hpy <- sample_hpy(J, rep(N, J), conc.top, dsct.top, 
                    rep(conc.local, J), rep(dsct.local, J))
Y <- s.hpy$abund

# Running parallel chains
library(parallel)
n.chain <- 4
Y.list <- vector("list", n.chain)
for (i in 1:n.chain) Y.list[[i]] <- Y

s.p <- mclapply(Y.list, b_sampler, p.shape = 1, p.scale = 1,
                n.iter=n.iter, n.burn=n.burn, mc.cores=n.chain)

save.image(paste0(out.dir, "/b_sampler_sim", sim_number, ".RData"))
