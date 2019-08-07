# Plotting traceplots from run_b_sampler.R

load("/mnt/GREENWOOD_SCRATCH/kevin.mcgregor/pitman_yor/battiston_sampler/b_sampler_test.RData")
plot.dir <- "/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/battiston_sampler/plots"

# Checking traceplots for parallel chains ####
col <- rainbow(n.chain)

# tables
pdf(paste0(plot.dir, "/table_traceplots_long.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  tmp <- unlist(lapply(s.p, function(x) {x$n.tab[,j]}))
  y.rge <- range(tmp, s.hpy$t.tab[j])
  plot(s.p[[1]]$n.tab[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$n.tab[,j], ylim=y.rge, col=col[[i]])
    abline(h=s.hpy$t.tab[j], col="black", lty=2, lwd=2)
  }
}
dev.off()

# top-level concentration param
pdf(paste0(plot.dir, "/gamma_traceplot_long.pdf"), height=8, width=12)
par(mfrow=c(1,1), mar=rep(2,4))
tmp <- unlist(lapply(s.p, function(x) {x$gamma}))
y.rge <- range(tmp, conc.top)
plot(s.p[[1]]$gamma, type="l", ylim=y.rge, col=col[1])
for (i in 1:n.chain) {
  lines(s.p[[i]]$gamma, ylim=y.rge, col=col[[i]])
}
abline(h=conc.top, col="black", lty=2, lwd=2)
dev.off()

# top-level discount param
pdf(paste0(plot.dir, "/alpha_traceplot_long.pdf"), height=8, width=12)
par(mfrow=c(1,1), mar=rep(2,4))
tmp <- unlist(lapply(s.p, function(x) {x$alpha}))
y.rge <- range(tmp, dsct.top)
plot(s.p[[1]]$alpha, type="l", ylim=y.rge, col=col[1])
for (i in 1:n.chain) {
  lines(s.p[[i]]$alpha, ylim=y.rge, col=col[[i]])
}
abline(h=dsct.top, col="black", lty=2, lwd=2)
dev.off()

# local-level concentration params
pdf(paste0(plot.dir, "/theta_traceplots_long.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  tmp <- unlist(lapply(s.p, function(x) {x$theta[,j]}))
  y.rge <- range(tmp, conc.local[j])
  plot(s.p[[1]]$theta[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$theta[,j], ylim=y.rge, col=col[[i]])
    abline(h=conc.local[j], col="black", lty=2, lwd=2)
  }
}
dev.off()

# local-level discount params
pdf(paste0(plot.dir, "/sigma_traceplots_long.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  tmp <- unlist(lapply(s.p, function(x) {x$sigma[,j]}))
  y.rge <- range(tmp, dsct.local[j])
  plot(s.p[[1]]$sigma[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$sigma[,j], ylim=y.rge, col=col[[i]])
    abline(h=dsct.local[j], col="black", lty=2, lwd=2)
  }
}
dev.off()
