# Plotting bias for Battiston sampler
# (Simulations run on Graham cluster)
library(ggplot2)

plot.dir <- "~/research/pitman_yor/battiston_sampler/plots/sim/"

load("/mnt/GREENWOOD_SCRATCH/kevin.mcgregor/pitman_yor/battiston_sampler/sim/results/sim_results.RData")

n.sim <- NCOL(alpha.est)
n.rep <- NROW(alpha.est)
J <- dim(theta)[3]

b.gamma <- (gamma.est-gamma)/gamma
b.alpha <- (alpha.est-alpha)/alpha
b.theta <- (theta.est-theta)/theta
b.sigma <- (sigma.est-sigma)/sigma
b.n.tab <- (n.tab.est-n.tab)/n.tab

# Plotting results
# Gamma (top level concentration)
top.df.gamma <- data.frame(val=c(b.gamma),
                     sim=factor(rep(1:5, each=n.rep)))

plot.title.gamma = "Top-level concentration"
gamma.plot = ggplot(top.df.gamma, aes(sim, val, fill=sim)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.gamma) +
  ylim(0,5) +
  xlab("Simulation") + ylab("(Normalized) Bias") + cowplot::panel_border()
ggsave(paste0(plot.dir, "/gamma_bias.pdf"), gamma.plot, "pdf",
       width=12, height=8)

# Alpha (top level discount)
top.df.alpha <- data.frame(val=c(b.alpha),
                     sim=factor(rep(1:5, each=n.rep)))

plot.title.alpha = "Top-level discount"
alpha.plot = ggplot(top.df.alpha, aes(sim, val, fill=sim)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.alpha) +
  #ylim(0,5) +
  xlab("Simulation") + ylab("(Normalized) Bias") + cowplot::panel_border()
ggsave(paste0(plot.dir, "/alpha_bias.pdf"), alpha.plot, "pdf",
       width=12, height=8)

# Theta (Local level concentrations)
bt.tmp <- aperm(b.theta, c(1,3,2))
local.df.theta <- data.frame(val=c(bt.tmp),
                       sim=factor(rep(1:5, each=n.rep*J)),
                       pop=factor(rep(rep(1:10, each=n.rep), n.sim)))

plot.title.local.theta = "Local-level concentration"
local.plot.theta = ggplot(local.df.theta, aes(sim, val, fill=pop)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.local.theta) +
  xlab("Simulation") + ylab("(Normalized) Bias") + cowplot::panel_border()
ggsave(paste0(plot.dir, "/theta_bias.pdf"), local.plot.theta, "pdf",
       width=12, height=8)

# Sigma (Local level discounts)
bs.tmp <- aperm(b.sigma, c(1,3,2))
local.df.sigma <- data.frame(val=c(bs.tmp),
                             sim=factor(rep(1:5, each=n.rep*J)),
                             pop=factor(rep(rep(1:10, each=n.rep), n.sim)))

plot.title.local.sigma = "Local-level discount"
local.plot.sigma = ggplot(local.df.sigma, aes(sim, val, fill=pop)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.local.sigma) +
  xlab("Simulation") + ylab("(Normalized) Bias") + cowplot::panel_border()
ggsave(paste0(plot.dir, "/sigma_bias.pdf"), local.plot.sigma, "pdf",
       width=12, height=8)

# Number of tables
bnt.tmp <- aperm(b.n.tab, c(1,3,2))
tab.df <- data.frame(val=c(bnt.tmp),
                     sim=factor(rep(1:5, each=n.rep*J)),
                     pop=factor(rep(rep(1:10, each=n.rep), n.sim)))

plot.title.table = "Local ancestors"
table.plot = ggplot(tab.df, aes(sim, val, fill=pop)) + geom_boxplot() + #facet_wrap( ~ z*n.otu, ncol=2) + 
  theme_bw() +
  theme(legend.position = "bottom",text = element_text(size=20), axis.text.x = element_text(size=16,angle=0,vjust=1,hjust=0.5)) + 
  labs(title=plot.title.table) +
  xlab("Population") + ylab("(Normalized) Bias") + cowplot::panel_border()
ggsave(paste0(plot.dir, "/table_bias.pdf"), table.plot, "pdf",
       width=12, height=8)

