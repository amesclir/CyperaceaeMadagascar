library(diversitree)

phy <- read.tree("tree.tree")
mydata2 <- read.csv ("mydata.csv")
states <- as.character(mydata2$Species)
names(states) <- mydata2$Area
setdiff(phy$tip.label, states)
setdiff(states,phy$tip.label)
setdiff(phy$tip.label, states)
setdiff(states,phy$tip.label)
phy$tip.state <- names(states)[match(phy$tip.label,states)]
names(phy$tip.state) <- states[match(phy$tip.label,states)] 

p <- starting.point.geosse(phy)
p

phy$tip.state[phy$tip.state == "A"] <- 1
phy$tip.state[phy$tip.state == "B"] <- 2
phy$tip.state[phy$tip.state == "AB"] <- 0
phy$tip.state <- as.numeric(phy$tip.state)
names(phy$tip.state) <- states[match(phy$tip.label,states)] 

lik1 <- make.geosse(phy, phy$tip.state, sampling.f = c(0.5,0.5,0.25))
lik2 <- constrain(lik1, sAB ~ 0)
lik3 <- constrain(lik1, sA ~ sB, xA ~ xB)

ml1 <- find.mle(lik1, p)
p <- coef(ml1)
ml2 <- find.mle(lik2, p[argnames(lik2)])
ml3 <- find.mle(lik3, p[argnames(lik3)])

round(rbind(full = coef(ml1),no.sAB = coef(ml2, TRUE),eq.div = coef(ml3, TRUE)), 3)

anova(ml1, no.sAB = ml2, eq.div = ml3)



p <- coef(ml2)
prior <- make.prior.exponential(1/2)

set.seed(1)
tmp <- mcmc(lik2, p, nsteps=1000, prior=prior, w=1, print.every=0)
w <- diff(sapply(tmp[2:7], quantile, c(0.025, 0.975)))

mcmc2 <- mcmc(lik2, p, nsteps=50000, prior=prior, w=w)

save.image(file="samples2.RData")
