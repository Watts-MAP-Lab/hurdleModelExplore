## May 29th 2025
sim.param.nitems <- c(7,14,21,28)
sim.param.ncates <- c(3,5)
sim.param.discri <- c(.6,2)
sim.param.2pl.spread <- c(3)
sim.param.sample <- c(15000)
sim.param.faccor <- c(.4,.8)
sim.param.difgrmF <- c(-3,-1)
sim.param.difgrmC <- c(1)
sim.param.dif2pl <- c(-2,1)
sim.param.discri2 <- c(.6,2)
sim.iter <- 1:50
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.discri2,sim.param.dif2pl, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCat", "discrim2pl", "diffSpread", "n", "facCor", "difGrmF","difGrmC","grmDiscrim","dif2PL","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
## NOw add the simulation permutation to these
rep.iter <- table(all.sim.vals$iter)[1]
all.sim.vals$simPerm <- rep(1:rep.iter, times = 50)