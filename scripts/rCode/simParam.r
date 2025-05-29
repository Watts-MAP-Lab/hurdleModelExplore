## May 29th 2025
sim.param.nitems <- c(5,10)
sim.param.ncates <- c(3,6)
sim.param.discri <- c(.8,1.8)
sim.param.2pl.spread <- c(3)
sim.param.sample <- c(15000)
sim.param.faccor <- c(.2,.8)
sim.param.difgrmF <- c(-2,0.5)
sim.param.difgrmC <- c(3)
sim.param.dif2pl <- c(-3,0.5)
sim.param.discri2 <- c(1.8,3)
sim.iter <- 1:40
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.discri2,sim.param.dif2pl, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCat", "discrim2pl", "diffSpread", "n", "facCor", "difGrmF","difGrmC","grmDiscrim","dif2PL","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)