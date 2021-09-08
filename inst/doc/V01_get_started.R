## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE, fig.width=7, fig.height=4---------------------
library(transfR)
data(Oudon)
obs <- as_transfr(st = Oudon$obs[,,1:3], hl = Oudon$hl[1:3]) #gauged catchments
sim <- as_transfr(st = Oudon$obs[,,4:6], hl = Oudon$hl[4:6]) #catchments considered as ungauged

## ---- echo=TRUE, message=FALSE, fig.width=7, fig.height=4---------------------
plot(x = obs, i = 1, attribute = "Qobs", format = "%b %d")
plot(obs$hl[[1]], axes = T, main = "Hydraulic length of gauged catchment 1 [m]",
     downsample=1,col=hcl.colors(n=20,palette="Blues"))

## ---- echo=TRUE, message=TRUE, fig.width=7, fig.height=4----------------------
obs <- velocity(obs, method = "loire2016")
obs$uc
sim <- velocity(sim, method = "loire2016")
sim$uc

## ---- echo=TRUE, message=FALSE, results='hide', fig.width=7, fig.height=4-----
obs <- uh(obs)
sim <- uh(sim)
plot(obs, i = 1, attribute = "uh")

## ---- echo=TRUE, message=FALSE, results='hide'--------------------------------
obs <- lagtime(obs)
obs <- rapriori(obs)

## ----Inversion, echo=TRUE, message=FALSE, results='hide'----------------------
obs <- inversion(obs, parallel = TRUE, cores=2)

## ---- fig.show='hold', echo=TRUE, message=FALSE, results='hide'---------------
mdist <- hdist(x = obs, y = sim, method = "rghosh", parallel = TRUE, cores=2)
sim <- mixr(obs = obs, sim = sim, mdist = mdist)

## ---- fig.show='hold', echo=TRUE, message=FALSE, results='hide', fig.width=7, fig.height=4----
sim <- convolution(sim)
plot(x = sim, i = 1, attribute = c("Qobs","Qsim"), 
     ylab = expression(paste("Discharge [",m^3/s,"]")),
     col = c("#a6bddb","#045a8d"), format = "%b %d")

