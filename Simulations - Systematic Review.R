
rm(list=ls())
#### Functions for simulation

# Branching process
bp <- function(gens = 100, init.size = 1, offspring, ...){  
  Z <- list()
  Z[[1]] <- init.size
  i <- 1 
  while(sum(Z[[i]]) > 0 && i <= gens) { 
    Z[[i+1]] <- offspring(sum(Z[[i]]), ...) 
    i <- i+1 
  } 
  return(Z)
} 

#### Imperfect observation - returns both imperfect and counterfactual perfect surveillance system
Missing <- function(true_R, true_k, num_chains = 2000, p1, p2, prob_cens){
  z <- replicate(num_chains,bp(offspring = rnbinom, mu = true_R, size = true_k)) 
  z.pass <- z; z.act <- z                               
  
  ## - - - - - - - - - - - - - - - - - - - - -
  ## Passive surveillance
  ## - - - - - - - - - - - - - - - - - - - - -
  for (i in 1:length(z.pass)){ 
    for (j in 1:length(z.pass[[i]])){
      for (k in 1:length(z.pass[[i]][[j]])){
        for (l in 1:length(z.pass[[i]][[j]][[k]])){
          if (runif(1) < (1 - p1)){                     
            z.pass[[i]][[j]][[k]] <- NA 
          }}}}}
  ## - - - - - - - - - - - - - - - - - - - - -
  ## Active case finding (only chains with at least one case observed by passive surveillance)
  ## - - - - - - - - - - - - - - - - - - - - -
  a <- z.pass
  for (i in 1:length(a)){
    a[[i]][[1]] <- NA
  }
  b <- cbind(a, sapply(a, function(x) all(is.na(unlist(x))))) # Determine if at least one case in the chain has been seen
  for (i in 1:length(z.pass)){
    if (b[i,2] == TRUE){  
      z.act[[i]] <- a[[i]] # Skip if cluster no cases are observed
    } else {              
      for (j in 1:length(z.pass[[i]])){
        for (k in 1:length(z.pass[[i]][[j]])){
          for (l in 1:length(z.pass[[i]][[j]][[k]])){
            if (is.na(z.pass[[i]][[j]][[k]])){      
              if(runif(1) <= (p2)){       #Active probability of being seen by case detection         
                z.act[[i]][[j]][[k]] <- z[[i]][[j]][[k]] #Reassign original value
              } else {z.act[[i]][[j]][[k]] <- z.pass[[i]][[j]][[k]]}
            }}}}}}
 
  # "Break" chains based on the position of missing cases in the chain
  l <- z.act  #dummy/temp data to not change z.act
  for (i in 1:length(l)){ 
    l[[i]][[1]] <- NULL #remove first position of the nested list so that it eases summing lengths (can't sum based on integer values in imperfect observations)
  }
  #"Break apart" the chains 
  t1 <- lapply(lapply(seq_along(l), function(nm) {split(l[[nm]], cumsum(sapply(l[[nm]], function(x) all(is.na(x)))))}), function(lstA) lapply(lstA,function(x) Filter(function(y) !all(is.na(y)), x)))
  t2 <- rapply(unlist(t1,recursive=FALSE),function(x) x[!is.na(x)], how="replace") #Remove NA values. 
  z.broken <- Filter(length,t2) #remove all with length 0 (missing/unobserved)
  
  ## - - - - - - - - - - - - - - - - - - - - -
  ## Censoring
  ## - - - - - - - - - - - - - - - - - - - - -
  z.cen <- z.broken                          # Initialize the censored list
  for (i in 1:length(z.broken)){             # Iterate through the list
    if (length(z.broken[[i]]) > 1) {           # List must have at least length of two (cant be censored if the index case isnt seen, then it is unobserved as above)
      if(runif(1) <= prob_cens){               # Stochastic process to determine if the nested list will be censored
        if(length(z.broken[[i]]) == 2){
          n <- 2} else {                     # A little trick to get over the issue of sample(2:2,1) returning values of 1 as well as 2
            n <- sample(2:length(z.broken[[i]]), 1)}   # Randomly determine what list position in the nested list will be the censor threshold
        z.cen[[i]][n:length(z.broken[[i]])] <- NA   # Fill all positions from n to the end of the nested list with NA
      }}}
  out_list <- lapply(z.cen, function(x) {    # Remove all nested list elements that contain NA 
    inds <- sapply(x, function(x) any(is.na(x)))
    if(any(inds)) x[seq_len(which.max(inds) - 1)] else x})
  cens <- numeric(length(out_list))
  true <- numeric(length(out_list))
  for (k in 1:length(out_list)){
    cens[k] <- sum(lengths(out_list[[k]])) # Get cluster size of censored clusters
    true[k] <- sum(lengths(z.broken[[k]])) # Get cluster size of uncensored (but imperfect obs) clusters
  }
  Y_cens <- data.frame(y.cens=cens, censor=ifelse(cens!=true,1,0)) #Create a censoring index (1=censored, 0=uncensored)

  
  #  z.fin <- rapply(z.act,function(x) x[!is.na(x)], how = "replace") 
  # j1 <- z.fin 
  # for (i in 1:length(j1)){
  #   j1[[i]][[1]] <- NULL
  # }
  # j2 <- integer()
  # for (k in 1:length(j1)){
  #   j2[[k]] <- sum(lengths(j1[[k]]))
  # }
  # y.imp <- j2[j2 > 0] 
   y.true <- unlist(lapply(z, function(x) sum(unlist(x)))) 
  return(list(y.true, Y_cens[,1]))
  }


##### Simulate Data
# manually encode parameter estimates for convenience
ests <- data.frame(study = c("Alvarez (2020)","Asare (2020)","Chee (2020)","Guerra-Assunção (2015)","Jiang (2019)","Macedo (2020)",
                             "Verza (2020)","Walker (2014)","Yang (2017)"),
                   R_hat = c(0.73, 0.21, 0.12, 0.48, 0.16, 0.56, 0.29, 0.10, 0.20),
                   k_hat = c(0.07, 0.15, 0.02, 0.24, 0.24, 0.19, 0.48, 0.06, 0.34))

num_sims <- 1000
num_chains <- 2000
p1 <- 0.60
p2 <- 0.20
p_cens <- 0.10

allsims.missing <- vector("list", length = nrow(ests))
allsims.true <- vector("list", length = nrow(ests))
for (i in 1:nrow(ests)){
  tmp <- replicate(num_sims, Missing(ests[i,"R_hat"], ests[i,"k_hat"], num_chains, p1, p2, p_cens))
  imperfect.indx <- seq(1 ,length(tmp), 2)
  perfect.indx <- seq(2 ,length(tmp), 2)
  
  allsims.missing[[i]] <- tmp[imperfect.indx]
  allsims.true[[i]] <- tmp[perfect.indx]
}

names(allsims.missing) <- names(allsims.true) <- ests[,1]

save(allsims.missing, file = paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_data_missing_",Sys.Date(),".RData"))
save(allsims.true, file = paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_data_true_",Sys.Date(),".RData"))

#### Infer Parameters

clust.likelihood <- function(Y,k,r0){
  nb.likelihood <- function(x){
    lgamma(k*x+(x-1))-(lgamma(k*x)+lgamma(x+1))+(x-1)*log(r0/k)-(k*x+(x-1))*log(1+r0/k)
  }
  calc_sizes <- unique(c(1,Y))
  likelihoods <- c()
  likelihoods[calc_sizes] <- sapply(calc_sizes, nb.likelihood)
  sexpl <- sum(exp(likelihoods),na.rm=TRUE)
  if (sexpl < 1) {
    maxl<-log(1 - sum(exp(likelihoods), na.rm=TRUE))
  } else {maxl <- -Inf}
  likelihoods <- c(likelihoods, maxl)
  cluster_likelihoods <- likelihoods[Y]
  result <- sum(cluster_likelihoods)
  return(result)
}

param.est <- function(y){
  ret <- c()
  if(is.null(ncol(y)) == FALSE){
    y <- y[ ,1]}
  r0_hat <- 1 - (1 / mean(y)) 
  if (r0_hat < 1e-9) {r0_hat <- 1e-9}
  t <- optim(par = c(0.01,0.01), lower = c(1e-9,1e-9), fn = function(x){-clust.likelihood(y, x[1], r0_hat)}, method = "L-BFGS-B") #Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows constraints; each variable can be given a lower and/or upper bound
  k_hat<-t$par[[1]]
  ret$k_hat <- k_hat; ret$r0_hat <- r0_hat
  ret <- unlist(ret)
  return(ret)
}

allsimests.miss <- allsimests.true <- vector("list", length = nrow(ests))

names(allsimests.miss) <- names(allsimests.true) <- ests[,1]

for (i in 1:nrow(ests)){
  allsimests.miss[[i]] <- data.frame(matrix(NA, ncol = 3, nrow = num_sims))
  allsimests.true[[i]] <- data.frame(matrix(NA, ncol = 3, nrow = num_sims))
}

for (i in 1:nrow(ests)){
  for (j in 1:num_sims){
    allsimests.miss[[i]][j, 1] <- ests[i,1]
    allsimests.miss[[i]][j, 2] <- param.est(allsims.missing[[i]][[j]])[[1]]
    allsimests.miss[[i]][j, 3] <- param.est(allsims.missing[[i]][[j]])[[2]]
    
    allsimests.true[[i]][j, 1] <- ests[i,1]
    allsimests.true[[i]][j, 2] <- param.est(allsims.true[[i]][[j]])[[1]]
    allsimests.true[[i]][j, 3] <- param.est(allsims.true[[i]][[j]])[[2]]
  }
}


save(allsimests.miss, file = paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_estimates_missing_",Sys.Date(),".RData"))
save(allsimests.true, file = paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_estimates_true_",Sys.Date(),".RData"))




# pch_vals <- 1:length(allsimests)

makeTransparent <- function(someColor, alpha = 100){
  if(alpha < 1){alpha <- alpha*100} else{alpha <- alpha}
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red = curcoldata[1], green = curcoldata[2], blue = curcoldata[3], alpha = alpha, maxColorValue = 255)})
}

colors <- ggsci::pal_npg(palette = c("nrc"), alpha = 1)(9)
tcolors <- makeTransparent(colors, alpha = 2)

pngscale <- 0.7
#pdf(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/R Figures/simulations_p1_",p1,"_p2_",p2,"_perfimperf", Sys.Date(),".pdf"), width = 8, height = 10)
png(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/R Figures/simulations_p1_",p1,"_p2_",p2,"_perfimperf", Sys.Date(),".png"), width = 800*pngscale, height = 1000*pngscale)
par(mfrow=c(2,1))
 # perfect surveillance
plot(c(0,1),c(0,1),type = 'n', xlab = '', ylab = '', axes = FALSE)
  for (i in 1:nrow(ests)){
    points(allsimests.miss[[i]][ ,3], allsimests.miss[[i]][ ,2], pch = 16, col = tcolors[i]) # recall [,3] is R and [,2] is k
    points(median(allsimests.miss[[i]][ ,3]), median(allsimests.miss[[i]][ ,2]), pch = 4, col = colors[i])
  }
  points(ests[ ,2], ests[ ,3], pch = 16, col = colors[1:9])
  axis(1) 
  axis(2)
  mtext("R estimates", 1, padj=4)
  mtext("k estimates", 2, padj=-4)
  legend("topleft",c("True simulated value", "Median MLE value"), pch=c(16,4), bty="n")
  legend("topright", ests[ ,1], pch = 16, col = colors[1:9], bty="n", cex = 0.75)
 # imperfect surveillance
plot(c(0,1),c(0,1),type = 'n', xlab = '', ylab = '', axes = FALSE)
  for (i in 1:nrow(ests)){
    points(allsimests.true[[i]][ ,3], allsimests.true[[i]][ ,2], pch = 16, col = tcolors[i]) # recall [,3] is R and [,2] is k
    points(median(allsimests.true[[i]][ ,3]), median(allsimests.true[[i]][ ,2]), pch = 4, col = colors[i])
  }
  points(ests[ ,2], ests[ ,3], pch = 16, col = colors[1:9])
  axis(1) 
  axis(2)
  mtext("R estimates",1, padj=4)
  mtext("k estimates",2, padj=-4)
  legend("topleft",c("True simulated value", "Median MLE value"), pch=c(16,4), bty="n")
  legend("topright", ests[,1], pch=16, col = colors[1:9], bty="n", cex = 0.75)
dev.off()



#### Sample Size figures

num_sims <- 1000
num_chains <- c(2000, 1000, 500, 100)
p1 <- 1
p2 <- 1
p_cens <- 0


R <- 0.50
k <- 0.15


nsims <- vector("list", length = length(num_chains))

for (i in 1:length(num_chains)){
  tmp <- replicate(num_sims, Missing(R, k, num_chains[i], 1, 1, 0))
  perfect.indx <- seq(2 ,length(tmp), 2)
  nsims[[i]] <- tmp[perfect.indx]
}

nsims.ests <-  vector("list", length = length(num_chains))

for (i in 1:length(num_chains)){
  nsims.ests[[i]] <- data.frame(matrix(NA, ncol = 3, nrow = num_sims))
}

for (i in 1:length(num_chains)){
  for (j in 1:num_sims){
    nsims.ests[[i]][j, 1] <- num_chains[i]
    nsims.ests[[i]][j, 2] <- param.est(nsims[[i]][[j]])[[1]]
    nsims.ests[[i]][j, 3] <- param.est(nsims[[i]][[j]])[[2]]
  }
}

axixx <- seq(0,1,0.1)
pchh <- 1:4

pdf(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/R Figures/sim_samplesize_", Sys.Date(),".pdf"), width = 6, height = 6)
plot(c(0.25,0.75), c(0,0.5),type = 'n', xlab = '', ylab = '', axes = FALSE)
for (i in 1:length(num_chains)){
  points(nsims.ests[[i]][ ,3], nsims.ests[[i]][ ,2], pch = 16, cex = 0.75, pch = pchh[i], col = tcolors[i]) # recall [,3] is R and [,2] is k
  points(median(nsims.ests[[i]][ ,3]), median(nsims.ests[[i]][ ,2]), pch = pchh[i], col = colors[i])
  segments(x0 = quantile(nsims.ests[[i]][ ,3], probs=0.25), x1 = quantile(nsims.ests[[i]][ ,3], probs = 0.75), y0 = median(nsims.ests[[i]][ ,2]), col = colors[i])
  segments(y0 = quantile(nsims.ests[[i]][ ,2], probs=0.25), y1 = quantile(nsims.ests[[i]][ ,2], probs=0.75), x0 = median(nsims.ests[[i]][ ,3]), col = colors[i])
}
points(R, k, pch = 16, col = "black", cex = 0.75)
axis(1, at = axixx, label = axixx) 
axis(2, at = axixx, label = axixx)
mtext("R estimates",1, padj=4)
mtext("k estimates",2, padj=-4)
legend("topleft",c("True simulated value", "N = 2000","N = 1000","N = 500","N = 100"), pch=c(16,pchh), bty="n", col = c("black", colors))

dev.off()



