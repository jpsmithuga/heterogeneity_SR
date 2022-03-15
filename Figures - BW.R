### FIGURE 2
ltyy <- c(1:3, 4, 6)
pdf(paste0(generalpath, "R Figures/Figure_surface_ests_dstb_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange_fig), 
               as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R', ylab="Dispersion parameter, k", 
               xlim = c(0.01,1.15), 
               ylim = log10(c(0.02,max(krange_fig))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                 points(allests[dstb,"R_hat"], log10(allests[dstb,"k_hat"]), pch = 19, cex = 1, 
                        col = "black")
                 
                 for (i in 1:length(ctlns95[dstb])){
                   lines(ctlns95[dstb][[i]][[2]], ctlns95[dstb][[i]][[3]], lwd = 1.5, lty = ltyy[i], col = "black")
                 }
               })
mtext(expression(bold("A")), adj = 0, cex = 1.5)
legend(0.35-adjx, 0.4+adjy, ncol = 3, legend = leg_ds, bty = "n", cex = 0.7)

 legend(0.605-adjx, 0.34+adjy, rep("", times = 5), 
        lty = ltyy, 
        pch = NA,
        lwd = 1.5, col = "black", bty="n")
 legend(0.80, -0.1+adjy, c("MLE","95% CR"),
        pch = c(19,NA),
        lty = c(NA, 1), lwd = 2, col = "black", bty = "n", cex = 0.85)
dev.off()



pdf(paste0(generalpath, "R Figures/Figure_surface_ests_drtb_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange_fig), 
               as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R', ylab="Dispersion parameter, k", 
               xlim = c(0.01,1.15), 
               ylim = log10(c(0.02,max(krange_fig))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                 points(allests[drtb,"R_hat"], log10(allests[drtb,"k_hat"]), pch = 19, cex = 1, 
                        col = "black")
                 
                 for (i in 1:length(ctlns95[drtb])){
                   lines(ctlns95[drtb][[i]][[2]], ctlns95[drtb][[i]][[3]], lwd = 1.5, lty = ltyy[i], col = "black")
                 }
               })
mtext(expression(bold("A")), adj = 0, cex = 1.5)
legend(0.35-adjx, 0.4+adjy, ncol = 3, legend = leg_dr, bty = "n", cex = 0.7)

legend(0.605-adjx, 0.34+adjy, rep("", times = 4), 
       lty = ltyy, 
       pch = NA,
       lwd = 1.5, col = "black", bty="n")
legend(0.80, -0.1+adjy, c("MLE","95% CR"),
       pch = c(19,NA),
       lty = c(NA, 1), lwd = 2, col = "black", bty = "n", cex = 0.85)
dev.off()



#pdf(paste0(generalpath, "R Figures/Figure_surface_ests_drtb_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange_fig), as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R',ylab="Dispersion parameter, k", 
               xlim=c(0.01, 1.15),
               ylim=log10(c(0.02,max(krange_fig))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                 points(allests[drtb,"R_hat"], log10(allests[drtb,"k_hat"]), pch = 19, cex = 1, col = colors[drtb])
                 
                 for (i in 1:length(ctlns95[drtb])){
                   lines(ctlns95[drtb][[i]][[2]], ctlns95[drtb][[i]][[3]], lwd = 2.5, lty = 1, col=colors[drtb][i])
                   lines(ctlns90[drtb][[i]][[2]], ctlns90[drtb][[i]][[3]], lwd = 2.5, lty = 3, col=colors[drtb][i])
                 }
               })
mtext(expression(bold("B")), adj = 0, cex = 1.5)
legend(0.45-adjx, 0.4+adjy, ncol = 3, legend = leg_dr, bty = "n", cex = 0.8)
legend(0.605-adjx, 0.34+adjy, rep("", times = 4), 
       lty = c(rep(1,length(allests[drtb,1]))), 
       pch = c(rep(19,length(allests[drtb,1]))),
       lwd = 1.5, col = c(colors[drtb]), bty="n")
legend(0.80, -0.1+adjy, c("MLE","95% CR", "90% CR"),
       pch = c(19,NA,NA),
       lty = c(NA, 1,3), lwd = 2, col = "black", bty = "n", cex = 0.85)
#dev.off()


###################################
# Figure 3
###################################
propinfection <- function(R, k, prop){
  xp <- qgamma(1 - prop, shape = k, rate = k/R) 
  tp <- 1 - pgamma(xp, shape = k+1, rate = k/R) 
  return(tp) 
}

ests <- data.frame(study = c("Alvarez (2020)","Asare (2020)","Chee (2020)","Guerra-Assunção (2015)","Jiang (2019)","Macedo (2020)",
                             "Verza (2020)","Walker (2014)","Yang (2017)"),
                   R_hat = c(0.73, 0.21, 0.12, 0.48, 0.16, 0.56, 0.29, 0.10, 0.20),
                   k_hat = c(0.07, 0.15, 0.02, 0.24, 0.24, 0.19, 0.48, 0.06, 0.34))


xx <- seq(0,1,0.001)
prop_infect <- list()
# ests[[7]][2,3] <- 10000000 #approximate infinity
for (i in 1:length(ests)){
  prop_infect[[i]] <- cbind(xx, 
                            propinfection(ests[[i]][1,1], ests[[i]][2,1], xx),
                            round(propinfection(ests[[i]][1,1], ests[[i]][2,1], xx),2))
}

see <- prop_infect[[1]]
see[findInterval(0.50, see[,2]),]

perc50 <- cbind(
  unlist(lapply(prop_infect, function(x) round(x[findInterval(0.5, x[,2]),1]*100))),
  unlist(lapply(prop_infect_lo, function(x) round(x[findInterval(0.5, x[,2]),1]*100))),
  unlist(lapply(prop_infect_hi, function(x) round(x[findInterval(0.5, x[,2]),1]*100)))
)
perc50 <- paste0(perc50[,1], "% (",perc50[,2], "-", perc50[,3], ")")

perc80 <- cbind(
  unlist(lapply(prop_infect, function(x) round(x[findInterval(0.8, x[,2]),1]*100))),
  unlist(lapply(prop_infect_lo, function(x) round(x[findInterval(0.8, x[,2]),1]*100))),
  unlist(lapply(prop_infect_hi, function(x) round(x[findInterval(0.8, x[,2]),1]*100)))
)
perc80 <- paste0(perc80[,1], "% (",perc80[,2], "-", perc80[,3], ")")


perc90 <- cbind(
  unlist(lapply(prop_infect, function(x) round(x[findInterval(0.9, x[,2]),1]*100))),
  unlist(lapply(prop_infect_lo, function(x) round(x[findInterval(0.9, x[,2]),1]*100))),
  unlist(lapply(prop_infect_hi, function(x) round(x[findInterval(0.9, x[,2]),1]*100)))
)
perc90 <- paste0(perc90[,1], "% (",perc90[,2], "-", perc90[,3], ")")



propests <- data.frame(study = c("Alvarez (2020)","Asare (2020)","Chee (2020)",
                                 "Guerra (2015)","Jiang (2019)","Macedo (2020)",
                                 "Verza (2020)","Walker (2014)","Yang (2017)"),
                       perc50 = perc50,
                       perc80 = perc80,
                       perc90 = perc90)


ds <- c(1:2, 4, 7:8)
dr <- setdiff(1:9, ds)

legx <- 0.25
linx <- legx - 0.05
legy <- 0.2
liny <- legy
cexx2 <- 0.80

pdf(paste0(generalpath, "R Figures/Figure_pro_infected_both_", Sys.Date(),".pdf"), width = 6*1.3, height = 10*1.3)
par(mfrow=c(2,1))

plot(xx, xx, type = 'n', xlab = '', ylab = '', axes = FALSE)
for (i in 1:length(prop_infect[dstb])){
  lines(prop_infect[dstb][[i]][,1], prop_infect[dstb][[i]][,2], lty=i, lwd=2, col="black")
  #    polygon(c(xx, rev(xx)), c(prop_infect_lo[dstb][[i]][,2], rev(prop_infect_hi[dstb][[i]][,2])), col = colorst[dstb][i], border = NA)
}
abline(h = c(0.50, 0.80, 0.90), lty = 3, col = "grey40")
axis(side = 1, at=seq(0,1,0.1), labels = seq(0, 100, 10))
axis(side = 2, at=seq(0,1,0.1))
text(c("50% of secondary transmission", "80% of secondary transmission", "90% of secondary transmission"),
     x=0.75, y=c(0.50, 0.80, 0.90)+0.02, col = "black", cex = 0.7)
mtext(side=1, 'Percent of Infectious Cases', padj=4)
mtext(side=2, 'Expected Proportion of Secondary Transmission', padj=-4)
mtext(expression(bold("A")), adj = 0.02, padj = -0.5, cex = 1.2)
legend(x = legx, y = legy, ncol = 4, cex = cexx2, bty = "n",
       legend = unlist(propests[ds,]))
legend(x = linx, y = liny, rep(NA, 5), lty = 1:5, bty = "n", cex = cexx2)

text(x = legx + 0.5, y = legy + 0.08, "Percent of infectious cases\nresponsible for secondary transmission, % (UI)", cex = cexx2)
text(x = legx + 0.25, y = legy + 0.01, "50%", cex = cexx2, font = 2)
text(x = legx + 0.44, y = legy + 0.01, "80%", cex = cexx2, font = 2)
text(x = legx + 0.63, y = legy + 0.01, "90%", cex = cexx2, font = 2)

plot(xx, xx, type = 'n', xlab = '', ylab = '', axes = FALSE)
for (i in 1:length(prop_infect[drtb])){
  lines(prop_infect[drtb][[i]][,1], prop_infect[drtb][[i]][,2], lty=i, lwd=2, col="black")
  #    polygon(c(xx, rev(xx)), c(prop_infect_lo[drtb][[i]][,2], rev(prop_infect_hi[drtb][[i]][,2])), col = colorst[drtb][i], border = NA)
}
abline(h = c(0.50, 0.80, 0.90), lty = 3, col = "grey40")
axis(side = 1, at=seq(0,1,0.1), labels = seq(0, 100, 10))
axis(side = 2, at=seq(0,1,0.1))
text(c("50% of secondary transmission", "80% of secondary transmission", "90% of secondary transmission"),
     x=0.75, y=c(0.50, 0.80, 0.90)+0.02, col = "black", cex = 0.7)
mtext(side=1, 'Percent of Infectious Cases', padj=4)
mtext(side=2, 'Expected Proportion of Secondary Transmission', padj=-4)
mtext(expression(bold("A")), adj = 0.02, padj = -0.5, cex = 1.2)
legend(x = legx, y = legy, ncol = 4, cex = cexx2, bty = "n",
       legend = unlist(propests[dr,]))
legend(x = linx, y = liny, rep(NA, 4), lty = 1:4, bty = "n", cex = cexx2)
text(x = legx + 0.5, y = legy + 0.08, "Percent of infectious cases\nresponsible for secondary transmission, % (UI)", cex = cexx2)
text(x = legx + 0.25, y = legy + 0.01, "50%", cex = cexx2, font = 2)
text(x = legx + 0.44, y = legy + 0.01, "80%", cex = cexx2, font = 2)
text(x = legx + 0.63, y = legy + 0.01, "90%", cex = cexx2, font = 2)
dev.off()


###################################
# Figure 5
###################################

datee <- "2021-11-09"

load(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_estimates_true_",datee,".RData"))
load(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data/simulated_estimates_missing_",datee,".RData"))

ests <- data.frame(study = c("Alvarez (2020)","Asare (2020)","Chee (2020)","Guerra-Assunção (2015)","Jiang (2019)","Macedo (2020)",
                             "Verza (2020)","Walker (2014)","Yang (2017)"),
                   R_hat = c(0.73, 0.21, 0.12, 0.48, 0.16, 0.56, 0.29, 0.10, 0.20),
                   k_hat = c(0.07, 0.15, 0.02, 0.24, 0.24, 0.19, 0.48, 0.06, 0.34))

num_sims <- 1000
num_chains <- 2000
p1 <- 0.60
p2 <- 0.20
p_cens <- 0.10
pchh <- rep(5, 9)# c(0, 1, 2, 4, 5, 6, 7, 13, 11)

initials <- c("Al", "As", "Ch", "GA", "Ji", "Ma", "Ve", "Wa", "Ya")

init.cex <- 0.6

#pngscale <- 0.7
pdf(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/R Figures/simulations_p1_",p1,"_p2_",p2,"_perfimperf_bw_", Sys.Date(),".pdf"), width = 8, height = 10)
#png(paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/R Figures/simulations_p1_",p1,"_p2_",p2,"_perfimperf", Sys.Date(),".png"), width = 800*pngscale, height = 1000*pngscale)
par(mfrow = c(2,1))
# perfect surveillance
plot(c(0,1),c(0,1),type = 'n', xlab = '', ylab = '', axes = TRUE, frame = FALSE)
for (i in 1:nrow(ests)){
  segments(
    x0 = quantile(allsimests.miss[[i]][ ,3], probs = 0.25), 
    x1 = quantile(allsimests.miss[[i]][ ,3], probs = 0.75),
    y0 = median(allsimests.miss[[i]][ ,2]), col = "black") # recall [,3] is R and [,2] is k
  segments(
    y0 = quantile(allsimests.miss[[i]][ ,2], probs = 0.25), 
    y1 = quantile(allsimests.miss[[i]][ ,2], probs = 0.75),
    x0 = median(allsimests.miss[[i]][ ,3]), col = "black") # recall [,3] is R and [,2] is k
  points(median(allsimests.miss[[i]][ ,3]), median(allsimests.miss[[i]][ ,2]), pch = pchh[i], col = "black")
}
text(x = ests[ ,2]+0.025, y = ests[ ,3], initials, cex = init.cex)
points(ests[ ,2], ests[ ,3], pch = 19, col = "black", cex = 0.75)
#axis(1, at = seq(0,1,0.1), labels = seq(0,1,0.1))
#axis(2, at = seq(0,1,0.1), labels = seq(0,1,0.1))
mtext("R estimates", 1, padj=4)
mtext("k estimates", 2, padj=-4)
legend("topleft",c("True simulated value", "Median MLE Value", "Interquartile Range"), pch=c(16, 5, 3), bty="n")
legend(x = 0.70, y = 1, ests[ ,1], col = "black", bty="n", cex = 0.75)
legend(x = 0.64, y = 1, initials, col = "black", bty="n", cex = 0.75)

# imperfect surveillance
# perfect surveillance
plot(c(0,1),c(0,1),type = 'n', xlab = '', ylab = '', axes = FALSE)
for (i in 1:nrow(ests)){
  segments(
    x0 = quantile(allsimests.true[[i]][ ,3], probs = 0.25), 
    x1 = quantile(allsimests.true[[i]][ ,3], probs = 0.75),
    y0 = median(allsimests.true[[i]][ ,2]), col = "black") # recall [,3] is R and [,2] is k
  segments(
    y0 = quantile(allsimests.true[[i]][ ,2], probs = 0.25), 
    y1 = quantile(allsimests.true[[i]][ ,2], probs = 0.75),
    x0 = median(allsimests.true[[i]][ ,3]), col = "black") # recall [,3] is R and [,2] is k
  points(median(allsimests.true[[i]][ ,3]), median(allsimests.true[[i]][ ,2]), pch = pchh[i], col = "black")
  segments(x0 = ests[i ,2], y0 = ests[i ,3],
           x1 = median(allsimests.true[[i]][ ,3]), y1 =median(allsimests.true[[i]][ ,2]),
           lty = 3)
  # text(x = median(allsimests.true[[i]][ ,3])+0.025, 
  #      y = median(allsimests.true[[i]][ ,2])+0.01, initials[i], cex = init.cex)
}
text(x = ests[ ,2]+0.02 , y = ests[ ,3], initials, cex = init.cex)
points(ests[ ,2], ests[ ,3], pch = 19, col = "black", cex = 0.75)
axis(1) 
axis(2)
mtext("R estimates", 1, padj=4)
mtext("k estimates", 2, padj=-4)
legend("topleft",c("True simulated value", "Median MLE Value", "Interquartile Range"), pch=c(16, 5, 3), bty="n")
legend(x = 0.70, y = 1, ests[ ,1], col = "black", bty="n", cex = 0.75)
legend(x = 0.64, y = 1, initials, col = "black", bty="n", cex = 0.75)

dev.off()










plot(c(0,1),c(0,1),type = 'n', xlab = '', ylab = '', axes = FALSE)
for (i in 1:nrow(ests)){
  segments(
    x0 = quantile(allsimests.true[[i]][ ,3], probs = 0.25), 
    x1 = quantile(allsimests.true[[i]][ ,3], probs = 0.75),
    y0 = median(allsimests.true[[i]][ ,2]),  col = "black") # recall [,3] is R and [,2] is k
  segments(
    y0 = quantile(allsimests.true[[i]][ ,2], probs = 0.25), 
    y1 = quantile(allsimests.true[[i]][ ,2], probs = 0.75),
    x0 = median(allsimests.true[[i]][ ,3]), col = "black") # recall [,3] is R and [,2] is k
  points(median(allsimests.true[[i]][ ,3]), median(allsimests.true[[i]][ ,2]), pch = pchh[i], col = "black")
}
points(ests[ ,2], ests[ ,3], pch = 16, col = "black")
axis(1) 
axis(2)
mtext("R estimates",1, padj = 4)
mtext("k estimates",2, padj = -4)
legend("topleft",c("True simulated value", "Interquartile Range"), pch=c(16,3), bty="n")
legend("topright", ests[ ,1], pch = pchh, col = "black", bty="n", cex = 0.75)
dev.off()

# 
# 
# for (i in 1:nrow(ests)){
#   points(allsimests.true[[i]][ ,3], allsimests.true[[i]][ ,2], pch = 16, col = tcolors[i]) # recall [,3] is R and [,2] is k
#   points(median(allsimests.true[[i]][ ,3]), median(allsimests.true[[i]][ ,2]), pch = 4, col = colors[i])
# }
# points(ests[ ,2], ests[ ,3], pch = 16, col = colors[1:9])
# axis(1) 
# axis(2)
# mtext("R estimates",1, padj=4)
# mtext("k estimates",2, padj=-4)
# legend("topleft",c("True simulated value", "Median MLE value"), pch=c(16,4), bty="n")
# legend("topright", ests[,1], pch=16, col = colors[1:9], bty="n", cex = 0.75)
# dev.off()