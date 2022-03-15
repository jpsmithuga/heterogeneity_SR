# rm(list=ls())

# set libraries
 library(ggplot2); library(cowplot)
# define paths
functionpath <- "~jonathansmith/Dropbox/CDC/KOPANYO/Heterogeneity (Kopanyo)/"
datpath <- "~jonathansmith/Dropbox/Emory/Systematic Review - k/Abstracted Cluster Data"
generalpath <- "~jonathansmith/Dropbox/Emory/Systematic Review - k/"

# Pull in likelihood functions
source(paste0(functionpath,"R Code/Kopanyo Heterogeneity - Likelihood Functions.R"))

#Pull in Data
allfiles <- list.files(path = datpath, pattern = "*.csv", full.names = TRUE)

# Calculations
  # Set R and k range
R1 <- 0.01; R2 <- 1.55; Rrange <- seq(R1, R2, by = 0.01)
k1 <- 0.01; k2 <- 3.05; krange <- c(seq(k1, k2, by = 0.01),100)

surflikes <- maxsurfs <- ests <- list()
studyname <- c()
for (i in 1:length(allfiles)){
   surflikes[[i]] <- surflike(read.csv(allfiles[i]), Rrange = Rrange, krange = krange)
   maxsurfs[[i]] <- surflikes[[i]] == max(surflikes[[i]])
   ests[[i]] <- calc_profile(surflikes[[i]], maxsurfs[[i]], Rrange = Rrange, krange = krange)
}

k_hat <- k_lo <- k_hi <- R_hat <- R_lo <- R_hi <- c()
for (i in 1:length(ests)){
  k_hat[i] <- ests[[i]][2,1]
  k_lo[i] <-  ests[[i]][2,2]
  k_hi[i] <-  ests[[i]][2,3]
  R_hat[i] <- ests[[i]][1,1]
  R_lo[i] <-  ests[[i]][1,2]
  R_hi[i] <-  ests[[i]][1,3]
  studyname[i] <- gsub(".*/(.*)\\..*", "\\1", allfiles[i], perl = TRUE)
}

tbtype <- ifelse(studyname %in% c("Chee_2021","Jiang_2019","Macedo_2019","Yang_2017"),"DR","ALL")

allests <- data.frame(studyname, tbtype, k_hat, k_lo, k_hi, R_hat, R_lo, R_hi)
allests[allests == 100] <- Inf

allests.forexport <- data.frame(studyname, tbtype,
                                k_val = paste0(k_hat," (", k_lo,"-",k_hi,")"),
                                R_val = paste0(R_hat," (", R_lo,"-",R_hi,")"))
write.csv(allests.forexport, paste0("~jonathansmith/Dropbox/Emory/Systematic Review - k/prelim_ests_",Sys.Date(),".csv"), row.names = FALSE)


# Assign colors
makeTransparent <- function(someColor, alpha = 100){
  if(alpha<1){alpha <- alpha*100} else{alpha <- alpha}
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

colors <- ggsci::pal_npg(palette = c("nrc"), alpha = 1)(9)

col_alvarez <- colors[1]
col_asare <- colors[2]
col_chee <- colors[3]
col_ga <- colors[4]
col_jiang <- colors[5]
col_macedo <- colors[6]
col_verza <- colors[7]
col_walker <- colors[8]
col_yang <- colors[9]

alphaForTransparentColor <- 60
colorst <- makeTransparent(colors, alpha = alphaForTransparentColor)
colt_alvarez <- makeTransparent(col_alvarez, alpha = alphaForTransparentColor)
colt_asare <- makeTransparent(col_asare, alpha = alphaForTransparentColor)
colt_chee <- makeTransparent(col_chee, alpha = alphaForTransparentColor)
colt_ga <- makeTransparent(col_ga, alpha = alphaForTransparentColor)
colt_jiang <- makeTransparent(col_jiang, alpha = alphaForTransparentColor)
colt_macedo <- makeTransparent(col_macedo, alpha = alphaForTransparentColor)
colt_verza <- makeTransparent(col_verza, alpha = alphaForTransparentColor)
colt_walker <- makeTransparent(col_walker, alpha = alphaForTransparentColor)
colt_yang <- makeTransparent(col_yang, alpha = alphaForTransparentColor)


########################### Figure

chi95 <- qchisq(0.95, df = 1)/2 
ctlns95 <- list()
for (i in 1:length(allfiles)){
  ctlns95[[i]] <- contourLines(x = Rrange, y = log10(krange), as.matrix(surflikes[[i]] - max(surflikes[[i]])),
                               levels = c(-chi95, chi95+0.001))
}
ctlns95 <- unlist(ctlns95, recursive=FALSE)

chi90 <- qchisq(0.90, df = 1)/2 
ctlns90 <- list()
for (i in 1:length(allfiles)){
  ctlns90[[i]] <- contourLines(x = Rrange, y = log10(krange), as.matrix(surflikes[[i]] - max(surflikes[[i]])),
                              levels = c(-chi90, chi90+0.001))
}
ctlns90 <- unlist(ctlns90, recursive=FALSE)

xtick <- seq(0, 1.4, 0.1)
#ytick <- c(0.01,0.05,0.2,0.1,0.5,1,2)
ytick <- c(0.01,0.05,0.1,0.5,1,5,10,50,100)

dstb <- allests$tbtype == "ALL"
drtb <- allests$tbtype == "DR"

all_names <- c("Alvarez (2020)","Asare (2020)","Chee (2020)","Guerra-Assunção (2015)","Jiang (2019)","Macedo (2020)",
               "Verza (2020)","Walker (2014)","Yang (2017)")
allests$labels <- all_names
dstb_names <- all_names[dstb]
drtb_names <- all_names[drtb]
str(allests)
allests

allests$income <- c("High","Lower-Middle","High","Low","Upper-Middle","High","Upper-Middle","High","Upper-Middle")
allests$incidence <- c(205, 165, 40, round(mean(c(87,124))), 65, 23,  23.7, 8.4,  65)

cor.test(allests$incidence[dstb], allests$k_hat[dstb])
cor.test(allests$incidence[drtb], allests$k_hat[drtb])
cor.test(allests$gdp[dstb], allests$k_hat[dstb])
cor.test(allests$gdp[drtb], allests$k_hat[drtb])

allests$gdp[dstb] <- c(43241.61783,
                       2328.534642,
                       625.2941292,
                       6796.844542,
                       40284.63846)/1000
allests$gdp[drtb] <- c(59797.75218,
                       10500.39562,
                       22439.87687,
                       10500.39562)/1000

allests[7,5] <- 1.55



padjj <- -3
cexx <- 0.75
pchh <- c(19, 2)

pdf(paste0(generalpath, "R Figures/Figure_inc_gdp_vs_k_combined_", Sys.Date(),".pdf"), width = 10, height = 10)

par(mfrow = c(2,2))

plot(allests$incidence[dstb], allests$k_hat[dstb], pch = pchh[1], xlim = c(0,300), ylim = c(0,1), 
       xlab = "Incidence (per 100,000)", ylab = NA, col = colors[dstb])
  segments(x0 = allests$incidence[dstb], y0 = allests$k_lo[dstb], y1 = allests$k_hi[dstb], col = colors[dstb])
  points(allests$incidence[drtb], allests$k_hat[drtb], pch = pchh[2], col = colors[drtb])
  segments(x0 = allests$incidence[drtb], y0 = allests$k_lo[drtb], y1 = allests$k_hi[drtb], col = colors[drtb])
  mtext("Dispersion parameter, k", side = 2, padj = padjj)
  mtext(expression(bold("A")), side = 3, adj = 0.02, padj = -0.5, cex = 1.2)
  legend("topright", c("All TB", dstb_names, NA, "DRTB Only", drtb_names), 
         bty="n", pch = c(NA,rep(pchh[1],5),NA,NA,rep(pchh[2],4)), col = c(NA,colors[dstb],NA,NA,colors[drtb]),
         text.font = c(2,rep(1,5),NA,2,1,1,1,1), cex = cexx)
  
plot(allests$incidence[dstb], allests$R_hat[dstb], pch = pchh[1], xlim = c(0,300), ylim = c(0,1), 
       xlab = "Incidence (per 100,000)", ylab = NA, col = colors[dstb])
  segments(x0 = allests$incidence[dstb], y0 = allests$R_lo[dstb], y1 = allests$R_hi[dstb], col = colors[dstb])
  points(allests$incidence[drtb], allests$R_hat[drtb], pch = pchh[2], col = colors[drtb])
  segments(x0 = allests$incidence[drtb], y0 = allests$R_lo[drtb], y1 = allests$R_hi[drtb], col = colors[drtb])
  mtext("Reproductive number, R", side = 2, padj = padjj)
  mtext(expression(bold("B")), side = 3, adj = 0.02, padj = -0.5, cex = 1.2)
  legend("topright", c("All TB", dstb_names, NA, "DRTB Only", drtb_names), 
         bty="n", pch = c(NA,rep(pchh[1],5),NA,NA,rep(pchh[2],4)), col = c(NA,colors[dstb],NA,NA,colors[drtb]),
         text.font = c(2,rep(1,5),NA,2,1,1,1,1), cex = cexx-0.1)


#gdp
plot(allests$gdp[dstb], allests$k_hat[dstb], pch = pchh[1], xlim = c(0,300), ylim = c(0,1), 
       xlab = "GDP (per 100,000 population)", ylab = NA, col = colors[dstb])
  segments(x0 = allests$gdp[dstb], y0 = allests$k_lo[dstb], y1 = allests$k_hi[dstb], col = colors[dstb])
  points(allests$gdp[drtb], allests$k_hat[drtb], pch = pchh[2], col = colors[drtb])
  segments(x0 = allests$gdp[drtb], y0 = allests$k_lo[drtb], y1 = allests$k_hi[drtb], col = colors[drtb])
  mtext("Dispersion parameter, k", side = 2, padj = padjj)
  mtext(expression(bold("C")), side = 3, adj = 0.02, padj = -0.5, cex = 1.2)
  legend("topright", c("All TB", dstb_names, NA, "DRTB Only", drtb_names), 
         bty="n", pch = c(NA,rep(pchh[1],5),NA,NA,rep(pchh[2],4)), col = c(NA,colors[dstb],NA,NA,colors[drtb]),
         text.font = c(2,rep(1,5),NA,2,1,1,1,1), cex = cexx)

  
plot(allests$gdp[dstb], allests$R_hat[dstb], pch = pchh[1], xlim = c(0,300), ylim = c(0,1), 
       xlab = "GDP (per 100,000 population)", ylab = NA, col = colors[dstb])
  segments(x0 = allests$gdp[dstb], y0 = allests$R_lo[dstb], y1 = allests$R_hi[dstb], col = colors[dstb])
  points(allests$gdp[drtb], allests$R_hat[drtb], pch = pchh[2], col = colors[drtb])
  segments(x0 = allests$gdp[drtb], y0 = allests$R_lo[drtb], y1 = allests$R_hi[drtb], col = colors[drtb])
  mtext("Reproductive number, R", side = 2, padj = padjj)
  mtext(expression(bold("D")), side = 3, adj = 0.02, padj = -0.5, cex = 1.2)
  legend("topright", c("All TB", dstb_names, NA, "DRTB Only", drtb_names), 
         bty="n", pch = c(NA,rep(pchh[1],5),NA,NA,rep(pchh[2],4)), col = c(NA,colors[dstb],NA,NA,colors[drtb]),
         text.font = c(2,rep(1,5),NA,2,1,1,1,1), cex = cexx)

dev.off()
  
  
## Original

pdf(paste0(generalpath, "R Figures/Figure_inc_gdp_vs_k_", Sys.Date(),".pdf"), width = 10, height = 10)
par(mfrow=c(2,2))
  
plot(allests$incidence[drtb], allests$k_hat[drtb], pch = 19, xlim = c(0,300), ylim = c(0,1), 
     xlab = "Incidence (per 100,000)", ylab = NA, col = colors[drtb])
  segments(x0 = allests$incidence[drtb], y0 = allests$k_lo[drtb], y1 = allests$k_hi[drtb], col = colors[drtb], lty = 1)
  segments(x0 = allests$incidence[drtb][2], y0 = allests$k_lo[drtb][2], y1 = allests$k_hi[drtb][2], col = colors[drtb][2])
  mtext(expression(hat(k)), side = 2, padj = -2.5)
  legend("topright", c(drtb_names, NA, "r = 0.57","p = 0.42"), bty="n", pch = c(rep(19, 4), NA, NA, NA), col = c(colors[drtb], NA, NA, NA))

plot(allests$gdp[dstb], allests$k_hat[dstb], pch = 19, xlim = c(0,100), ylim = c(0,1), 
       xlab = "Per Capita GDP (USD, in thousands)", ylab = NA, col = colors[dstb])
  segments(x0=allests$gdp[dstb], y0 = allests$k_lo[dstb], y1=allests$k_hi[dstb], col = colors[dstb])
  mtext(expression(hat(k)), side = 2, padj = -2.5)
  legend("topright", c(dstb_names, NA, "r = -0.65","p = 0.24"), bty="n", pch = c(rep(19, 5), NA, NA, NA), col = c(colors[dstb], NA, NA, NA))

plot(allests$gdp[drtb], allests$k_hat[drtb], pch = 19, xlim = c(0,100), ylim = c(0,1), 
       xlab =  "Per Capita GDP (USD, in thousands)", ylab = NA, col = colors[drtb])
  segments(x0 = allests$gdp[drtb], y0 = allests$k_lo[drtb], y1 = allests$k_hi[drtb], col = colors[drtb], lty = 1)
  segments(x0 = allests$gdp[drtb][2], y0 = allests$k_lo[drtb][2], y1 = allests$k_hi[drtb][2], col = colors[drtb][2])
  mtext(expression(hat(k)), side = 2, padj = -2.5)
  legend("topright", c(drtb_names, NA, "r = -0.94","p = 0.06"), bty="n", pch = c(rep(19, 4), NA, NA, NA), col = c(colors[drtb], NA, NA, NA))
  
    
  
  # plot(allests$gdp[dstb], allests$k_hat[dstb], pch = 19, xlim = c(0,100), ylim = c(0,1), xlab = "Per Capita GDP (USD, in thousands)", ylab = expression(hat(k)))
  # text(8, 0.9, "r= -0.65\np=0.24" )
  # plot(allests$gdp[drtb], allests$k_hat[drtb], pch = 19, xlim = c(0,100), ylim = c(0,1), xlab = "Per Capita GDP (USD, in thousands)", ylab = expression(hat(k)))
  # text(8, 0.9, "r= -0.94\np=0.06" )
  # 
dev.off()


summary(allests$R_hat)
summary(allests$k_hat)


allests[drtb,1]

all_k <- ggplot(data = allests, aes(y = reorder(labels, k_hat), x = k_hat, xmin = k_lo, xmax = k_hi, color = tbtype)) + #, linetype = tbtype)) +
            geom_point() + xlim(c(0,1.55)) +
            geom_errorbarh(height=0.05) +
            scale_color_manual(values = c("black","darkred")) +
            labs(title='', x = "Dispersion parameter, k", y = "") +
            geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
            theme_classic() + theme(legend.position = "none")

all_R <- ggplot(data = allests, aes(y = reorder(labels, R_hat), x = R_hat, xmin = R_lo, xmax = R_hi, color = tbtype)) + #, linetype = tbtype)) +
            geom_point() + xlim(c(0,1.55)) +
            geom_errorbarh(height=0.05) +
            scale_color_manual(values = c("black","darkred")) +
            labs(title='', x = "Reproductive number, R", y = "") +
            geom_vline(xintercept = 1, color='black', linetype='dashed', alpha=.5) +
            theme_classic() + theme(legend.position = "none")

Rk_grid <- plot_grid(all_k, all_R, labels = "AUTO", nrow = 2)

save_plot(paste0(generalpath,"R Figures/forest_all_",Sys.Date(),".png"), Rk_grid, ncol = 1, nrow = 2)
save_plot(paste0(generalpath,"R Figures/forest_all_",Sys.Date(),".pdf"), Rk_grid, ncol = 1, nrow = 2)


ds_k <- ggplot(data = allests[dstb,], aes(y=labels, x=k_hat, xmin=k_lo, xmax=k_hi)) +
          geom_point() + xlim(c(0,1.55)) +
          geom_errorbarh(height=0.05) +
          labs(title='', x=expression(hat(k)), y = '') +
          geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
          theme_classic()
ds_R <- ggplot(data = allests[dstb,], aes(y=labels, x=R_hat, xmin=R_lo, xmax=R_hi)) +
          geom_point() + xlim(c(0,1.55)) +
          geom_errorbarh(height=0.05) +
          labs(title='', x=expression(hat(R)), y = '') +
          geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
          theme_classic()

dr_k <- ggplot(data = allests[drtb,], aes(y=labels, x=k_hat, xmin=k_lo, xmax=k_hi)) +
          geom_point() + xlim(c(0,1.55)) +
          geom_errorbarh(height=0.05) +
          #scale_y_continuous(breaks=1:nrow(allests[drtb,]), labels=allests[drtb,]$labels) +
          labs(title='', x=expression(hat(k)), y = '') +
          geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
          theme_classic()
dr_R <- ggplot(data = allests[drtb,], aes(y=labels, x=R_hat, xmin=R_lo, xmax=R_hi)) +
          geom_point() + xlim(c(0,1.55)) +
          geom_errorbarh(height=0.05) +
          labs(title='', x=expression(hat(R)), y = '') +
          geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
          theme_classic()


plot_grid(ds_k, ds_R, dr_k, dr_R, labels = c("A","B","C","D"))

krange_fig <- krange
krange_fig[krange_fig>4] <- 3.06



see <- as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]]))

legcol_dstb <- legend(0.35, 0.4, rep("", times = 5), 
                           lty = c(rep(1,length(allests[dstb,1])), NA, 1,3), 
                           pch = c(rep(19,length(allests[dstb,1])), NA, NA, NA),
                           lwd = 1.5, col = c(colors[dstb], NA, "black", "black"), bty="n")
dstb_names[3] <- "Guerra (2015)"

# allests[7,7] <- "0.20"
# allests[8,6] <- "0.10"
# allests[8,8] <- "0.20"
# allests[9,6] <- "0.20"
leg_ds <- c(NA,paste0(dstb_names, "         "), 
             expression(bold(hat(R) ~ "(95% CI)")),
             paste0(allests[dstb,6], " (",allests[dstb, 7],"-",allests[dstb, 8],")"),
             expression(bold(hat(k) ~ "(95% CI)")),
             paste0(allests[dstb,3], " (",allests[dstb, 4],"-",allests[dstb, 5],")"))

leg_dr <- c(NA, paste0(drtb_names, "         "), 
            expression(bold(hat(R) ~ "(95% CI)")),
            paste0(allests[drtb,6], " (",allests[drtb, 7],"-",allests[drtb, 8],")"),
            expression(bold(hat(k) ~ "(95% CI)")),
            paste0(allests[drtb,3], " (",allests[drtb, 4],"-",allests[drtb, 5],")"))


adjx <- 0.05
adjy <- 0.08
pdf(paste0(generalpath, "R Figures/Figure_surface_ests_dstb_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange_fig), 
               as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R', ylab="Dispersion parameter, k", 
               xlim = c(0.01,1.15), 
               ylim = log10(c(0.02,max(krange_fig))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                 points(allests[dstb,"R_hat"], log10(allests[dstb,"k_hat"]), pch = 19, cex = 1, col = colors[dstb])
                 
                 for (i in 1:length(ctlns95[dstb])){
                   lines(ctlns95[dstb][[i]][[2]], ctlns95[dstb][[i]][[3]], lwd = 2.5, lty = 1, col=colors[dstb][i])
                   lines(ctlns90[dstb][[i]][[2]], ctlns90[dstb][[i]][[3]], lwd = 2.5, lty = 3, col=colors[dstb][i])
                 }
               })
mtext(expression(bold("A")), adj = 0, cex = 1.5)
legend(0.45-adjx, 0.4+adjy, ncol = 3, legend = leg_ds, bty = "n", cex = 0.8)
legend(0.605-adjx, 0.34+adjy, rep("", times = 5), 
       lty = c(rep(1,length(allests[dstb,1]))), 
       pch = c(rep(19,length(allests[dstb,1]))),
       lwd = 1.5, col = c(colors[dstb]), bty="n")
legend(0.80, -0.1+adjy, c("MLE","95% CR", "90% CR"),
       pch = c(19,NA,NA),
       lty = c(NA, 1,3), lwd = 2, col = "black", bty = "n", cex = 0.85)
dev.off()



pdf(paste0(generalpath, "R Figures/Figure_surface_ests_drtb_", Sys.Date(),".pdf"), width = 10, height = 8)
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
dev.off()



##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##' proportion responsible for ongoing transmission
##' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#' Function allows to set any proportion of ongoing transmission, when "prop = 1.0", estimates the proportion 
#' responsible for all (100%) of secondary infections
propinfection <- function(R, k, prop){
  xp <- qgamma(1 - prop, shape = k, rate = k/R) 
  tp <- 1 - pgamma(xp, shape = k+1, rate = k/R) 
  return(tp) 
}

xx <- seq(0,1,0.001)
prop_infect <- prop_infect_lo <- prop_infect_hi <- list()
plot50x <- plot50x.hi <- plot50x.lo <- c()
plot80x <- plot80x.hi <- plot80x.lo <- c()
plot90x <- plot90x.hi <- plot90x.lo <- c()

 # ests[[7]][2,3] <- 10000000 #approximate infinity
for (i in 1:length(ests)){
  prop_infect[[i]] <- cbind(xx,propinfection(ests[[i]][1,1], ests[[i]][2,1], xx))
  prop_infect_lo[[i]] <- cbind(xx,propinfection(ests[[i]][1,2], ests[[i]][2,2], xx))
  prop_infect_hi[[i]] <- cbind(xx,propinfection(ests[[i]][1,3], ests[[i]][2,3], xx))
  
  plot50x[i] <- prop_infect[[i]][min(which(round(prop_infect[[i]][,2],2)==0.50)),1]
  plot80x[i] <- prop_infect[[i]][min(which(round(prop_infect[[i]][,2],2)==0.80)),1]
  plot90x[i] <- prop_infect[[i]][min(which(round(prop_infect[[i]][,2],2)==0.90)),1]
  
  plot50x.lo[i] <- prop_infect_lo[[i]][min(which(round(prop_infect_lo[[i]][,2],2)==0.50)),1]
  plot50x.hi[i] <- prop_infect_hi[[i]][min(which(round(prop_infect_hi[[i]][,2],2)==0.50)),1]
  plot80x.lo[i] <- prop_infect_lo[[i]][min(which(round(prop_infect_lo[[i]][,2],2)==0.80)),1]
  plot80x.hi[i] <- prop_infect_hi[[i]][min(which(round(prop_infect_hi[[i]][,2],2)==0.80)),1]
  plot90x.lo[i] <- prop_infect_lo[[i]][min(which(round(prop_infect_lo[[i]][,2],2)==0.80)),1]
  plot90x.hi[i] <- prop_infect_hi[[i]][min(which(round(prop_infect_hi[[i]][,2],2)==0.80)),1]
  
}

round(cbind(plot50x[dstb], plot50x.lo[dstb], plot50x.hi[dstb])*100)
round(cbind(plot50x[drtb], plot50x.lo[drtb], plot50x.hi[drtb])*100)

round(cbind(plot80x[dstb], plot80x.lo[dstb], plot80x.hi[dstb])*100)
round(cbind(plot80x[drtb], plot80x.lo[drtb], plot80x.hi[drtb])*100)

round(cbind(plot90x[dstb], plot90x.lo[dstb], plot90x.hi[dstb])*100)
round(cbind(plot90x[drtb], plot90x.lo[drtb], plot90x.hi[drtb])*100)

prop_infect_hi[[3]]
studyname
round(plot80x[dstb],2)*100
round(plot80x.lo[dstb],2)*100
round(plot80x.hi[dstb],2)*100

round(plot80x[drtb],2)*100
round(plot80x.lo[drtb],2)*100
round(plot80x.hi[drtb],2)*100

summary(round(plot80x,2)*100)
### Plots ###

library(plotrix)

studyname_char <- strsplit(studyname, "_")
studyname_fortable <- c()
for (i in 1:length(studyname_char)){
  studyname_fortable[i] <- paste0(studyname_char[[i]][1], " (", studyname_char[[i]][2], ")")
}

fifty <- round(cbind(plot50x, plot50x.lo, plot50x.hi)*100)
eighty <- round(cbind(plot80x, plot80x.lo, plot80x.hi)*100)
ninety <- round(cbind(plot90x, plot90x.lo, plot90x.hi)*100)

plottable <- data.frame(Study = studyname_fortable, 
                        color = "",
                        `50%` = paste0(fifty[,1], "% (", fifty[,2],"-",fifty[,3],")"), 
                        `80%` = paste0(eighty[,1], "% (", eighty[,2],"-",eighty[,3],")"), 
                        `90%` = paste0(ninety[,1], "% (", ninety[,2],"-",ninety[,3],")"))
plottable[3,3] <- "1% (1-2)"

colnames(plottable) <- c(NA, "       ", "50%", "80%", "90%")
rownames(plottable) <- plottable[,1]

prop_leg <- c(NA, paste0(dstb_names, "    "), 
         expression(bold("50%")), plottable[dstb,3], 
         expression(bold("80%")), plottable[dstb,4], 
         expression(bold("90%")), plottable[dstb,5])
prop_leg_dr <- c(NA, paste0(drtb_names, "    "), 
              expression(bold("50%")), plottable[drtb,3], 
              expression(bold("80%")), plottable[drtb,4], 
              expression(bold("90%")), plottable[drtb,5])

library(plotrix)

xadj2 <- 0.10
yadj2 <- 0.10




#susceptible

png(paste0(generalpath, "R Figures/Figure_pro_infected_both_", Sys.Date(),".png"), width = 600, height = 1000)
#pdf(paste0(generalpath, "R Figures/Figure_pro_infected_both_", Sys.Date(),".pdf"), width = 6*1.3, height = 10*1.3)
par(mfrow=c(2,1))
#pdf(paste0(generalpath, "R Figures/Figure_prop_A-only_", Sys.Date(),".pdf"), width = 6*1.3, height = 5*1.3)

plot(xx, xx, type = 'n', xlab = '', ylab = '', axes = FALSE)
  for (i in 1:length(prop_infect[dstb])){
    lines(prop_infect[dstb][[i]][,1], prop_infect[dstb][[i]][,2], lty=1, lwd=3, col=colors[dstb][i])
#    polygon(c(xx, rev(xx)), c(prop_infect_lo[dstb][[i]][,2], rev(prop_infect_hi[dstb][[i]][,2])), col = colorst[dstb][i], border = NA)
  }
  abline(h = c(0.50, 0.80, 0.90), lty=3, col="darkgrey")
  axis(side = 1, at=seq(0,1,0.1), labels = seq(0, 100, 10))
  axis(side = 2, at=seq(0,1,0.1))
  text(c("50% of secondary transmission", "80% of secondary transmission", "90% of secondary transmission"),
       x=0.75, y=c(0.50, 0.80, 0.90)+0.02, col = "grey23", cex = 0.8)
  mtext(side=1, 'Percent of Infectious Cases', padj=4)
  mtext(side=2, 'Expected Proportion of Secondary Transmission', padj=-4)
  mtext(expression(bold("A")), adj = 0.02, padj = -0.5, cex = 1.2)
  
  legend(0.3-xadj2, 0.4-yadj2, legend = prop_leg[1:6], bty = "n", cex = 0.8)
  legend(0.6-xadj2, 0.4-yadj2, legend = prop_leg[7:12], bty = "n", cex = 0.8)
  legend(0.75-xadj2, 0.4-yadj2, legend = prop_leg[13:18], bty = "n", cex = 0.8)
  legend(0.92-xadj2, 0.4-yadj2, legend = prop_leg[19:24], bty = "n", cex = 0.8)
  legend(0.52-xadj2, 0.375-yadj2,rep("", times = 5),#c(dstb_names), 
         lty = c(rep(1,length(prop_infect[dstb]))), 
         lwd = 5, col = c(colors[dstb]), bty="n", cex = 0.87)
  #brackets(x1 = 0.99, x2 = 0.99, y1 = 0.3, y2 = 0, type = 1, lwd = 1)
  text(0.85-xadj2,0.45-yadj2,"Percent of infectious cases\nresponsible for secondary transmission, % (UI)",cex = 0.8)

# dev.off()

#resistant
#pdf(paste0(generalpath, "R Figures/Figure_prop_B-only_", Sys.Date(),".pdf"), width = 6*1.3, height = 5*1.3)
plot(xx, xx, type = 'n', xlab = '', ylab = '', axes = FALSE)
  for (i in 1:length(prop_infect[drtb])){
    lines(prop_infect[drtb][[i]][,1], prop_infect[drtb][[i]][,2], lty=1, lwd=3, col=colors[drtb][i])
    #polygon(c(xx, rev(xx)), c(prop_infect_lo[drtb][[i]][,2], rev(prop_infect_hi[drtb][[i]][,2])), col = colorst[drtb][i], border = NA)
    #segments(x0 = plot80x[drtb][i], y0 = 0, y1 = plot80y[drtb][i], lty = 2, col = colors[drtb][i])
  }
  abline(h=c(0.50, 0.80, 0.90), lty=3, col="darkgrey")
  axis(side = 1, at=seq(0,1,0.1), labels = seq(0, 100, 10))
  axis(side = 2, at=seq(0,1,0.1))
  text(c("50% of secondary transmission", "80% of secondary transmission", "90% of secondary transmission"),
       x=0.75, y=c(0.50, 0.80, 0.90)+0.02, col = "grey23", cex = 0.8)
  mtext(side = 1, 'Percent of Infectious Cases', padj=4)
  mtext(side = 2, 'Expected Proportion of Secondary Transmission', padj=-4)
  mtext(expression(bold("B")), adj = 0.02, padj = -0.5, cex = 1.2)
  
  legend(0.3-xadj2, 0.4-yadj2, legend = prop_leg_dr[1:5], bty = "n", cex = 0.8)
  legend(0.6-xadj2, 0.4-yadj2, legend = prop_leg_dr[6:10], bty = "n", cex = 0.8)
  legend(0.75-xadj2, 0.4-yadj2, legend = prop_leg_dr[11:15], bty = "n", cex = 0.8)
  legend(0.92-xadj2, 0.4-yadj2, legend = prop_leg_dr[16:20], bty = "n", cex = 0.8)
  legend(0.52-xadj2, 0.375-yadj2,rep("", times = 4),#c(dstb_names), 
         lty = c(rep(1,length(prop_infect[drtb]))), 
         lwd = 5, col = c(colors[drtb]), bty="n", cex = 0.87)
  
  text(0.85-xadj2,0.45-yadj2,"Percent of infectious cases\nresponsible for secondary transmission, % (UI)",cex = 0.8)
  #segments(x0 = 0.55, x1 = 0.95, y0 = 0.35)
dev.off()






khats <- c()
for (i in 1:length(ests)){
  khats[i] <- ests[[i]][2,1]
}
Rhats <- c()
for (i in 1:length(ests)){
  Rhats[i] <- ests[[i]][1,1]
}



pdf(paste0(generalpath, "R Figures/Figure_prop_bt_paramest_", Sys.Date(),".pdf"), width = 12, height = 8)
par(mfrow=c(1,2))
plot(c(0,1),c(0,1), type='n', ylab = "", xlab = "",axes = FALSE)
    for (i in 1:length(plot80x[drtb])){
      points(plot80x[drtb][i],ests[drtb][[i]][2,1], pch=17)  
    }
    for (i in 1:length(plot80x[dstb])){
      points(plot80x[dstb][i],ests[dstb][[i]][2,1], pch=19)  
    }
  text(0.9,0.1,"r = 0.98\np < 0.001")
  axis(1)
  axis(2)  
  mtext(side=1, 'Proportion of cases responsible for 80% or transmission', padj=4)
  mtext(side=2, 'Paramter estimate (k)', padj=-4)
  legend("topright", c("All TB Surveillance", "Drug Resistant TB Surveillance"),
       pch=c(19,17), bty="n")

plot(c(0,1),c(0,1), type='n', ylab = "", xlab = "",axes = FALSE)
    for (i in 1:length(plot80x[drtb])){
      points(plot80x[drtb][i],ests[drtb][[i]][1,1], pch=2)  
    }
    for (i in 1:length(plot80x[dstb])){
      points(plot80x[dstb][i],ests[dstb][[i]][1,1], pch=1)  
    }
text(0.9,0.1,"r = 0.01\np = 0.98")
  axis(1)
  axis(2)  
  mtext(side=1, 'Proportion of cases responsible for 80% or transmission', padj=4)
  mtext(side=2, 'Paramter estimate (R)', padj=-4)
  legend("topright", c("All TB Surveillance", "Drug Resistant TB Surveillance"),
         pch=c(1,2), bty="n")

dev.off()

cor.test(plot80x,khats)
cor.test(plot80x,Rhats)









### Cluster descriptive data
# median/IQR
sumstats <- list()
for (i in 1:length(allfiles)){
  temp <- read.csv(allfiles[i])
  temp2 <- temp[temp[,1]>1,1]
  sumstats[[i]] <- summary(temp2)
}

sumstats[9]

ga <- read.csv(allfiles[])
i <- 1
num_isolated <- prop_isolated <- c()
for (i in 1:length(allfiles)){
  temp <- read.csv(allfiles[i])
  temp2 <- temp[temp[,1]==1,1]
  num_isolated[i] <- length(temp2)
  prop_isolated[i] <- round(sum(temp2)/sum(temp[,1])*100,0)
}
prop_isolated
prop_clustered <- 100-prop_isolated
summary(prop_clustered)
num_isolated



isolated <- paste0(num_isolated," (", prop_isolated,")")
isolated[drtb]





##### probability of seeing a cluster of size Y
prob <- function(Y,R,k) { 
  p <- exp(lgamma(k*Y+Y-1) - lgamma(k*Y) - lgamma(Y+1) + (Y-1) * log(R/k) - (k*Y+Y-1) * log(1+R/k))
  return(p)
}

Y <- c(5, 10, 25, 50)
Rs <- seq(0.02,1.0,0.005)
ks <-  10^seq(-2.3,0.5,0.005) 
P05 <- P10 <- P25 <- P50 <- matrix(NA, ncol = length(ks), nrow = length(Rs))

# This calculates the probability of observing a cluster of at least size Y for all values across Rs and ks
for(i in 1:length(Rs)) {
  for(j in 1:length(ks)) {
    P05[i,j] <- 1 - sum(prob(1:(Y[1] - 1), Rs[i], ks[j]))
    P10[i,j] <- 1 - sum(prob(1:(Y[2] - 1), Rs[i], ks[j]))
    P25[i,j] <- 1 - sum(prob(1:(Y[3] - 1), Rs[i], ks[j]))
    P50[i,j] <- 1 - sum(prob(1:(Y[4] - 1), Rs[i], ks[j]))
  }
}

lbreaks <- c(-1e-10,1e-5,1e-3,5e-3,1e-2,2e-2,5e-2,1e-1,1.5e-1,1)

image.scale <- function(z, zlim, col = heat.colors(12),breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

pdf(paste0(generalpath, "R Figures/Figure_heatmap_", Sys.Date(),".pdf"), width = 8, height = 8)
#tiff(paste0(generalpath, "R Figures/Figure_heatmap_", Sys.Date(),".tiff"), width = 8, height = 8, units = "in", res = 300)
## use below if legend is requested
# layout(matrix(c(1,2), nrow=1, ncol=2), widths = c(4,1), heights = c(1,1))
# layout.show(4)

par(mfrow = c(2, 2))
image(Rs, ks, P05, xlab = bquote("Reproduction Number," ~ italic("R")), ylab = bquote("Dispersion parameter," ~ italic("k") ~ "(log)"),
      log="y", col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks, zlim = lbreaks)
contour(Rs, ks, P05, levels = lbreaks, labcex = 0.8, add = T)
contour(Rs, ks, P05, levels = lbreaks, labcex = 0.8, add = T)
mtext(paste0("Y = ",Y[1]), adj = 1)
mtext(expression(bold("A")), adj = 0.02, padj = -0.5, cex = 1.2)
#text(x=0.12+0.05, y=0.02-.003, paste0("Chee et al"), col="purple")
#text(x=0.20+0.05, y=0.34-.05, paste0("Yang et al"), col="purple")
for (i in 1:length(ests[dstb])){
  points(ests[dstb][[i]][1,1],ests[dstb][[i]][2,1], pch=19, cex=1.5)
}
for (i in 1:length(ests[drtb])){
  points(ests[drtb][[i]][1,1],ests[drtb][[i]][2,1], pch=2, cex=1.5)
}
#points(x=0.12, y=0.02, pch=24, col="purple")
#points(x=0.20, y=0.34, pch=24, col="purple")
 #legend("topleft", c("All TB Surveillance (MLE)", "Drug Resistant TB Surveillance (MLE)"), pch = c(19, 2), bty="n")


image(Rs, ks, P10, xlab = bquote("Reproduction Number," ~ italic("R")), ylab = bquote("Dispersion parameter," ~ italic("k") ~ "(log)"),
      log="y", col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks, zlim = lbreaks)
  contour(Rs, ks, P10, levels = lbreaks, labcex = 0.8, add = T)
  contour(Rs, ks, P10, levels = lbreaks, labcex = 0.8, add = T)
  mtext(paste0("Y = ", Y[2]), adj = 1)
  mtext(expression(bold("B")), adj = 0.02, padj = -0.5, cex = 1.2)
  #text(x=0.12+0.05, y=0.02-.003, paste0("Chee et al"), col="purple")
  #text(x=0.20+0.05, y=0.34-.05, paste0("Yang et al"), col="purple")
  for (i in 1:length(ests[dstb])){
   points(ests[dstb][[i]][1,1],ests[dstb][[i]][2,1], pch=19, cex=1.5)
  }
  for (i in 1:length(ests[drtb])){
    points(ests[drtb][[i]][1,1],ests[drtb][[i]][2,1], pch=2, cex=1.5)
  }
  #points(x=0.12, y=0.02, pch=24, col="purple")
  #points(x=0.20, y=0.34, pch=24, col="purple")
  #legend("topleft", c("All TB Surveillance (MLE)", "Drug Resistant TB Surveillance (MLE)"), pch = c(19, 2), bty="n")

  
image(Rs, ks, P25, xlab = bquote("Reproduction Number," ~ italic("R")), ylab = bquote("Dispersion parameter," ~ italic("k") ~ "(log)"),
        log="y", col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks, zlim = lbreaks)
  contour(Rs, ks, P25, levels = lbreaks, labcex = 0.8, add = T)
  contour(Rs, ks, P25, levels = lbreaks, labcex = 0.8, add = T)
  mtext(paste0("Y = ", Y[3]), adj = 1)
  mtext(expression(bold("C")), adj = 0.02, padj = -0.5, cex = 1.2)
  #text(x=0.12+0.05, y=0.02-.003, paste0("Chee et al"), col="purple")
  #text(x=0.20+0.05, y=0.34-.05, paste0("Yang et al"), col="purple")
  for (i in 1:length(ests[dstb])){
    points(ests[dstb][[i]][1,1],ests[dstb][[i]][2,1], pch = 19, cex=1.5)
  }
  for (i in 1:length(ests[drtb])){
    points(ests[drtb][[i]][1,1],ests[drtb][[i]][2,1], pch = 2, cex=1.5)
  }
  #points(x=0.12, y=0.02, pch=24, col="purple")
  #points(x=0.20, y=0.34, pch=24, col="purple")
  #legend("topleft", c("All TB Surveillance (MLE)", "Drug Resistant TB Surveillance (MLE)"), pch = c(19, 2), bty="n")

  
image(Rs, ks, P50, xlab = bquote("Reproduction Number," ~ italic("R")), ylab = bquote("Dispersion parameter," ~ italic("k") ~ "(log)"),
        log="y", col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks, zlim = lbreaks)
  contour(Rs, ks, P50, levels = lbreaks, labcex = 0.8, add = T)
  contour(Rs, ks, P50, levels = lbreaks, labcex = 0.8, add = T)
  mtext(paste0("Y = ", Y[4]), adj = 1)
  mtext(expression(bold("D")), adj = 0.02, padj = -0.5, cex = 1.2)
  #text(x=0.12+0.05, y=0.02-.003, paste0("Chee et al"), col="purple")
  #text(x=0.20+0.05, y=0.34-.05, paste0("Yang et al"), col="purple")
  for (i in 1:length(ests[dstb])){
    points(ests[dstb][[i]][1,1],ests[dstb][[i]][2,1], pch=19, cex=1.5)
  }
  for (i in 1:length(ests[drtb])){
    points(ests[drtb][[i]][1,1],ests[drtb][[i]][2,1], pch=2, cex=1.5)
  }
  #points(x=0.12, y=0.02, pch=24, col="purple")
  #points(x=0.20, y=0.34, pch=24, col="purple")
  #legend("topleft", c("All TB Surveillance (MLE)", "Drug Resistant TB Surveillance (MLE)"), pch = c(19, 2), bty="n")
  
  
## Legend color scale
# par(mar=c(4.5,0.5,4,4))
# options(scipen = 8)
# image.scale(log(lbreaks), col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks, horiz = FALSE,  ylim = range(0,0.15), yaxt="n")
# axis(4, at = round(lbreaks,2), las=2)

dev.off()


colnames(P50) <- ks
rownames(P50) <- Rs

Rtemp <- ktemp <-  c()
for (i in 1:length(ests)){
  Rtemp[i] <- ests[[i]][1,1]
  ktemp[i] <- ests[[i]][2,1]
}
temp <- data.frame(studyname=studyname, R=Rtemp, k=ktemp)


sigdig <- 8
for (i in 1:nrow(temp)){
  temp[i,4] <- round(1 - sum(prob(1:(Y[1]-1), temp[i,2], temp[i, 3])), sigdig)
  temp[i,5] <- round(1 - sum(prob(1:(Y[2]-1), temp[i,2], temp[i, 3])), sigdig)
  temp[i,6] <- round(1 - sum(prob(1:(Y[3]-1), temp[i,2], temp[i, 3])), sigdig)
  temp[i,7] <- round(1 - sum(prob(1:(Y[4]-1), temp[i,2], temp[i, 3])), sigdig)
}
names(temp)[4:7] <- c("Y≥5", "Y≥10", "Y≥25", "Y≥50")

write.csv(rbind(temp[dstb,-c(2,3)],temp[drtb,-c(2,3)]), paste0(generalpath, "R Figures/prob_table_", Sys.Date(),".csv"), row.names = FALSE)


0.00001/0.00000000001

cheeP <- 1 - sum(prob(1:(Y[4]-1), 0.12, 0.02))
YangP <- 1 - sum(prob(1:(Y[4]-1), 0.20, 0.34))

cheeP/YangP

### Forrest Plot
textfortable_ds <- allests[allests$tbtype=="ALL",c(1,3)]
textfortable_dr <- allests[allests$tbtype=="DR",c(1,3)]


### Figure

chi95 <- qchisq(0.95, df = 1)/2 
ctlns95 <- list()
for (i in 1:length(allfiles)){
  ctlns95[[i]] <- contourLines(x = Rrange, y = log10(krange), as.matrix(surflikes[[i]] - max(surflikes[[i]])),
                          levels = c(-chi95, chi95+0.001))
}
ctlns95 <- lapply(ctlns95, unlist, recursive = FALSE)

xtick <- seq(0,1.4,0.1)
#ytick <- c(0.01,0.05,0.2,0.1,0.5,1,2)
ytick <- c(0.01,0.05,0.1,0.5,1,5,10,50,100)

colors <- ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)
dstb <- allests$tbtype == "ALL"
drtb <- allests$tbtype == "DR"

pdf(paste0(generalpath, "R Figures/Figure_surface_ests_ds_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange), as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R',ylab="Dispersion parameter, k", 
               xlim=c(0.01,1.15),
               ylim=log10(c(0.02,max(krange))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick);
                 
                points(allests[dstb,"R_hat"], log10(allests[dstb,"k_hat"]), pch = 19, cex = 1, col = colors[dstb])
                 
                for (i in 1:length(ctlns95[dstb])){
                  lines(ctlns95[dstb][[i]][[2]], ctlns95[dstb][[i]][[3]], lwd = 2.5, lty = 1, col=colors[dstb][i])
                }
})
legend(0.55, 0.4, c(allests[dstb,1],NA,"95% CI"), 
       lty = c(rep(1,length(allests[dstb,1])), NA, 1), 
       pch = c(rep(19,length(allests[dstb,1])), NA, NA),
       lwd = 1.5, col = c(colors[dstb], NA, "black"), bty="n")

dev.off()
                 
pdf(paste0(generalpath, "R Figures/Figure_surface_ests_dr_", Sys.Date(),".pdf"), width = 10, height = 8)
filled.contour(x = Rrange, y = log10(krange), as.matrix(maxsurfs[[1]] - max(maxsurfs[[1]])), col='white',levels=seq(-3,0,0.5), 
               xlab = 'Reproduction number, R',ylab="Dispersion parameter, k", 
               xlim=c(0.01,1.15),
               ylim=log10(c(0.02,max(krange))),
               plot.axes = {
                 axis(1, at=xtick, label=xtick);
                 axis(2, at=log10(ytick), label=ytick)
                 
                 points(allests[drtb,"R_hat"], log10(allests[drtb,"k_hat"]), pch = 19, cex = 1, col = colors[drtb])
                 
                 for (i in 1:length(ctlns95[drtb])){
                   lines(ctlns95[drtb][[i]][[2]], ctlns95[drtb][[i]][[3]], lwd = 2.5, lty = 1, col=colors[drtb][i])
                 }
               })
legend(0.55, 0.4, c(allests[drtb,1],NA,"95% CI"), lty = c(rep(1,length(allests[drtb,1])), NA, 1), 
       pch = c(rep(19,length(allests[drtb,1])), NA, NA),lwd = 1.5, col = c(colors[drtb], NA, "black"), bty="n")

dev.off()


 
#  for (i in 1:length(allfiles)){
#    assign(gsub(".*/(.*)\\..*","\\1", allfiles[i], perl = TRUE), ests[[i]])   
#  }
# 


