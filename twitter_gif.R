# rm(list = ls())

filepth <- "~jonathansmith/Downloads/GIF/"
  
propinfection <- function(R, k, prop){
  xp <- qgamma(1 - prop, shape = k, rate = k/R) 
  tp <- 1 - pgamma(xp, shape = k+1, rate = k/R) 
  return(tp) 
}
a <- 1:9
dstb <- c(1,2,4,7:8)
drtb <- setdiff(a, dstb)


ests <- data.frame(study = c("Alvarez (2020)","Asare (2020)",
                             "Chee (2020)",
                             "Guerra-Assunção (2015)",
                             "Jiang (2019)", "Macedo (2020)",
                             "Verza (2020)","Walker (2014)",
                             "Yang (2017)"),
                   R_hat = c(0.73, 0.21, 0.12, 0.48, 0.16, 0.56, 0.29, 0.10, 0.20),
                   k_hat = c(0.07, 0.15, 0.02, 0.24, 0.24, 0.19, 0.48, 0.06, 0.34),
                   tbtype = c("S", "S", "R", "S", "R", "R", "S", "S", "R"))

ests <- ests[order(ests$k_hat),]
estss <- ests[,c(2:3)]
xx <- seq(0,1,0.001)
prop_infect <- list()
# ests[[7]][2,3] <- 10000000 #approximate infinity
for (i in 1:nrow(estss)){
  prop_infect[[i]] <- cbind(xx, 
                            propinfection(estss[i,1], estss[i,2], xx),
                            round(propinfection(estss[i,1], estss[i,2], xx),2))
}
legx <- 0.25
linx <- legx - 0.05
legy <- 0.2
liny <- legy
cexx2 <- 0.80

# output images
for(j in seq(0, 1000, 10)){
  png(paste0(filepth, "gif", j,".png"))
      
  plot(xx, xx, type = 'n', xlab = '', ylab = '', xaxs = "i", yaxs="i", axes = FALSE, xlim = c())
    axis(side = 1, at = seq(0,1,0.1), labels = seq(0, 100, 10))
    axis(side = 2, at = seq(0,1,0.1), las = 1)
    mtext(side = 1, 'Percent of Infectious Cases', padj = 4)
    mtext(side = 2, 'Expected Proportion of Secondary Transmission', padj = -4)
    mtext("Expected Proportion of TB Transmission\nAttributed to a Given Proportion of Infectious Cases",
          padj = -1)
  legend("bottomright", c(ests[,1], NA, "Drug Susceptible TB", "Drug Resistant TB"), 
         lty = c(1,1:3,2,4,3,4,5, NA, 1,1), bty = "n", cex = 1, lwd = c(rep(2, 9), NA, 3,3),
         col = c("maroon", "black", "black", "black",
                 "maroon","black","maroon", "maroon","black", NA, "black", "maroon"))
  ltyy <- c(1,1:3,2,4,3,4,5)
  coll <- c("maroon", "black", "black", "black",
          "maroon","black","maroon", "maroon","black")

for (i in 1:length(prop_infect[dstb])){
  lines(prop_infect[dstb][[i]][1:j,1], prop_infect[dstb][[i]][1:j,2], lty = i, lwd = 2, col = "black")
   }
for (i in 1:length(prop_infect[drtb])){
  lines(prop_infect[drtb][[i]][1:j, 1], prop_infect[drtb][[i]][1:j, 2], lty = i, lwd = 2, col = "maroon")
  }
dev.off()
}


imgs <- list.files(filepth, full.names = TRUE)
rlst <- c(1, 2, 14, 25, 36, 47, 58, 69, 80, 91)
thou <- 4
imgs <- c(imgs[rlst], imgs[setdiff(1:101, c(rlst, thou))], imgs[thou])
img_list <- lapply(imgs, magick::image_read)

## join the images together
img_joined <- magick::image_join(img_list)

## animate at 2 frames per second
img_animated <- magick::image_animate(img_joined, fps = 50)

## view animated image
img_animated

image_write(image = img_animated,
            path = "~jonathansmith/Dropbox/AJE Technical Review/twitter_gif.gif")