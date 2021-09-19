library(gstat)
library(readr)
library(geoR)
library(rstudioapi)

raw.vgm <- read_csv("~/location to the variogram/", 
                    col_names = FALSE)

names(raw.vgm) <- c("np", "dist","gamma","dir.hor","dir.ver")
raw.vgm[] <- lapply(raw.vgm, as.numeric)
class(raw.vgm) = c("gstatVariogram", "data.frame")
plot(raw.vgm)

fit <- fit.variogram(raw.vgm,fit.method=6, vgm(psill = 0.00008,"Gau",nugget=0.00001))
plot(raw.vgm, fit)
theta<- cbind(fit$psill[1],1*fit$psill[2],fit$range[2])

fit2 <- fit.variogram(raw.vgm,fit.method=7 ,vgm("Exp"))
theta <- cbind(fit$psill[2],fit$range[2],fit$psill[1])
plot(raw.vgm, fit)
write.table(theta, file = "theta.csv",sep = ",", row.names = FALSE, col.names=FALSE, qmethod = "double")

