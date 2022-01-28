###########################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Simulation code

## Date created: 12/01/2022
## Last updated: 24/01/2022

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

library(dplyr)
library(INLA) 
library(raster)
library(maptools)
library(ggplot2)
library(gtools)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

set.seed(2022)

## Data
#Read in data file containing coordinates for observation 
#locations
data <- read.csv(paste0(
  "Simulation dataset/",
  "file_name.csv"), 
  head = T) 

data %>% 
    
  ## remove latlon = NA
  filter(LONGNUM != is.na(LONGNUM) |
           LATNUM != is.na(LATNUM)) %>%
  
  ## remove locations where <= 1 child was sampled
  filter(total > 1) -> data

data_coords <- as.matrix(cbind(
  data$LONGNUM, data$LATNUM))

#Example raster data at 1 km
raster <- raster(paste0(
  "Simulation dataset/",
  "n28_viirs_nightlights.tif"))

#Aggregate to 5 km
raster_agg <- aggregate(raster, fact = 5, fun = mean)
raster_val <- as.matrix(getValues(raster_agg))
pred_coords_all <- coordinates(raster_agg)


ind <- apply(raster_val, 1, function(x) any(is.na(x)))
miss <- which(ind == T)
nonmiss <- which(ind == F)

pred_coords <- pred_coords_all[nonmiss, ]


## Meshes
shp_ng  <- readShapePoly(paste0(
  "Simulation dataset/",
  "gadm36_NGA_0.shp"))

shp_df <- fortify(shp_ng)
shp_bnd <- cbind(shp_df$long, shp_df$lat)
c_bnd <- as.matrix(shp_bnd)

# bnd <- inla.nonconvex.hull(c.bnd, convex=-0.1)

meshfit <- inla.mesh.2d(
  loc.domain = c_bnd, 
  max.edge = c(0.2, 0.6),
  cutoff = 0.2)

plot(meshfit)
points(data_coords[,1], data_coords[,2], col = "red")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Priors

## Matern smoothness parameter, redundant 
## as nu = 1 implies alpha = 2.

nu <- 1 
alpha <- 2

## Option 1: r0 = 0.48 # 5% of the extent of 
## Nigeria in the north-south direction (i.e. 
## 0.05*(ymax-ymin) ).

## Option 2: r0 = 2.55 # Q1 of the euclidean
## distance matrix of the data locations.
r0 <- as.numeric(summary(dist(data_coords))[2]) 

sigma0 <- 1


## Matern SPDE model object

spde <- inla.spde2.pcmatern(
  mesh = meshfit, 
  alpha = alpha, 
  prior.range = c(r0, NA), ## fixing spatial range 
  prior.sigma = c(sigma0, NA)) ## fixing spatial SD

kappa0 <- sqrt(8*nu)/r0 ## SPDE scaling parameter
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0) ## SPDE var parameter 

## SPDE precision matrix

Q_stat <- inla.spde2.precision(
  theta = c(log(tau0),log(kappa0)),
  spde = spde)

sample_stat <- as.vector(inla.qsample(
  n = 1, Q = Q_stat, seed = 123))

## A matrix for prediction 

Apred <- inla.spde.make.A(
  loc = pred_coords,
  mesh = meshfit)

Spred <- as.vector(Apred %*% sample_stat) 

## A matrix for observation

Apoints <- inla.spde.make.A(
  loc = data_coords,
  mesh = meshfit)

Spoints <- as.vector(Apoints %*% sample_stat)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Covariate matrix

n <- nrow(data_coords)
m <- nrow(pred_coords)
z <- cbind(1,rnorm(n+m), rgamma(n+m,1,1), rt(n+m,2))

## Regression parameters 
beta <- c(0.5, 0.8, 0.8, 0.2)

## Generate the probabilities

p  <- numeric(n+m)
SS <- c(Spoints, Spred)

## Nugget effect

sigma_e <- 1 
ee <- rnorm(n+m, 0, sigma_e)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## p_1(s)
for (i in 1:(n+m)) p[i] <- inv.logit(
    z[i,1]*beta[1] + z[i,2]*beta[2] + 
    z[i,3]*beta[3] + z[i,4]*beta[4] + 
    SS[i] + ee[i])

p.dat.1  <- p[1:n]
p.pred.1 <- p[(n+1):(n+m)] 

## p_2(s) - note the change on the logit 
## scale by subtracting a factor of 1.3

for (i in 1:(n+m)) p[i] <- inv.logit(
    z[i,1]*beta[1] + z[i,2]*beta[2] + 
    z[i,3]*beta[3] + z[i,4]*beta[4] + 
    SS[i] + ee[i] - 1.3) 

p.dat.2  <- p[1:n]
p.pred.2 <- p[(n+1):(n+m)] 

summary(p.pred.1-p.pred.2)
summary(p.pred.2)

## p_3(s) - note the change on the logit 
## scale by subtracting a factor of 2.5

for (i in 1:(n+m)) p[i] <- inv.logit(
    z[i,1]*beta[1] + z[i,2]*beta[2] + 
    z[i,3]*beta[3] + z[i,4]*beta[4] + 
    SS[i] + ee[i] - 2.5)

p.dat.3  <- p[1:n]
p.pred.3 <- p[(n+1):(n+m)] 

summary(p.pred.2-p.pred.3)
summary(p.pred.3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Reconstruct simulated surfaces and save these

## p_1(s)

ll <- 1:length(ind)
ll[nonmiss] <- p.pred.1
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "prob_p1_sim.tif", overwrite = TRUE)

## p_2(s)

ll <- 1:length(ind)
ll[nonmiss] <- p.pred.2
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "prob_p2_sim.tif", overwrite = TRUE)

## p_3(s)

ll <- 1:length(ind)
ll[nonmiss] <- p.pred.3
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "prob_p3_sim.tif", overwrite = TRUE)


## Construct surfaces for the simulated covariates

## X2

ll <- 1:length(ind)
ll[nonmiss] <- z[(n+1):(n+m),2]
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "sim_cov_x2.tif", overwrite = TRUE)

## X3

ll <- 1:length(ind)
ll[nonmiss] <- z[(n+1):(n+m),3]
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "sim_cov_x3.tif", overwrite = TRUE)

## X4

ll <- 1:length(ind)
ll[nonmiss] <- z[(n+1):(n+m),4]
ll[miss] <- NA

rr.out <- raster(raster_agg)
values(rr.out) <- ll
plot(rr.out)

writeRaster(rr.out, "sim_cov_x4.tif", overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Different sample sizes

## data_coords is in lon/lat
data.frame(cbind(z[1:n,], data_coords)) %>%
  rename(LONGNUM = X5, LATNUM = X6) -> cov_dat

write.csv(cov_dat, "Sim_cov_dat.csv")

## 2-10 (discrete uniform distribution)
n.dat <- sample(2:10, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_10.csv")

#2-20 (discrete uniform distribution)
n.dat <- sample(2:20, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_20.csv")

#2-30 (discrete uniform distribution)
n.dat <- sample(2:30, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_30.csv")

#2-40 (discrete uniform distribution)
n.dat <- sample(2:40, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_40.csv")

#2-50 (discrete uniform distribution)
n.dat <- sample(2:50, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_50.csv")


#2-60 (discrete uniform distribution)
n.dat <- sample(2:60, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_60.csv")

#2-70 (discrete uniform distribution)
n.dat <- sample(2:70, n, replace = T)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_70.csv")

#2-80 (discrete uniform distribution)
n.dat <- sample(2:80, n, replace = TRUE)
y.dat.1 <- round(n.dat*p.dat.1, 0)
y.dat.2 <- round(n.dat*p.dat.2, 0)
y.dat.3 <- round(n.dat*p.dat.3, 0)

dd <- data.frame(n.dat, y.dat.1, p.dat.1, y.dat.2, 
                 p.dat.2, y.dat.3, p.dat.3)

head(dd)
write.csv(dd, "Sim_dat_2_80.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
###########################################################
