# R code for "Amazonian hyperdominance exceeds expectations of neutral theory"
# Uses publicly available ATDN data from ter Steege et al 2020 (file = "terSteege_2019_bias_corrected_data.csv") 
    # https://www.nature.com/articles/s41598-020-66686-3
# All simulation code found here: https://github.com/tfmilton/Hyperdominance

library (sads)
library (vegan)
library (untb)
library (Rmisc)
library (Hmisc)
library (dgof)
library (KSgeneral)
library (rgr)
library (breheny)
library (bitops)
library (RCurl)
library (breheny)
library (SADISA)
library (viridis)

########## Rand.neutral simulations and analysis ##########

ATDN <- read.csv ("terSteege_2019_bias_corrected_data.csv", stringsAsFactors = F)
 
ATDN_tots <- ATDN$n.ind #make a vector
names (ATDN_tots) <- ATDN$Accepted.species
ATDN_tots <- ATDN_tots [order(ATDN_tots, decreasing = T)]
ATDN_census <- census(ATDN_tots) #Make census
optimal.theta(ATDN_census) # 692.8512
no.of.ind(ATDN_census) # 979614
no.of.spp(ATDN_census) # 5027
no.of.singletons(ATDN_census) # 506
# To compare theta from fitting mzsm
mz_dist <- fitmzsm(x = ATDN_tots)
coef (mz_dist) # 692.2475

# Use rand.neutral to simulate neutral communities 
my.theta <- optimal.theta(ATDN_census)
totinds <- no.of.ind(ATDN_census)
nsims <- 1000
rand.neut_sim <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim <- append (rand.neut_sim, a)
  #cat (ii)
}
saveRDS (rand.neut_sim, file = "rand.neut_sim_TS_2020.rds")

rand.neut_sim2 <- rand.neut_sim
## Make all list elements the same length.
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
## Convert the list into a dataframe.
rand.neut.TS_2020_DF <- matrix(unlist(rand.neut_sim2),
                       ncol = max.len,
                       nrow = nsims,
                       byrow = T )
rand.neut.TS_2020_DF[is.na(rand.neut.TS_2020_DF) ] <- 0
colnames(rand.neut.TS_2020_DF) <- 1:ncol(rand.neut.TS_2020_DF)
write.csv(rand.neut.TS_2020_DF) #save

#
RN_ATDN <- read.csv("rand.neut.TS_2020_DF.csv", stringsAsFactors = F)
max_abund <- max (RN_ATDN[,1])

## Find simulation envelope
CI_025 <- numeric (length (ncol (RN_ATDN)))
CI_975 <- numeric (length (ncol (RN_ATDN)))
mean_RN <- numeric (length (ncol (RN_ATDN)))
for (ii in 1:ncol (RN_ATDN)){
  temp_vec <- RN_ATDN [, ii]
  CI_025 [ii] <- quantile (temp_vec, 0.025)
  CI_975 [ii] <- quantile (temp_vec, 0.975)
  mean_RN [ii] <- mean (temp_vec, na.rm = T)
}

# theta uncertainty
optimal.theta(ATDN_census, like = 2) #2 likelihood units (recommended)
#lower      mle    upper 
#672.0320 692.8512 714.1306

### High theta
my.theta <- 714.1306
totinds <- no.of.ind(ATDN_census)
nsims <- 100
rand.neut_sim <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim <- append (rand.neut_sim, a)
  cat (ii)
}
saveRDS (rand.neut_sim)

rand.neut_sim2 <- rand.neut_sim
## Make all list elements the same length.
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
## Convert the list into a dataframe.
rand_neut_hightheta <- matrix(unlist(rand.neut_sim2),
                              ncol = max.len,
                              nrow = nsims,
                              byrow = T )
rand_neut_hightheta[is.na(rand_neut_hightheta) ] <- 0
colnames(rand_neut_hightheta) <- 1:ncol(rand_neut_hightheta)

write.csv(rand_neut_hightheta)
rand_neut_hightheta <- read.csv ("rand_neut_hightheta_DF.csv")

### Low theta
my.theta <- 672.0320
totinds <- no.of.ind(ATDN_census)
nsims <- 100
rand.neut_sim_low <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim_low <- append (rand.neut_sim_low, a)
  cat (ii, "\t")
}

saveRDS (rand.neut_sim_low)

rand.neut_sim2 <- rand.neut_sim_low
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
rand_neut_lowtheta <- matrix(unlist(rand.neut_sim2),
                             ncol = max.len,
                             nrow = nsims,
                             byrow = T )
rand_neut_lowtheta[is.na(rand_neut_lowtheta) ] <- 0
colnames(rand_neut_lowtheta) <- 1:ncol(rand_neut_lowtheta)

write.csv(rand_neut_lowtheta)

rand_neut_lowtheta <- read.csv ("rand_neut_lowtheta_DF.csv")

### optimal theta
CI_025_opt <- numeric (length (ncol (rn_full_ATDN)))
CI_975_opt <- numeric (length (ncol (rn_full_ATDN)))
mean_opt <- numeric (length (ncol (rn_full_ATDN)))
for (ii in 1:ncol (rn_full_ATDN)){
  temp_vec <- rn_full_ATDN [, ii]
  CI_025_opt [ii] <- quantile (temp_vec, 0.025)
  CI_975_opt [ii] <- quantile (temp_vec, 0.975)
  mean_opt [ii] <- mean (temp_vec)
}

### high theta
rand_neut_hightheta <- as.data.frame (rand_neut_hightheta)
CI_025_high <- numeric (length (ncol (rand_neut_hightheta)))
CI_975_high <- numeric (length (ncol (rand_neut_hightheta)))
mean_high <- numeric (length (ncol (rand_neut_hightheta)))
for (ii in 1:ncol (rand_neut_hightheta)){
  temp_vec <- rand_neut_hightheta [, ii]
  CI_025_high [ii] <- quantile (temp_vec, 0.025)
  CI_975_high [ii] <- quantile (temp_vec, 0.975)
  mean_high [ii] <- mean (temp_vec)
}

### Low theta
rand_neut_lowtheta <- as.data.frame (rand_neut_lowtheta)
CI_025_low <- numeric (length (ncol (rand_neut_lowtheta)))
CI_975_low <- numeric (length (ncol (rand_neut_lowtheta)))
mean_low <- numeric (length (ncol (rand_neut_lowtheta)))
for (ii in 1:ncol (rand_neut_lowtheta)){
  temp_vec <- rand_neut_lowtheta [, ii]
  CI_025_low [ii] <- quantile (temp_vec, 0.025)
  CI_975_low [ii] <- quantile (temp_vec, 0.975)
  mean_low [ii] <- mean (temp_vec)
}

sim_CI_RAD_thetaRange <- data.frame (empir = ATDN_tots, mean_opt = mean_opt [1:length (ATDN_tots)], 
                                     CI_025_opt = CI_025_opt[1:length (ATDN_tots)],
                                     CI_975_opt = CI_975_opt [1:length (ATDN_tots)], 
                                     mean_lowT = mean_low [1:length (ATDN_tots)], CI_025_lowT = CI_025_low[1:length (ATDN_tots)],
                                     CI_975_lowT = CI_975_low [1:length (ATDN_tots)],
                                     mean_highT = mean_high [1:length (ATDN_tots)], CI_025_highT = CI_025_high[1:length (ATDN_tots)],
                                     CI_975_highT = CI_975_high [1:length (ATDN_tots)],
                                     row.names = 1:length (ATDN_tots))

# test if each rank is above, within, or below neutral simulation envelope
comp.each.rank <- numeric (length (ATDN_tots))
comp.value <- numeric (length (ATDN_tots))
for (ii in 1:length (ATDN_tots)){
  temp_rank <- ATDN_tots[ii]
  if (temp_rank < CI_025_opt[ii]){
    comp.each.rank[ii] <- "below"
    comp.value [ii] <- temp_rank - CI_025_opt[ii]
  }else if (temp_rank > CI_975_opt [ii]){
    comp.each.rank [ii] <- "above"
    comp.value [ii] <- temp_rank - CI_975_opt[ii]
  }else{
    comp.each.rank [ii] <- "in"
    comp.value [ii] <- 0
  }
}
sim_CI_RAD_thetaRange$comp_opt <- comp.each.rank
sim_CI_RAD_thetaRange$col <- ifelse(sim_CI_RAD_thetaRange$comp_opt == "above", "#FF0000", 
                                    ifelse(sim_CI_RAD_thetaRange$comp_opt == "below", "#0000FF", "#000000"))

# test each rank with expanded SE including theta range
comp.each.rank_Trange <- numeric (length (ATDN_tots))
for (ii in 1:length (ATDN_tots)){
  temp_rank <- ATDN_tots[ii]
  if (temp_rank < CI_025_high[ii] & temp_rank < CI_025_low [ii]){
    comp.each.rank_Trange[ii] <- "below"
  }else if (temp_rank > CI_975_low [ii] & temp_rank > CI_975_high [ii]){
    comp.each.rank_Trange [ii] <- "above"
  }else{
    comp.each.rank_Trange [ii] <- "in"
  }
}
sim_CI_RAD_thetaRange$comp_Trange <- comp.each.rank_Trange
sim_CI_RAD_thetaRange$col_Trange <- ifelse(sim_CI_RAD_thetaRange$comp_Trange == "above", "#FF0000", 
                                           ifelse(sim_CI_RAD_thetaRange$comp_Trange == "below", "#0000FF", "#000000"))

### Figure 1 ###
windows (4, 4)
par(mar = c(4,4,1,1))
plot (sim_CI_RAD_thetaRange$mean_opt, type = "l", xlab = "Rank", ylab = "Abundance", lwd = 2,
      ylim = c(1, max(ATDN_tots)), log = "y")
polygon(x = c(1:4902, rev(1:4902)),
        y = c(sim_CI_RAD_thetaRange$CI_975_opt [1:4902], rev(sim_CI_RAD_thetaRange$CI_025_opt [1:4902])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
points (sim_CI_RAD_thetaRange$empir, 
        col = ifelse (sim_CI_RAD_thetaRange$comp_opt == "above", "red", ifelse(sim_CI_RAD_thetaRange$comp_opt == "below", "blue","black")), pch = 20) 
legend ("bottomleft", legend = c("neutral sims+SE" ,"above neutral", "within neutral", "below neutral"), 
        col = c("grey", "red", "black", "blue"), lwd = c(7, 2, 2, 2), cex = 0.7, bty = "n")
legend ("bottomleft", legend = c(NA, NA, NA, NA), lwd = c(2, NA, NA, NA), cex = 0.7, bty = "n")
#Inset
par (cex = 0.6,lwd = 0.5, mar = c(1,3,1.6,1.6), fig = c(0.6,0.999,0.6,0.999), new = T) #oma=c(3.2,2,1,1),
plot (sim_CI_RAD_thetaRange$mean_opt, type = "l", xlab = "Rank", ylab = "Abundance", lwd = 2,
      ylim = c(1, 18000), xlim = c(0, 50), cex.axis = 1.2)
polygon(x = c(1:50, rev(1:50)),
        y = c(sim_CI_RAD_thetaRange$CI_975_opt [1:50], rev(sim_CI_RAD_thetaRange$CI_025_opt [1:50])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
points (sim_CI_RAD_thetaRange$empir, 
        col = ifelse (sim_CI_RAD_thetaRange$comp_opt == "above", "red", ifelse(sim_CI_RAD_thetaRange$comp_opt == "below", "blue","black")), pch = 19) 


### Figure S1 -- with theta range ###
windows (4, 4)
par(mar = c(4,4,1,1))
plot (sim_CI_RAD_thetaRange$mean_opt, type = "l", xlab = "Rank", ylab = "Abundance", lwd = 2,
      ylim = c(1, max(ATDN_tots)), log = "y")
polygon(x = c(1:4777, rev(1:4777)),
        y = c(sim_CI_RAD_thetaRange$CI_975_lowT [1:4777], rev(sim_CI_RAD_thetaRange$CI_025_lowT [1:4777])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
polygon(x = c(1:4777, rev(1:4777)),
        y = c(sim_CI_RAD_thetaRange$CI_975_highT [1:4777], rev(sim_CI_RAD_thetaRange$CI_025_highT [1:4777])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
points (sim_CI_RAD_thetaRange$empir, 
        col = ifelse (sim_CI_RAD_thetaRange$comp_Trange == "above", "red", ifelse(sim_CI_RAD_thetaRange$comp_Trange == "below", "blue","black")), pch = 20) 
legend ("bottomleft", legend = c("neutral sims+SE" ,"above neutral", "within neutral", "below neutral"), 
        col = c("grey", "red", "black", "blue"), lwd = c(7, 2, 2, 2), cex = 0.7, bty = "n")
legend ("bottomleft", legend = c(NA, NA, NA, NA), lwd = c(2, NA, NA, NA), cex = 0.7, bty = "n")
#Inset
par (cex = 0.6,lwd = 0.5, mar = c(1,3,1.6,1.6), fig = c(0.6,0.999,0.6,0.999), new = T) #oma=c(3.2,2,1,1),
plot (sim_CI_RAD_thetaRange$mean_opt, type = "l", xlab = "Rank", ylab = "Abundance", lwd = 2,
      ylim = c(1, 18000), xlim = c(0, 50), cex.axis = 1.2)
polygon(x = c(1:50, rev(1:50)),
        y = c(sim_CI_RAD_thetaRange$CI_975_lowT [1:50], rev(sim_CI_RAD_thetaRange$CI_025_lowT [1:50])),
        col =  adjustcolor("black", alpha.f = 0.15), border = NA)
polygon(x = c(1:50, rev(1:50)),
        y = c(sim_CI_RAD_thetaRange$CI_975_highT [1:50], rev(sim_CI_RAD_thetaRange$CI_025_highT [1:50])),
        col =  adjustcolor("black", alpha.f = 0.15), border = NA)
points (sim_CI_RAD_thetaRange$empir, 
        col = ifelse (sim_CI_RAD_thetaRange$comp_Trange == "above", "red", ifelse(sim_CI_RAD_thetaRange$comp_Trange == "below", "blue","black")), pch = 19)


## find proportion of hyperdominant species, diversity indices, Preston bins 

# hyperdominant species
sims <- nrow (RN_ATDN)
cumu50 <- numeric (length (sims))
cumu50_prop <- numeric (length (sims))
for (ii in 1:sims){
  temp_sim <- as.numeric (RN_ATDN [ii, ])
  temp_sim <- temp_sim [temp_sim >0]
  cumu50 [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  cumu50_prop [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1) / length (temp_sim)
}
hist (cumu50, breaks = 100, xlab = "# species", xlim = c(210, 325),
      main = "species accumulating 50% of inds")
lines (x = c(((length(cumsum(ATDN_tots)[cumsum(ATDN_tots) < sum(ATDN_tots)/2])) + 1),
             (length(cumsum(ATDN_tots)[cumsum(ATDN_tots) < sum(ATDN_tots)/2])) + 1),
       y = c(0, 100), col = "red", lwd = 3)
range (cumu50)   #220 to 306
range (cumu50_prop) #0.04431910 to 0.05952149
# the number of species that make up 50% of inds in ATDN data is less than the lowest number in the simulated SADs 

## Diversity indices
div_shan <- diversity (ATDN_tots, index = "shannon")
div_simp <- diversity (ATDN_tots, index = "simpson")
div_invsimp <- diversity (ATDN_tots, index = "invsimpson")

# Find the range of shannon indices for simulated data

sim_shan <- numeric (length (sims))
sim_simp <- numeric (length (sims))
sim_invsimp <- numeric (length (sims))
for (ii in 1:sims){
  temp_sim <- as.numeric (RN_ATDN[ii,])
  sim_shan [ii] <- diversity (temp_sim, index = "shannon")
  sim_simp [ii] <- diversity (temp_sim, index = "simpson")
  sim_invsimp[ii] <- diversity(temp_sim, index = "invsimpson")
}

## Figure S7
#Shannon
windows (6.8, 2)
par (mfrow = c(1, 3))
par (mar = c(4,4,2,3))
hist (sim_shan, breaks = 50, xlim = c(6.95, max (sim_shan)), main = "Shannon", xlab = " ")
lines(x = c(div_shan, div_shan), y = c(0, 200), lwd=3, col= "red")
#Simpson
hist (sim_simp, breaks = 30, xlim = c(.9975, max(sim_simp)), main = "Simpson", xlab = " ")
lines(x = c(div_simp, div_simp), y = c(0, 200), lwd=3, col= "red")
#Inverse Simpson
hist (sim_invsimp, breaks = 50, xlim = c(400, max(sim_invsimp)), main = "Inv Simpson", xlab = " ")
lines(x = c(div_invsimp, div_invsimp), y = c(0, 200), lwd=3, col= "red")


## Species richness
SR <- numeric (nrow (RN_ATDN))
for (ii in 1:nrow (RN_ATDN)){
  temp_vec <- as.vector (RN_ATDN [ii, ])
  SR [ii] <- length (temp_vec [temp_vec > 0])
}

## Preston octaves
octav(ATDN_tots)
ATDN.OC <- preston (untb::count(ATDN_tots))
rn_octs <- matrix (data = NA, nrow = nrow (RN_ATDN), ncol = 15)
colnames(rn_octs) <- names (preston (ATDN_tots, n = 15))
for (ii in 1:nrow (RN_ATDN)){
  rn_octs [ii, ] <- preston (RN_ATDN [ii, ], n = 15)
}

# mean and SE for all octaves
mean_rn <- numeric (15)
CI_upper <- numeric (15)
CI_lower <- numeric (15)
for (ii in 1:15){
  temp <- rn_octs[,ii]
  mean_rn [ii] <- mean (temp)
  CI_upper [ii] <- quantile (temp, .975)[[1]]
  CI_lower [ii] <- quantile (temp, .025)[[1]]
}
oct.mean.CI.DF <- data.frame (bin = names (preston (ATDN_tots, n = 15)), 
                              mean = mean_rn, CI_upper = CI_upper, CI_lower = CI_lower,
                              ATDN = as.vector (ATDN.OC))

## Figure 3
windows (6.8, 4)
par (mfrow = c(3, 5), omi = c(.2, .2, .2, 0))
par (mar = c(2.5, 2.5, 2.5, 1))
layout.show(15)
hist (rn_octs[,1], xlim = c(506, max (rn_octs[,1])), main = "1", xlab = "", ylab = "")
lines (x = c(506, 506), y = c(0, 300), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[1], oct.mean.CI.DF$CI_lower [1]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[1], oct.mean.CI.DF$CI_upper [1]), y = c(0, 300), lwd = 3)

hist (rn_octs[,2], main = "2", xlab = "",  ylab = "")
lines (x = c(322, 322), y = c(0, 300), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[2], oct.mean.CI.DF$CI_lower [2]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[2], oct.mean.CI.DF$CI_upper [2]), y = c(0, 300), lwd = 3)

hist (rn_octs[,3], xlim = c(350, max(rn_octs[,3])), main = "3-4", xlab = "", ylab = "")
lines (x = c(353, 353), y = c(0, 200), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[3], oct.mean.CI.DF$CI_lower [3]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[3], oct.mean.CI.DF$CI_upper [3]), y = c(0, 300), lwd = 3)

hist (rn_octs[,4], main = "5-8", xlab = "",  ylab = "")
lines (x = c(434, 434), y = c(0, 200), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[4], oct.mean.CI.DF$CI_lower [4]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[4], oct.mean.CI.DF$CI_upper [4]), y = c(0, 300), lwd = 3)

hist (rn_octs[,5], main = "9-16", xlab = "",  ylab = "")
lines (x = c(488, 488), y = c(0, 200), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[5], oct.mean.CI.DF$CI_lower [5]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[5], oct.mean.CI.DF$CI_upper [5]), y = c(0, 300), lwd = 3)

hist (rn_octs[,6], xlim = c(min (rn_octs[,6]), 600), main = "17-32", xlab = "",  ylab = "")
lines (x = c(579, 579), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[6], oct.mean.CI.DF$CI_lower [6]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[6], oct.mean.CI.DF$CI_upper [6]), y = c(0, 300), lwd = 3)

hist (rn_octs[,7], xlim = c(min (rn_octs[,7]), 600), main = "33-64", xlab = "",  ylab = "")
lines (x = c(574, 574), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[7], oct.mean.CI.DF$CI_lower [7]), y = c(0, 400), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[7], oct.mean.CI.DF$CI_upper [7]), y = c(0, 400), lwd = 3)

hist (rn_octs[,8], main = "65-128", xlab = "",  ylab = "")
lines (x = c(483, 483), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[8], oct.mean.CI.DF$CI_lower [8]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[8], oct.mean.CI.DF$CI_upper [8]), y = c(0, 300), lwd = 3)

hist (rn_octs[,9], main = "129-256", xlab = "",  ylab = "")
lines (x = c(449, 449), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[9], oct.mean.CI.DF$CI_lower [9]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[9], oct.mean.CI.DF$CI_upper [9]), y = c(0, 300), lwd = 3)

hist (rn_octs[,10], main = "257-512", xlab = "",  ylab = "")
lines (x = c(374, 374), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[10], oct.mean.CI.DF$CI_lower [10]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[10], oct.mean.CI.DF$CI_upper [10]), y = c(0, 300), lwd = 3)

hist (rn_octs[,11], main = "513-1024            ", xlab = "",  ylab = "")
lines (x = c(254, 254), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[11], oct.mean.CI.DF$CI_lower [11]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[11], oct.mean.CI.DF$CI_upper [11]), y = c(0, 300), lwd = 3)

hist (rn_octs[,12], xlim = c(125, max (rn_octs[,12])), main = "1025-2048", xlab = "",  ylab = "")
lines (x = c(135, 135), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[12], oct.mean.CI.DF$CI_lower [12]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[12], oct.mean.CI.DF$CI_upper [12]), y = c(0, 300), lwd = 3)

hist (rn_octs[,13], main = "2049-4096", xlab = "",  ylab = "")
lines (x = c(52, 52), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[13], oct.mean.CI.DF$CI_lower [13]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[13], oct.mean.CI.DF$CI_upper [13]), y = c(0, 300), lwd = 3)

hist (rn_octs[,14], main = "4097-8192", xlab = "",  ylab = "")
lines (x = c(14, 14), y = c(0, 400), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[14], oct.mean.CI.DF$CI_lower [14]), y = c(0, 300), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[14], oct.mean.CI.DF$CI_upper [14]), y = c(0, 300), lwd = 3)

hist (rn_octs[,15], xlim = c(0, 6), main = "8193-Inf", xlab = "",  ylab = "")
lines (x = c(6, 6), y = c(0, 1000), col = "red", lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_lower[15], oct.mean.CI.DF$CI_lower [15]), y = c(0, 800), lwd = 3)
lines (x = c(oct.mean.CI.DF$CI_upper[15], oct.mean.CI.DF$CI_upper [15]), y = c(0, 800), lwd = 3)

mtext ("Frequency", side = 2, outer = T, at = 0.5)
mtext ("N species", side = 1, outer = T, at = 0.5)
mtext ("Preston Octave Bins", side = 3, outer = T, at = 0.5)


## Full RAD Kolmogorav-Smirnov comparison, max distance

#ecdf for each simulation
maxA <- max(RN_ATDN$X1)
output_ecdfs <- matrix (NA, nrow (RN_ATDN), maxA)
for (ii in 1:nrow (RN_ATDN)){
  temp <- as.vector (t(RN_ATDN [ii,]))
  tempecdf <- ecdf (temp)
  output_ecdfs [ii, ] <- tempecdf (1:maxA)
}
mean_ecdf <- colMeans(output_ecdfs)
plot (mean_ecdf)

mean4step <- c(0, mean_ecdf)
meanRN_ecdf <- stepfun(x = 1:14444, y = mean4step)
plot (meanRN_ecdf)
head (mean_ecdf)
dgof::ks.test(ATDN_tots, meanRN_ecdf)
#D = 0.092383

## Find  ecdf SE
CIL_ECDF <- numeric (ncol (output_ecdfs))
CIU_ECDF <- numeric (ncol (output_ecdfs))
for (ii in 1:ncol (output_ecdfs)){
  CIL_ECDF [ii] <- quantile (output_ecdfs [, ii], 0.025)
  CIU_ECDF [ii] <- quantile (output_ecdfs [, ii], 0.975)
}

CIU4step <- c(0, CIU_ECDF)
RNCIU_ecdf <- stepfun(x = 1:14444, y = CIU4step)

CIL4step <- c(0, CIL_ECDF)
RNCIL_ecdf <- stepfun(x = 1:14444, y = CIL4step)

## max difference
outRN <- meanRN_ecdf(1:14444)
outATDN <- ATDNecdf(1:14444)
absdifs <- abs(outRN - outATDN)
maxdif <- max(absdifs)
#0.09238295  - same as D-stat

# max dif position
maxPOS <- which.max(absdifs)
# 6

### Compare each RN sample to mean RN ecdf to get distribution of K-S D statistics 
KS_D_stat <- numeric (nrow (RN_ATDN))
for (ii in 1:nrow (RN_ATDN)){
  temp_vec <- as.vector(t(RN_ATDN [ii,]))
  temp_ks <- dgof::ks.test (temp_vec, meanRN_ecdf)
  KS_D_stat [ii] <- temp_ks$statistic
}
rnd <- dgof::ks.test (ATDN_tots, meanRN_ecdf)
quantile (KS_D_stat, c(.025, .975))
# 2.5%       97.5% 
# 0.00589216 0.02621688
range (KS_D_stat)
# 0.004324492 0.039571294


## Figure S3
windows (6.8, 3)
par (mfrow = c(1, 2))
par (mar = c(4, 4, 1, 2))
plot(ATDNecdf, verticals=T, log = "x", do.points = F, ylab = "Cumulative Rel. Abund.",
     las=1, col="red", lwd=2, ylim = c(0, 1), xlab = "Abundance", main = "", xlim = c(1, 16101))
plot(meanRN_ecdf, verticals=T, do.points=F, cex.lab=1.2,
     cex.axis=1.3, col="black", lwd=2, ylim = c(0, 1), add = T)
polygon.step(0:14444, CIL4step [1:14444], CIU4step [1:14444], col = "gray") 
plot(meanRN_ecdf, verticals=T, do.points=F, cex.lab=1.2,
     cex.axis=1.3, col="black", lwd=2, ylim = c(0, 1), add = T)
plot (ATDNecdf, verticals = T, do.points = F, col = "red", add = T)
lines (x = c(maxPOS, maxPOS), y = c(0, 1), lty = 3) #lwd = 2)
legend(200, 0.3, legend=c("ATDN", "mean_sim", "maxDif"), col=c("red", "black", "black"), 
       lty=c(1,1,3), bty="n", cex = 0.7)

hist (KS_D_stat, breaks = 50, xlim = c(.003, 0.1), ylim = c(0, 110), xlab = "D statistic", main = " ")
lines (x = c(maxdif, maxdif), y = c(0, 100), col = "red", lwd = 3)
lines (x = c(quantile (KS_D_stat, .025), quantile (KS_D_stat, .025)), y = c(0, 100), lwd = 3)
lines (x = c(quantile (KS_D_stat, .975), quantile (KS_D_stat, .975)), y = c(0, 100), lwd = 3)

## chi-square tests
ATDN.OC <- preston (untb::count(ATDN_tots))
rn_octs <- matrix (data = NA, nrow = nrow (RN_ATDN), ncol = 14)
colnames(rn_octs) <- names (preston (ATDN_tots, n = 14))
for (ii in 1:nrow (RN_ATDN)){
  rn_octs [ii, ] <- preston (RN_ATDN [ii, ], n = 14)
}
mean_prest_rn <- colMeans (rn_octs)

# Compare each RN SAD sample to mean RN SAD
chi_stats <- numeric (length (nrow (RN_ATDN)))
p_values <- numeric (length (nrow (RN_ATDN)))
expt <- 0
resids <- 0
stdres <- 0
for (ii in 1:nrow (rn_octs)){
  tempchi <- chisq.test(x = rn_octs [ii, ], p = mean_prest_rn, rescale.p = T)
  chi_stats [ii] <- tempchi$statistic
  p_values [ii] <- tempchi$p.value
  tempexp <- tempchi$expected
  expt <- list (expt, tempexp)
  tempresids <- tempchi$residuals
  resids <- list (resids, tempresids)
  tempstdres <- tempchi$stdres
  stdres <- list (stdres, tempstdres)
}
hist (chi_stats)

## Compare empirical to mean simulated Preston
data_preston <- preston (ATDN_census, n = 14)
data_chisq <- chisq.test(x = data_preston, p = mean_prest_rn, rescale.p = T, simulate.p.value = T)

#Figure 
windows (4, 3)
par (mar = c(4, 4, 1, 1))
hist (chi_stats, breaks = 20, xlim = c(0, 150), main = " ", xlab = "X-squared")
lines (x = c(data_chisq$statistic, data_chisq$statistic), y = c(0, 250), col = "red", lwd = 3)


########## Spatial neutral model analysis ##########
# Random and aggregated sampling from spatially-explicit neutral communities 
# simulated using C code on github: https://github.com/tfmilton/Hyperdominance


### Randomly sampled metacommunity 

## Global dispersal
sampRADs_glob_rand <- matrix (NA, 500, 200) # put the RADs for each sample
SRs_glob_rand <- numeric (500) # SR of each sample
for (ii in 1:500){
  lscp <- read.table (paste(file = "GL_compAMlscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  lscp <- as.matrix (lscp)
  samp <- sample (lscp, 10000, replace = F)
  tc <- as.data.frame (as.count (samp))
  sampRADs_glob_rand [ii, 1:nrow (tc)] <- tc$Freq
  SRs_glob_rand [ii] <- nrow (tc)
}

# order
for (ii in 1:nrow (sampRADs_glob_rand)){
  sampRADs_glob_rand [ii, ] <- sort  (sampRADs_glob_rand [ii, ], decreasing = T, na.last = T)
}

max (SRs_glob_rand)  
sampRADs_glob_rand <- sampRADs_glob_rand [, 1:54] 
sampRADs_glob_rand [is.na (sampRADs_glob_rand)] <- 0

Mean_rad_glob_rand <- numeric (ncol (sampRADs_glob_rand))
SE_hi_rad_glob_rand <- numeric (ncol (sampRADs_glob_rand))
SE_lo_rad_glob_rand <- numeric (ncol (sampRADs_glob_rand))
for (ii in 1:ncol (sampRADs_glob_rand)){
  Mean_rad_glob_rand [ii] <- mean (sampRADs_glob_rand [, ii])
  SE_hi_rad_glob_rand [ii] <- STDERR (sampRADs_glob_rand [, ii])[[1]]
  SE_lo_rad_glob_rand [ii] <- STDERR (sampRADs_glob_rand [, ii])[[3]]
}


## Dispersal limited, sigma = 5
sampRADs_gaus5_rand <- matrix (NA, 500, 200) 
SRs_gaus5_rand <- numeric (500) 
for (ii in 1:500){
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  lscp <- as.matrix (lscp)
  samp <- sample (lscp, 10000, replace = F)
  tc <- as.data.frame (as.count (samp))
  sampRADs_gaus5_rand [ii, 1:nrow (tc)] <- tc$Freq
  SRs_gaus5_rand [ii] <- nrow (tc)
}
max (SRs_gaus5_rand)
sampRADs_gaus5_rand <- sampRADs_gaus5_rand [, 1:61]
sampRADs_gaus5_rand [is.na (sampRADs_gaus5_rand)] <- 0

Mean_rad_gaus5_rand <- numeric (ncol (sampRADs_gaus5_rand))
SE_hi_rad_gaus5_rand <- numeric (ncol (sampRADs_gaus5_rand))
SE_lo_rad_gaus5_rand <- numeric (ncol (sampRADs_gaus5_rand))
for (ii in 1:ncol (sampRADs_gaus5_rand)){
  Mean_rad_gaus5_rand [ii] <- mean (sampRADs_gaus5_rand [, ii])
  SE_hi_rad_gaus5_rand [ii] <- STDERR (sampRADs_gaus5_rand [, ii])[[1]]
  SE_lo_rad_gaus5_rand [ii] <- STDERR (sampRADs_gaus5_rand [, ii])[[3]]
}


## Dispersal limited, sigma = 2
sampRADs_gaus2_rand <- matrix (NA, 500, 200) 
SRs_gaus2_rand <- numeric (500) 
for (ii in 1:500){
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  lscp <- as.matrix (lscp)
  samp <- sample (lscp, 10000, replace = F)
  tc <- as.data.frame (as.count (samp))
  sampRADs_gaus2_rand [ii, 1:nrow (tc)] <- tc$Freq
  SRs_gaus2_rand [ii] <- nrow (tc)
}

max (SRs_gaus2_rand)
sampRADs_gaus2_rand <- sampRADs_gaus2_rand [, 1:61]
sampRADs_gaus2_rand [is.na (sampRADs_gaus2_rand)] <- 0

Mean_rad_gaus2_rand <- numeric (ncol (sampRADs_gaus2_rand))
SE_hi_rad_gaus2_rand <- numeric (ncol (sampRADs_gaus2_rand))
SE_lo_rad_gaus2_rand <- numeric (ncol (sampRADs_gaus2_rand))
for (ii in 1:ncol (sampRADs_gaus2_rand)){
  Mean_rad_gaus2_rand [ii] <- mean (sampRADs_gaus2_rand [, ii])
  SE_hi_rad_gaus2_rand [ii] <- STDERR (sampRADs_gaus2_rand [, ii])[[1]]
  SE_lo_rad_gaus2_rand [ii] <- STDERR (sampRADs_gaus2_rand [, ii])[[3]]
}

## Dispersal limited, sigma = 2, richness same as global
sampRADs_gaus2_tuneDown_rand <- matrix (NA, 500, 200) 
SRs_gaus2_tuneDown_rand <- numeric (500) 
for (ii in 1:500){
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  lscp <- as.matrix (lscp)
  samp <- sample (lscp, 10000, replace = F)
  tc <- as.data.frame (as.count (samp))
  sampRADs_gaus2_tuneDown_rand [ii, 1:nrow (tc)] <- tc$Freq
  SRs_gaus2_tuneDown_rand [ii] <- nrow (tc)
}

max (SRs_gaus2_tuneDown_rand)
sampRADs_gaus2_tuneDown_rand <- sampRADs_gaus2_tuneDown_rand [, 1:56] #was 53
sampRADs_gaus2_tuneDown_rand [is.na (sampRADs_gaus2_tuneDown_rand)] <- 0

Mean_rad_gaus_TD_rand <- numeric (ncol (sampRADs_gaus2_tuneDown_rand))
SE_hi_rad_gaus_TD_rand <- numeric (ncol (sampRADs_gaus2_tuneDown_rand))
SE_lo_rad_gaus_TD_rand <- numeric (ncol (sampRADs_gaus2_tuneDown_rand))
for (ii in 1:ncol (sampRADs_gaus2_tuneDown_rand)){
  Mean_rad_gaus_TD_rand [ii] <- mean (sampRADs_gaus2_tuneDown_rand [, ii])
  SE_hi_rad_gaus_TD_rand [ii] <- STDERR (sampRADs_gaus2_tuneDown_rand [, ii])[[1]]
  SE_lo_rad_gaus_TD_rand [ii] <- STDERR (sampRADs_gaus2_tuneDown_rand [, ii])[[3]]
}

### Aggregated sampling -- dispersed plots  
sampseq <- seq(1, 241, by = 60) # sampling "plots"

## Global dispersal
sampRADs_glob <- matrix (NA, 500, 200) # put the RADs for each sample
SRs_glob <- numeric (500) # SR of each sample
for (ii in 1:500){
  samp <- matrix (NA, 20, 1) # empty matrix that I will bind samples to
  lscp <- read.table (paste(file = "GL_compAMlscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  
  for (jj in sampseq){
    for  (kk in sampseq){
      tempsamp <- lscp [jj:(jj+19), kk:(kk+19)]
      samp <- cbind (samp, tempsamp) # these will end up being 20x500
    }
  }
  samp <- samp [, -1] # delete first NA col
  temprad <- as.data.frame(table(unlist(samp))) # make a table of the sample, i.e. the RAD
  SRs_glob [ii] <- nrow (temprad) # SR is the rows of the RAD table
  sampRADs_glob [ii, 1:nrow (temprad)] <- temprad [, 2] # Save RAD
  #  cat (ii, "\t")
}

# order
for (ii in 1:nrow (sampRADs_glob)){
  sampRADs_glob [ii, ] <- sort  (sampRADs_glob [ii, ], decreasing = T, na.last = T)
}

max (SRs_glob)
sampRADs_glob <- sampRADs_glob [, 1:55]
sampRADs_glob [is.na(sampRADs_glob)] <- 0

Mean_rad_glob <- numeric (ncol (sampRADs_glob))
SE_hi_rad_glob <- numeric (ncol (sampRADs_glob))
SE_lo_rad_glob <- numeric (ncol (sampRADs_glob))
for (ii in 1:ncol (sampRADs_glob)){
  Mean_rad_glob [ii] <- mean (sampRADs_glob [, ii])
  SE_hi_rad_glob [ii] <- STDERR (sampRADs_glob [, ii])[[1]]
  SE_lo_rad_glob [ii] <- STDERR (sampRADs_glob [, ii])[[3]]
}

## Gaussian, sigma = 5
sampRADs_gaus5 <- matrix (NA, 500, 200) 
SRs_gaus5 <- numeric (500) 
for (ii in 1:500){
  samp <- matrix (NA, 20, 1) 
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  
  for (jj in sampseq){
    for  (kk in sampseq){
      tempsamp <- lscp [jj:(jj+19), kk:(kk+19)]
      samp <- cbind (samp, tempsamp) 
    }
  }
  samp <- samp [, -1] 
  temprad <- as.data.frame(table(unlist(samp))) 
  SRs_gaus5 [ii] <- nrow (temprad) 
  sampRADs_gaus5 [ii, 1:nrow (temprad)] <- temprad [, 2] 
}

for (ii in 1:nrow (sampRADs_gaus5)){
  sampRADs_gaus5 [ii, ] <- sort  (sampRADs_gaus5 [ii, ], decreasing = T, na.last = T)
}

max (SRs_gaus5)
sampRADs_gaus5 <- sampRADs_gaus5 [, 1:55]
sampRADs_gaus5 [is.na(sampRADs_gaus5)] <- 0

Mean_rad_gaus5 <- numeric (ncol (sampRADs_gaus5))
SE_hi_rad_gaus5 <- numeric (ncol (sampRADs_gaus5))
SE_lo_rad_gaus5 <- numeric (ncol (sampRADs_gaus5))
for (ii in 1:ncol (sampRADs_gaus5)){
  Mean_rad_gaus5 [ii] <- mean (sampRADs_gaus5 [, ii])
  SE_hi_rad_gaus5 [ii] <- STDERR (sampRADs_gaus5 [, ii])[[1]]
  SE_lo_rad_gaus5 [ii] <- STDERR (sampRADs_gaus5 [, ii])[[3]]
}


## Gaussian, sigma = 2
sampRADs_gaus2 <- matrix (NA, 500, 200) 
SRs_gaus2 <- numeric (500) 
for (ii in 1:500){
  samp <- matrix (NA, 20, 1) 
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  
  for (jj in sampseq){
    for  (kk in sampseq){
      tempsamp <- lscp [jj:(jj+19), kk:(kk+19)]
      samp <- cbind (samp, tempsamp) 
    }
  }
  samp <- samp [, -1] 
  temprad <- as.data.frame(table(unlist(samp))) 
  SRs_gaus2 [ii] <- nrow (temprad) 
  sampRADs_gaus2 [ii, 1:nrow (temprad)] <- temprad [, 2] 
  # cat (ii, "\t")
}

for (ii in 1:nrow (sampRADs_gaus2)){
  sampRADs_gaus2 [ii, ] <- sort  (sampRADs_gaus2 [ii, ], decreasing = T, na.last = T)
}

max (SRs_gaus2)
sampRADs_gaus2 <- sampRADs_gaus2 [, 1:51]
sampRADs_gaus2 [is.na(sampRADs_gaus2)] <- 0

Mean_rad_gaus2 <- numeric (ncol (sampRADs_gaus2))
SE_hi_rad_gaus2 <- numeric (ncol (sampRADs_gaus2))
SE_lo_rad_gaus2 <- numeric (ncol (sampRADs_gaus2))
for (ii in 1:ncol (sampRADs_gaus2)){
  Mean_rad_gaus2 [ii] <- mean (sampRADs_gaus2 [, ii])
  SE_hi_rad_gaus2 [ii] <- STDERR (sampRADs_gaus2 [, ii])[[1]]
  SE_lo_rad_gaus2 [ii] <- STDERR (sampRADs_gaus2 [, ii])[[3]]
}


## Gausian, sigma = 2, global SR
sampRADs_gaus_TU <- matrix (NA, 500, 200) 
SRs_gaus_TU <- numeric (500) 
for (ii in 1:500){
  samp <- matrix (NA, 20, 1) 
  lscp <- read.table (paste(file = "Gaus_jm105lscp_1.0e+06_1_", ii, ".dat", sep = ""), 
                      fill = T, header = F)
  
  for (jj in sampseq){
    for  (kk in sampseq){
      tempsamp <- lscp [jj:(jj+19), kk:(kk+19)]
      samp <- cbind (samp, tempsamp) 
    }
  }
  samp <- samp [, -1] 
  temprad <- as.data.frame(table(unlist(samp))) 
  SRs_gaus_TU [ii] <- nrow (temprad) 
  sampRADs_gaus_TU [ii, 1:nrow (temprad)] <- temprad [, 2] 
  #  cat (ii, "\t")
}

for (ii in 1:nrow (sampRADs_gaus_TU)){
  sampRADs_gaus_TU [ii, ] <- sort  (sampRADs_gaus_TU [ii, ], decreasing = T, na.last = T)
}

max (SRs_gaus_TU)
sampRADs_gaus_TU <- sampRADs_gaus_TU [, 1:53]
sampRADs_gaus_TU [is.na (sampRADs_gaus_TU)] <- 0


Mean_rad_gaus_TU <- numeric (ncol (sampRADs_gaus_TU))
SE_hi_rad_gaus_TU <- numeric (ncol (sampRADs_gaus_TU))
SE_lo_rad_gaus_TU <- numeric (ncol (sampRADs_gaus_TU))
for (ii in 1:ncol (sampRADs_gaus_TU)){
  Mean_rad_gaus_TU [ii] <- mean (sampRADs_gaus_TU [, ii])
  SE_hi_rad_gaus_TU [ii] <- STDERR (sampRADs_gaus_TU [, ii])[[1]]
  SE_lo_rad_gaus_TU [ii] <- STDERR (sampRADs_gaus_TU [, ii])[[3]]
}


### Hyperdominance

## Global aggregated
sims <- nrow (sampRADs_glob)
glob_samp_cumu50 <- numeric (sims)
glob_samp_cumu50_Prop <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_glob [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  glob_samp_cumu50 [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  glob_samp_cumu50_Prop [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1)/length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_glob_samp <- STDERR(glob_samp_cumu50)
MSE_glob_samp_Prop <- STDERR (glob_samp_cumu50_Prop)
CIs_glob_samp_prop <- quantile(glob_samp_cumu50_Prop, c(0.025, 0.975))

## Global random 
sims <- nrow (sampRADs_glob_rand)
glob_samp_cumu50_rand <- numeric (sims)
glob_samp_cumu50_Prop_rand <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_glob_rand [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  glob_samp_cumu50_rand [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  glob_samp_cumu50_Prop_rand [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1)/length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_glob_samp_rand <- STDERR(glob_samp_cumu50_rand)
MSE_glob_samp_Prop_rand <- STDERR (glob_samp_cumu50_Prop_rand)
CIs_glob_samp_prop_rand <- quantile(glob_samp_cumu50_Prop_rand, c(0.025, 0.975))

## sigma = 5 aggregated
gaus_samp_sig5_cumu50 <- numeric (sims)
gaus_samp_sig5_cumu50_Prop <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_gaus5 [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig5_cumu50 [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig5_cumu50_Prop [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1) / length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_gaus_samp_sig5 <- STDERR(gaus_samp_sig5_cumu50)
MSE_gaus_samp_sig5_Prop <- STDERR (gaus_samp_sig5_cumu50_Prop)
CIs_gaus_samp_sig5_Prop <- quantile(gaus_samp_sig5_cumu50_Prop, c(0.025, 0.975))


## sigma = 5 random 
gaus_samp_sig5_cumu50_rand <- numeric (sims)
gaus_samp_sig5_cumu50_Prop_rand <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_gaus5_rand [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig5_cumu50_rand [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig5_cumu50_Prop_rand [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1) / length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_gaus_samp_sig5_rand <- STDERR(gaus_samp_sig5_cumu50_rand)
MSE_gaus_samp_sig5_Prop_rand <- STDERR (gaus_samp_sig5_cumu50_Prop_rand)
CIs_gaus_samp_sig5_Prop_rand <- quantile(gaus_samp_sig5_cumu50_Prop_rand, c(0.025, 0.975))


## sigma = 2, aggregated
gaus_samp_sig2_cumu50 <- numeric (sims)
gaus_samp_sig2_cumu50_Prop <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_gaus2 [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig2_cumu50 [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig2_cumu50_Prop [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1) / length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}

MSE_gaus_samp_sig2 <- STDERR (gaus_samp_sig2_cumu50)
MSE_gaus_samp_sig2_Prop <- STDERR (gaus_samp_sig2_cumu50_Prop)
CIs_gaus_samp_sig2_Prop <- quantile(gaus_samp_sig2_cumu50_Prop, c(0.025, 0.975))

## sigma = 2, random
gaus_samp_sig2_cumu50_rand <- numeric (sims)
gaus_samp_sig2_cumu50_Prop_rand <- numeric (sims)
for (ii in 1:sims){
  temp_sim <- as.numeric (sampRADs_gaus2_rand [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig2_cumu50_rand [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig2_cumu50_Prop_rand [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1) / length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_gaus_samp_sig2_rand <- STDERR (gaus_samp_sig2_cumu50_rand)
MSE_gaus_samp_sig2_Prop_rand <- STDERR (gaus_samp_sig2_cumu50_Prop_rand)
CIs_gaus_samp_sig2_Prop_rand <- quantile(gaus_samp_sig2_cumu50_Prop_rand, c(0.025, 0.975))


## sigma = 2, global SR, aggregated
gaus_samp_sig2_cumu50_T <- numeric (sims)
gaus_samp_sig2_cumu50_T_Prop <- numeric (sims)
for (ii in 1:500){
  temp_sim <- as.numeric (sampRADs_gaus_TU [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig2_cumu50_T [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig2_cumu50_T_Prop [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1)/length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_gaus_samp_sig2_T <- STDERR (gaus_samp_sig2_cumu50_T)
MSE_gaus_samp_sig2_T_Prop <- STDERR (gaus_samp_sig2_cumu50_T_Prop)
CIs_gaus_samp_sig2_T_Prop <- quantile(gaus_samp_sig2_cumu50_T_Prop, c(0.025, 0.975))


## sigma = 2, global SR, random 
gaus_samp_sig2_cumu50_tuneDown_rand <- numeric (sims)
gaus_samp_sig2_cumu50_tuneDown_Prop_rand <- numeric (sims)
for (ii in 1:500){
  temp_sim <- as.numeric (sampRADs_gaus2_tuneDown_rand [ii, ])
  temp_sim <- temp_sim[!is.na(temp_sim)]
  temp_sim <- temp_sim [temp_sim > 0]
  gaus_samp_sig2_cumu50_tuneDown_rand [ii] <- (length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1
  gaus_samp_sig2_cumu50_tuneDown_Prop_rand [ii] <- ((length(cumsum(temp_sim)[cumsum(temp_sim) < sum(temp_sim)/2])) + 1)/length (temp_sim)
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}
MSE_gaus_samp_sig2_tuneDown_rand <- STDERR (gaus_samp_sig2_cumu50_tuneDown_rand)
MSE_gaus_samp_sig2_tuneDown_rand_prop <- STDERR (gaus_samp_sig2_cumu50_tuneDown_Prop_rand)
CIs_gaus_samp_sig2_tuneDown_rand_prop <- quantile(gaus_samp_sig2_cumu50_tuneDown_Prop_rand, c(0.025, 0.975))


## Figure 4
windows (6.8, 6)
par (mfrow = c(2, 2))
par (oma = c(1, 1, 1, 1))
par(mar = c(3,3,1,1))
# Random sample RAD
plot (Mean_rad_gaus5_rand, type = "l", lwd = 1, col = "purple", xlim = c(1, 62), log = "y", ylab = " ", cex.axis = 1.2, cex.lab = 1.2,
      xlab = "Rank", main = "Random Sample")#,  ylab = "Abundance"
mtext("Abundance", side = 2, line = 3)
polygon(x = c(1:60, rev(1:60)),
        y = c(SE_hi_rad_gaus5_rand [1:60], rev(SE_lo_rad_gaus5_rand [1:60])),
        col =  adjustcolor("purple", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus2_rand, col = "forestgreen", type = "l", lwd = 1)
polygon(x = c(1:60, rev(1:60)),
        y = c(SE_hi_rad_gaus2_rand [1:60], rev(SE_lo_rad_gaus2_rand [1:60])),
        col =  adjustcolor("forestgreen", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus_TD_rand, col = "orange", type = "l", lwd = 1)
polygon(x = c(1:53, rev(1:53)),
        y = c(SE_hi_rad_gaus_TD_rand [1:53], rev(SE_lo_rad_gaus_TD_rand [1:53])),
        col =  adjustcolor("orange", alpha.f = 0.30), border = NA)
points (Mean_rad_glob_rand, type = "l", lwd = 1)
polygon(x = c(1:53, rev(1:53)),
        y = c(SE_hi_rad_glob_rand [1:53], rev(SE_lo_rad_glob_rand [1:53])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
mtext("Rank", side = 1, line = 2.3)

# Aggregated sample RAD
plot (Mean_rad_gaus5, log = "y", type = "l", xlab = "Rank", col = "purple",  ylab = " ", cex.axis = 1.2, cex.lab = 1.2,
      lwd = 1, main = "Aggregated Sample", xlim = c(1, 62)) #ylim = c(.005, max (SE_hi_rad_gaus5))
polygon(x = c(1:53, rev(1:53)),
        y = c(SE_hi_rad_gaus5 [1:53], rev(SE_lo_rad_gaus5 [1:53])),
        col =  adjustcolor("purple", alpha.f = 0.30), border = NA)
points (Mean_rad_glob, type = "l", lwd = 1)
polygon(x = c(1:54, rev(1:54)),
        y = c(SE_hi_rad_glob [1:54], rev(SE_lo_rad_glob [1:54])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus_TU, type = "l", col = "orange", lwd = 1)
polygon(x = c(1:53, rev(1:53)),
        y = c(SE_hi_rad_gaus_TU [1:53], rev(SE_lo_rad_gaus_TU [1:53])),
        col =  adjustcolor("orange", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus2, type = "l", col = "forestgreen", lwd = 1)
polygon(x = c(1:48, rev(1:48)),
        y = c(SE_hi_rad_gaus2 [1:48], rev(SE_lo_rad_gaus2 [1:48])),
        col =  adjustcolor("forestgreen", alpha.f = 0.30), border = NA)
mtext("Rank", side = 1, line = 2.3)

# Random sample hyperdominants
plot (x = 1, y = MSE_glob_samp_Prop_rand [2], pch = 19 , xaxt = "n",  xlim = c(0.9, 1.4), ylab = "", cex.axis = 1.2, cex.lab = 1.2,
      ylim = c(.07, .11), xlab  = "", main = " ")
axis (1, at = c(1, 1.1, 1.2, 1.3), labels = c("", "", "", ""), las = 2)
points (x = 1.1, y = MSE_gaus_samp_sig5_Prop_rand [2], pch = 19, col = "purple")
points (x = 1.2, y = MSE_gaus_samp_sig2_Prop_rand [2], pch = 19, col = "forestgreen")
points (x = 1.3, y = MSE_gaus_samp_sig2_tuneDown_rand_prop [2], pch = 19, col = "orange")

arrows(1, MSE_glob_samp_Prop_rand [1], 1, MSE_glob_samp_Prop_rand [3], length = 0, lwd = 2)
arrows(1.1, MSE_gaus_samp_sig5_Prop_rand[1], 1.1, MSE_gaus_samp_sig5_Prop_rand [3], length = 0, col = "purple", lwd = 2)
arrows(1.2, MSE_gaus_samp_sig2_Prop_rand [1], 1.2, MSE_gaus_samp_sig2_Prop_rand [3], length = 0, col = "forestgreen", lwd = 2)
arrows (1.3, MSE_gaus_samp_sig2_tuneDown_rand_prop [1], 1.3, MSE_gaus_samp_sig2_tuneDown_rand_prop [3], length = 0, col = "orange", lwd = 2)
mtext("Proportion Hyperdominant", side = 2, line = 3)

# Aggregated sample hyperdominants
plot (x = 1, y = MSE_glob_samp_Prop [2], pch = 19 , xaxt = "n", xlim = c(0.9, 1.4), ylim = c(.07, .11), xlab  = "", 
      ylab = " ", cex.axis = 1.2, cex.lab = 1.2,) 
axis (1, at = c(1, 1.1, 1.2, 1.3), labels = c("", "", "", ""), las = 2)
points (x = 1.1, y = MSE_gaus_samp_sig5_Prop [2], pch = 19, col = "purple")
points (x = 1.2, y = MSE_gaus_samp_sig2_Prop [2], pch = 19, col = "forestgreen")
points (x = 1.3, y = MSE_gaus_samp_sig2_T_Prop [2], pch = 19, col = "orange")

arrows(1, MSE_glob_samp_Prop [1], 1, MSE_glob_samp_Prop [3], length = 0, lwd = 2)
arrows(1.1, MSE_gaus_samp_sig5_Prop[1], 1.1, MSE_gaus_samp_sig5_Prop [3], length = 0, col = "purple", lwd = 2)
arrows(1.2, MSE_gaus_samp_sig2_Prop [1], 1.2, MSE_gaus_samp_sig2_Prop [3], length = 0, col = "forestgreen", lwd = 2)
arrows (1.3, MSE_gaus_samp_sig2_T_Prop [1], 1.3, MSE_gaus_samp_sig2_T_Prop [3], length = 0, col = "orange", lwd = 2)

# Random sample inset
par (cex = 0.7,lwd = 0.5, fig = c(0.067, 0.3, 0.573, 0.8), new = T)
plot (Mean_rad_gaus5_rand, type = "l", lwd = 1, col = "purple", xlim = c(1, 4), xlab = " ", ylab = " ",
      ylim = c(800, 3200), cex.axis = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus5_rand [1:4], rev(SE_lo_rad_gaus5_rand [1:4])),
        col =  adjustcolor("purple", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus2_rand, col = "forestgreen", type = "l", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus2_rand [1:4], rev(SE_lo_rad_gaus2_rand [1:4])),
        col =  adjustcolor("forestgreen", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus_TD_rand, col = "orange", type = "l", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus_TD_rand [1:4], rev(SE_lo_rad_gaus_TD_rand [1:4])),
        col =  adjustcolor("orange", alpha.f = 0.30), border = NA)
points (Mean_rad_glob_rand, type = "l", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_glob_rand [1:4], rev(SE_lo_rad_glob_rand [1:4])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)

# Aggregated sample inset
par (cex = 0.7,lwd = 0.5, fig = c(0.567, 0.79, 0.573, 0.8), new = T) 
plot (Mean_rad_gaus5, type = "l", xlab = " ", ylab = " ", col = "purple",
      lwd = 1, xlim = c(1, 4), ylim = c(800, 3200), cex.axis = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus5 [1:4], rev(SE_lo_rad_gaus5 [1:4])),
        col =  adjustcolor("purple", alpha.f = 0.30), border = NA)
points (Mean_rad_glob, type = "l", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_glob [1:4], rev(SE_lo_rad_glob [1:4])),
        col =  adjustcolor("black", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus_TU, type = "l", col = "orange", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus_TU [1:4], rev(SE_lo_rad_gaus_TU [1:4])),
        col =  adjustcolor("orange", alpha.f = 0.30), border = NA)
points (Mean_rad_gaus2, type = "l", col = "forestgreen", lwd = 1)
polygon(x = c(1:4, rev(1:4)),
        y = c(SE_hi_rad_gaus2 [1:4], rev(SE_lo_rad_gaus2 [1:4])),
        col =  adjustcolor("forestgreen", alpha.f = 0.30), border = NA)
#

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', c("Global dispersal","Disperal dist = 25m", "Dispersal dist = 10m", "Dispersal dist = 10m, global SR"), 
       col = c("black", "purple", "forestgreen", "orange"),lwd = 3, xpd = TRUE, horiz = F, cex = 1.2,
       ncol = 2, seg.len=1, bty = 'n')



########## Speciation mode ##########
# modified SADISA_sim function at the bottom of script

ATDN <- read.csv("terSteege_2019_bias_corrected_data.csv", stringsAsFactors = F)
abund <- ATDN$n.ind

### Protracted
# More protracted casee: phi = 7000, theta = 8000, dispersal = 10e6 
PR_phi7000_theta8180 <- SADISA_sim2(parsmc = c(8180, 7000), ii = 10e6, jj = 979614, model = c("pr", "dl"), mult = "single", nsim = 100) 

phi7000theta8180mat <- matrix(data = NA, nrow = 100, ncol = 5210)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (PR_phi7000_theta8180 [[ii]])), decreasing = T)
  length <- length (temp)
  phi7000theta8180mat [ii, 1:length] <- temp
}

write.csv(phi7000theta8180mat)
PRO_phi7000 <- read.csv (file = "phi7000theta8180mat.csv", stringsAsFactors = F)

# RAD and hyperdominance
SRproM <- c()
for (ii in 1:nrow (PRO_phi7000)){
  temp <- PRO_phi7000 [ii, ]
  temp <- temp [!is.na (temp)]
  SRproM [ii] <- length (temp)
}

hypePM <- c()
prop_hypPM <- c()
for (ii in 1:nrow (PRO_phi7000)){
  temp <- PRO_phi7000 [1, ]
  temp <- temp [!is.na (temp)]
  hypePM [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypPM [ii] <- hypePM [ii]/SRproM [ii]
}
ProM_hyp_mean <- mean(prop_hypPM)
ProM_hyp_SE_low <- quantile(prop_hypPM, 0.025)
ProM_hyp_SE_high <- quantile(prop_hypPM, 0.925)


## RAD mean and SE
ProM_mean <- numeric (ncol (PRO_phi7000))
ProM_CI_025 <- numeric (ncol (PRO_phi7000))
ProM_CI_975 <- numeric (ncol (PRO_phi7000))
for (ii in 1:ncol (PRO_phi7000)){
  temp_vec <- PRO_phi7000 [, ii]
  ProM_mean [ii] <- mean (temp_vec, na.rm = T)
  ProM_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  ProM_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}


# Less protracted case: phi = 10e5, theta = 782, dispersal = 10e6  
phi_1e6_t782 <- SADISA_sim2(parsmc = c(782, 1000000), ii = 10e6, jj = 979614, model = c("pr", "dl"), mult = "single", nsim = 100)
Pro_lessnew_mat <- matrix(data = NA, nrow = 100, ncol = 5249)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (phi_1e6_t782 [[ii]])), decreasing = T)
  length <- length (temp)
  Pro_lessnew_mat [ii, 1:length] <- temp
}

write.csv(Pro_lessnew_mat)
PRO_less <- read.csv ("Pro_less_new_t782_mat.csv", stringsAsFactors = F)

# RAD and hyperdominance
SRproL <- c()
for (ii in 1:nrow (PRO_phi1e6)){
  temp <- PRO_phi1e6 [ii, ]
  temp <- temp [!is.na (temp)]
  SRproL [ii] <- length (temp)
}

hypePL <- c()
prop_hypPL <- c()
for (ii in 1:nrow (PRO_less)){
  temp <- PRO_phi1e6 [1, ]
  temp <- temp [!is.na (temp)]
  hypePL [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypPL [ii] <- hypePL [ii]/SRproL [ii]
}
PL_hyp_mean <- mean(prop_hypPL)
PL_hyp_SE_low <- quantile(prop_hypPL, 0.025)
PL_hyp_SE_high <- quantile(prop_hypPL, 0.925)

## RAD mean and SE
ProL_mean <- numeric (ncol (PRO_phi1e6))
ProL_CI_025 <- numeric (ncol (PRO_phi1e6))
ProL_CI_975 <- numeric (ncol (PRO_phi1e6))
for (ii in 1:ncol (PRO_phi1e6)){
  temp_vec <- PRO_phi1e6 [, ii]
  ProL_mean [ii] <- mean (temp_vec, na.rm = T)
  ProL_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  ProL_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}


### Abundance-independent ("DD")
# Less abundance-independent: alpha = 0.2, theta = 237, dispersal = 10e6 
DD_theta237_a0.2 <- SADISA_sim2(parsmc = c(237, 0.2), ii = 10e6, jj = 979614, model = c("dd", "dl"), mult = "single", nsim = 100)

DD_theta237_a0.2mat <- matrix(data = NA, nrow = 100, ncol = 5181)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (DD_theta237_a0.2 [[ii]])), decreasing = T)
  length <- length (temp)
  DD_theta237_a0.2mat [ii, 1:length] <- temp
}

write.csv(DD_theta237_a0.2mat)
DD_theta237_a0.2 <- read.csv ("DD_theta6237_a0.2.csv", header = T, stringsAsFactors = F)

## hyperdominance and RAD
SRDD_l <- c()
for (ii in 1:nrow (DD_theta237_a0.2)){
  temp <- DD_theta237_a0.2 [ii, ]
  temp <- temp [!is.na (temp)]
  SRDD_l [ii] <- length (temp)
}

hypeDD_l <- c()
prop_hypDD_l <- c()
for (ii in 1:nrow (DD_theta237_a0.2)){
  temp <- DD_theta237_a0.2 [ii, ]
  temp <- temp [!is.na (temp)]
  hypeDD_l [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypDD_l [ii] <- hypeDD_l [ii]/SRDD_l [ii]
}
DD_l_meanprop <- mean (prop_hypDD_l) 
DD_l_prop_ciHigh <- quantile (prop_hypDD_l, 0.975)
DD_l_prop_cilow <- mean (prop_hypDD_l, 0.025)

DD_l_mean <- numeric (ncol (DD_theta237_a0.2))
DD_l_CI_025 <- numeric (ncol (DD_theta237_a0.2))
DD_l_CI_975 <- numeric (ncol (DD_theta237_a0.2))
for (ii in 1:ncol (DD_theta237_a0.2)){
  temp_vec <- DD_theta237_a0.2 [, ii]
  DD_l_mean [ii] <- mean (temp_vec, na.rm = T)
  DD_l_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  DD_l_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}


# More abundance-independent: alpha = 0.5, theta = 7, dispersal = 10e6 
DD_theta7_a0.5 <- SADISA_sim2(parsmc = c(7, 0.5), ii = 10e6, jj = 979614, model = c("dd", "dl"), mult = "single", nsim = 100)
DD_theta7_a0.5mat <- matrix(data = NA, nrow = 100, ncol = 5284)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (DD_theta7_a0.5 [[ii]])), decreasing = T)
  length <- length (temp)
  DD_theta7_a0.5mat [ii, 1:length] <- temp
}

write.csv(DD_theta7_a0.5mat)
DD_theta7_a0.5 <- read.csv ("DD_theta7_a0.5.csv", header = T, stringsAsFactors = F)
#RAD and hyperdominance
SRDD_m <- c()
for (ii in 1:nrow (DD_theta7_a0.5)){
  temp <- DD_theta7_a0.5 [ii, ]
  temp <- temp [!is.na (temp)]
  SRDD_m [ii] <- length (temp)
}

hypeDD_m <- c()
prop_hypDD_m <- c()
for (ii in 1:nrow (DD_theta7_a0.5)){
  temp <- DD_theta7_a0.5 [1, ]
  temp <- temp [!is.na (temp)]
  hypeDD_m [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypDD_m [ii] <- hypeDD_m [ii]/SRDD_m [ii]
}
DD_m_hyp_mean <- mean(prop_hypDD_m)
DD_m_hyp_SE_low <- quantile(prop_hypDD_m, 0.025)
DD_m_hyp_SE_high <- quantile(prop_hypDD_m, 0.925)

DD_m_mean <- numeric (ncol (DD_theta7_a0.5))
DD_m_CI_025 <- numeric (ncol (DD_theta7_a0.5))
DD_m_CI_975 <- numeric (ncol (DD_theta7_a0.5))
for (ii in 1:ncol (DD_theta7_a0.5)){
  temp_vec <- DD_theta7_a0.5 [, ii]
  DD_m_mean [ii] <- mean (temp_vec, na.rm = T)
  DD_m_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  DD_m_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}


### Point mutation
PM_baseline <- SADISA_sim2(parsmc = c(692.8512), ii = 10e6, jj = 979614, model = c("pm", "dl"), mult = "single", nsim = 100)

PM_baselinemat <- matrix(data = NA, nrow = 100, ncol = 5188)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (PM_baseline [[ii]])), decreasing = T)
  length <- length (temp)
  PM_baselinemat [ii, 1:length] <- temp
}

write.csv(PM_baselinemat)
PM_baseline <- read.csv ("PM_baselineMat.csv", header = T, stringsAsFactors = F)

# RAD and hyperdominance
SR_point <- c()
for (ii in 1:nrow (PM_baseline)){
  temp <- PM_baseline [ii, ]
  temp <- temp [!is.na (temp)]
  SR_point [ii] <- length (temp)
}

hype_point <- c()
prop_hyp_point <- c()
for (ii in 1:nrow (PM_baseline)){
  temp <- PM_baseline [ii, ]
  temp <- temp [!is.na (temp)]
  hype_point [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hyp_point [ii] <- hype_point [ii]/SR_point [ii]
}
mean_propHype_PM <- mean (prop_hyp_point)
CIhigh_prophype_PM <- quantile(prop_hyp_point, 0.975)
CIlow_prophype_PM <- quantile(prop_hyp_point, 0.025)

Point_mean <- numeric (ncol (PM_baseline))
Point_CI_025 <- numeric (ncol (PM_baseline))
Point_CI_975 <- numeric (ncol (PM_baseline))
for (ii in 1:ncol (PM_baseline)){
  temp_vec <- PM_baseline [, ii]
  Point_mean [ii] <- mean (temp_vec, na.rm = T)
  Point_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  Point_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

### ML best fit parameters 

# protracted
model <- c("pr", "dl")
init_params <- c(8000, 7000, 10e6) #theta, phi, I
id_pars <- c(1, 1, 0) # optimize theta and alpha but not I
label_pars <- c(1, 2, 3) 
Pr_test <- SADISA_ML(abund = abund, initpars = init_params, idpars = id_pars, labelpars = label_pars, model = model)

#Parameters after likelihood maximization:
#7.425198e+02 2.102163e+06 1.000000e+07
#Maximum loglikelihood:
#  -2305.203

# simulate
probf <- SADISA_sim2(parsmc = c(743, 2100000), ii = 10e6, jj = 979614, model = c("pr", "dl"), mult = "single", nsim = 100)

Pro_BF_mat <- matrix(data = NA, nrow = 100, ncol = 5246)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (probf [[ii]])), decreasing = T)
  length <- length (temp)
  Pro_BF_mat [ii, 1:length] <- temp
}

write.csv(Pro_BF_mat)
Pro_BF_mat <- read.csv (file = "Pro_BF_mat.csv", stringsAsFactors = F)

#RAD and hyperdominance
SRproBF <- c()
for (ii in 1:nrow (Pro_BF_mat)){
  temp <- Pro_BF_mat [ii, ]
  temp <- temp [!is.na (temp)]
  SRproBF [ii] <- length (temp)
}
hypeProBF <- c()
prop_hypProBF <- c()
for (ii in 1:nrow (Pro_BF_mat)){
  temp <- Pro_BF_mat [1, ]
  temp <- temp [!is.na (temp)]
  hypeProBF [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypProBF [ii] <- hypeProBF [ii]/SRproBF [ii]
}
ProBF_hyp_mean <- mean(prop_hypProBF)
ProBF_hyp_SE_low <- quantile(prop_hypProBF, 0.025)
ProBF_hyp_SE_high <- quantile(prop_hypProBF, 0.925)

ProBF_mean <- numeric (ncol (Pro_BF_mat))
ProBF_CI_025 <- numeric (ncol (Pro_BF_mat))
ProBF_CI_975 <- numeric (ncol (Pro_BF_mat))
for (ii in 1:ncol (Pro_BF_mat)){
  temp_vec <- Pro_BF_mat [, ii]
  ProBF_mean [ii] <- mean (temp_vec, na.rm = T)
  ProBF_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  ProBF_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

# Abundance independent best fit 
model <- c("dd", "dl")
init_params <- c(100, 0.2, 10e6) #theta, alpha, I
id_pars <- c(1, 1, 0) 
label_pars <- c(1, 2, 3) 
ML_test2 <- SADISA_ML(abund = abund, initpars = init_params, idpars = id_pars, labelpars = label_pars, model = model)

## Parameters after likelihood maximization:
#[1]  8.147916e+02 -3.856195e-02  1.000000e+07

#### NOTE: best fit alpha is negative -- this is NOT abundance-independent speciation
# It probably means negative ecological density-dependence given the relationship between dd speciation and ecological dd 

BestFitParamsSim <- SADISA_sim2(parsmc = c(814.7916, -0.03856195), ii = 10e6, jj = 979614, model = c("dd", "dl"), mult = "single", nsim = 100)

BF_mat <- matrix(data = NA, nrow = 100, ncol = 5192)
for (ii in 1:100){
  temp <- sort(as.numeric(unlist (BestFitParamsSim [[ii]])), decreasing = T)
  length <- length (temp)
  BF_mat [ii, 1:length] <- temp
}
write.csv(BF_mat)
DD_BF_Mat <- read.csv ("BestFit_mat.csv", header = T, stringsAsFactors = F)

# RAD and hyperdominance
SRDD_bf <- c()
for (ii in 1:nrow (DD_BF_Mat)){
  temp <- DD_BF_Mat [ii, ]
  temp <- temp [!is.na (temp)]
  SRDD_bf [ii] <- length (temp)
}

hypeDD_bf <- c()
prop_hypDD_bf <- c()
for (ii in 1:nrow (DD_BF_Mat)){
  temp <- DD_BF_Mat [ii, ]
  temp <- temp [!is.na (temp)]
  hypeDD_bf [ii] <- (length (cumsum(temp)[cumsum(temp) < sum (temp)/2])) + 1
  prop_hypDD_bf [ii] <- hypeDD_bf [ii]/SRDD_bf [ii]
}
DD_bf_meanprop <- mean (prop_hypDD_bf) 
DD_bf_prop_ciHigh <- quantile (prop_hypDD_bf, 0.975)
DD_bf_prop_cilow <- mean (prop_hypDD_bf, 0.025)

DD_bf_mean <- numeric (ncol (DD_BF_Mat))
DD_bf_CI_025 <- numeric (ncol (DD_BF_Mat))
DD_bf_CI_975 <- numeric (ncol (DD_BF_Mat))
for (ii in 1:ncol (DD_BF_Mat)){
  temp_vec <- DD_BF_Mat [, ii]
  DD_bf_mean [ii] <- mean (temp_vec, na.rm = T)
  DD_bf_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  DD_bf_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

### Fig. 5
windows (6.8, 6)
vir <- viridis(7)
par (mfrow = c(2, 2))
# two scenario RAD
par (mar = c(4, 4, 1, 1))
plot (DD_m_mean, type = "l", col = vir[5], log = "xy", xlab = "Rank", ylab = "Abundance",
      lwd = 2, ylim = c(1, max(DD_m_CI_975)))
polygon(x = c(1:5000, rev(1:5000)),
        y = c(DD_m_CI_975 [1:5000], rev(DD_m_CI_025 [1:5000])),
        col =  adjustcolor(vir[5], alpha.f = 0.30), border = NA)
#DD less
points (DD_l_mean, type = "l", col = vir[6], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(DD_l_CI_975 [1:5000], rev(DD_l_CI_025 [1:5000])),
        col =  adjustcolor(vir[6], alpha.f = 0.30), border = NA)
#Pro more
points (ProM_mean, type = "l", col = vir[2], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(ProM_CI_975 [1:5000], rev(ProM_CI_025 [1:5000])),
        col =  adjustcolor(vir[2], alpha.f = 0.30), border = NA)
#Pro less
points (ProL_mean, type = "l", col = vir[3], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(ProL_CI_975 [1:5000], rev(ProL_CI_025 [1:5000])),
        col =  adjustcolor(vir[3], alpha.f = 0.30), border = NA)
# point
points (Point_mean, type = "l", col = vir[7], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(Point_CI_975 [1:5000], rev(Point_CI_025 [1:5000])),
        col =  adjustcolor(vir[7], alpha.f = 0.30), border = NA)
#empirical
points (abund, type = "l", lwd = 2, col = "black")# vir[1])

#best fit RAD
vir <- viridis (n = 4)
par (mar = c(4, 4, 1, 2))
vir <- viridis (n = 4)
plot (abund, type = "l", lwd = 2, col = "black", log = "xy", xlab = "Rank", ylab = " ",
      ylim = c(1, max(DD_m_CI_975)))
#DD beft fit
points (DD_bf_mean, type = "l", col = vir[3], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(DD_bf_CI_975 [1:5000], rev(DD_bf_CI_025 [1:5000])),
        col =  adjustcolor(vir[2], alpha.f = 0.30), border = NA)
#Pro best fit
points (ProBF_mean, type = "l", col = vir[2], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(ProBF_CI_975 [1:5000], rev(ProBF_CI_025 [1:5000])),
        col =  adjustcolor(vir[3], alpha.f = 0.30), border = NA)
# point
points (Point_mean, type = "l", col = vir[4], lwd = 2)
polygon(x = c(1:5000, rev(1:5000)),
        y = c(Point_CI_975 [1:5000], rev(Point_CI_025 [1:5000])),
        col =  adjustcolor(vir[4], alpha.f = 0.30), border = NA)

#two scenario proportion hyperdominance
emp_prophyp <- (length((cumsum(abund)[cumsum(abund) < sum (abund)/2]))+1)/5027
vir <- viridis(7)
par (mar = c(6, 4, 1, 1))
plot (x = 1, y = mean_propHype_PM, pch = 19 , xaxt = "n", ylab = "Proportion hyperdominants", xlim = c(0.9, 1.5), ylim = c(0, .20), 
      xlab  = "", col = vir[7])
axis (1, at = c(1, 1.1, 1.2, 1.3, 1.4), labels = c("Point", "AI_less", "AI_more", "Pro_less", "Pro_more"), las = 2)
abline(emp_prophyp, 0, col = "black", lwd = 2)
points (x = 1.1, y =DD_l_meanprop, pch = 19, col = vir[6])
points (x = 1.2, y =DD_m_hyp_mean, pch = 19, col = vir[5])
points (x = 1.3, y = PL_hyp_mean, pch = 19, col = vir[3])
points (x = 1.4, y = ProM_hyp_mean, pch = 19, col = vir[2])

arrows(1, CIhigh_prophype_PM, 1, CIlow_prophype_PM, length = 0, lwd = 2, col = vir[7])
arrows(1.1, DD_l_prop_ciHigh, 1.1, DD_l_prop_cilow, length = 0, col = vir[6], lwd = 2)
arrows(1.2, DD_m_hyp_SE_high, 1.2, DD_m_hyp_SE_low, length = 0, col = vir[5], lwd = 2)
arrows (1.3, PL_hyp_mean, 1.3, PL_hyp_SE_low, length = 0, col = vir[3], lwd = 2)
arrows (1.4, ProM_hyp_SE_high, 1.4, ProM_hyp_SE_low, length = 0, col = vir[2], lwd = 2)
legend ("topleft", c("Empirical", "Point", "AI_less", "AI_more", "Pro_less", "Pro_more"),
        col = c("black", vir[7], vir[6], vir[5], vir[3], vir[2]), lwd = c(2, 2, 2, 2, 2, 2), bty = "n")

#best fit proportion hyperdominance
par (mar = c(6, 4, 1, 1))
vir <- viridis (n = 4)
plot (x = 1, y = mean_propHype_PM, pch = 19 , xaxt = "n", ylab = "", xlim = c(0.95, 1.25), ylim = c(0, .20), 
      xlab  = "", col = vir[4])
axis (1, at = c(1, 1.1, 1.2), labels = c("Point", "AI_BF", "Pro_BF"), las = 2)
abline(emp_prophyp, 0, col = "black", lwd = 2)
points (x = 1.1, y =DD_bf_meanprop, pch = 19, col = vir[3])
points (x = 1.2, y = ProBF_hyp_mean, pch = 19, col = vir[2])

arrows(1, CIhigh_prophype_PM, 1, CIlow_prophype_PM, length = 0, lwd = 2, col = vir[4])
arrows(1.1, DD_bf_prop_ciHigh, 1.1, DD_bf_prop_cilow, length = 0, col = vir[3], lwd = 2)
arrows (1.2, ProBF_hyp_SE_high, 1.2, ProBF_hyp_SE_low, length = 0, col = vir[2], lwd = 2)
legend ("topleft", c("Empirical", "Point", "AI_bestFit", "Pro_bestFit"),
        col = c("black", vir[4], vir[3], vir[2]), lwd = c(2, 2, 2, 2), bty = "n")


### Log likelihood

#point mutation
PM <- SADISA_loglik(abund = abund, pars = c(692, 10e6), model = c("pm", "dl"), mult = "single")
PM #-2334.264

#less protracted
PRless <- SADISA_loglik(abund = abund, pars = c(782, 10e5, 10e6), model = c("pr", "dl"), mult = "single")
PRless #-2322.269

#more protracted
PRMore <- SADISA_loglik(abund = abund, pars = c(8180, 7000, 10e6), model = c("pr", "dl"), mult = "single")
PRMore # -6248.794

## Less abundance dependent: 
DDless <- SADISA_loglik(abund = abund, pars = c(237, 0.2, 10e6), model = c("dd", "dl"), mult = "single")
DDless # -2881.538

## more abundance dependent 
DDmore <- SADISA_loglik(abund = abund, pars = c(7, 0.5, 10e6), model = c("dd", "dl"), mult = "single")
DDmore # -5770.208

## Density-dependent with alpha -1 and best fit theta -- should be RF
RFdd <- SADISA_loglik(abund = abund, pars = c(5054, -1, 10e6), model = c("dd", "dl"), mult = "single")
RFdd # -6582.856

## Random fission
RF <- SADISA_loglik(abund = abund, pars = c(5048, 10e6), model = c("rf", "dl"), mult = single)
RF #  -6582.884 

## Protracted best fit parameters
PRBF <- SADISA_loglik(abund = abund, pars = c(7.425198e+02, 2.102163e+06, 1.000000e+07), model = c("pr", "dl"), mult = "single")
PRBF #-2305.203

# Density dependent best fit (negative alpha -- between point and RF)
DDBF <- SADISA_loglik(abund = abund, pars = c(815, -0.3, 10e6), model = c("dd", "dl"), mult = "single")
DDBF #-4298.475


########## Error robustness ##########

### Make new data
splt_abund25 <- numeric (length (ATDN_tots)+213)
a <- 1
for (ii in 1:213){
  qrt <- round ((ATDN_tots [ii]*.25), digits = 0)
  splt_abund25 [a] <- qrt
  splt_abund25 [a +1] <- ATDN_tots [ii] - qrt
  a <- a + 2
}
splt_abund25 [427:5240] <- ATDN_tots [214:5027]
splt_abund25 <- sort(splt_abund25, decreasing = T)

# mZSM
ATDN.mzsm <- fitmzsm (splt_abund25) ##
coef(ATDN.mzsm)
# theta = 726.6956

# Simulate 1000 times
my.theta <- 726.6956
totinds <- sum (splt_abund25)
nsims <- 1000
rand.neut_sim <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim <- append (rand.neut_sim, a)
  cat (ii, "\t")
}
saveRDS (rand.neut_sim)
rand.neut_sim2 <- rand.neut_sim
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
rand.neut.TS_2020_SPLIT_DF <- matrix(unlist(rand.neut_sim2),
                                     ncol = max.len,
                                     nrow = nsims,
                                     byrow = T )
rand.neut.TS_2020_SPLIT_DF[is.na(rand.neut.TS_2020_SPLIT_DF) ] <- 0
colnames(rand.neut.TS_2020_SPLIT_DF) <- 1:ncol(rand.neut.TS_2020_SPLIT_DF)
write.csv(rand.neut.TS_2020_SPLIT_DF)


rn_ATDN_SPLIT <- read.csv("rand.neut.ATDN_SPLIT_DF.csv")

CI_025 <- numeric (length (ncol (rn_ATDN_SPLIT)))
CI_975 <- numeric (length (ncol (rn_ATDN_SPLIT)))
mean_25R <- numeric (length (ncol (rn_ATDN_SPLIT)))
for (ii in 1:ncol (rn_ATDN_SPLIT)){
  temp_vec <- rn_ATDN_SPLIT [, ii]
  CI_025 [ii] <- quantile (temp_vec, 0.025)
  CI_975 [ii] <- quantile (temp_vec, 0.975)
  mean_25R [ii] <- mean (temp_vec)
}

sim_CI_RAD <- data.frame (empir = splt_abund25, mean_25R = mean_25R [1:length(splt_abund25)], CI_Low = CI_025[1:length (splt_abund25)],
                          CI_Up = CI_975 [1:length (splt_abund25)], 
                          row.names = 1:length (splt_abund25))

# Above or below 95 SE
comp.each.rank <- numeric (length (splt_abund25))
comp.value <- numeric (length (splt_abund25))
for (ii in 1:length (splt_abund25)){
  temp_rank <- splt_abund25[ii]
  if (temp_rank < CI_025[ii]){
    comp.each.rank[ii] <- "below"
    comp.value [ii] <- temp_rank - CI_025[ii]
  }else if (temp_rank > CI_975 [ii]){
    comp.each.rank [ii] <- "above"
    comp.value [ii] <- temp_rank - CI_975[ii]
  }else{
    comp.each.rank [ii] <- "in"
    comp.value [ii] <- 0
  }
}

sim_CI_RAD$comp <- comp.each.rank
sim_CI_RAD$compvalue <- comp.value
sim_CI_RAD$col <- ifelse(sim_CI_RAD$comp == "above", "#FF0000", 
                         ifelse(sim_CI_RAD$comp == "below", "#0000FF", "#000000"))

#plot RAD

plot (log (sim_CI_RAD$CI_Up), ylim = c(0, log (max (splt_abund25))), type = "l",
      xlab = "rank", ylab = "log abundance", lwd = 3)
lines (log (sim_CI_RAD$CI_Low), lwd = 3)
points (log(sim_CI_RAD$empir), 
        col = ifelse (sim_CI_RAD$compvalue >0, "red", ifelse(sim_CI_RAD$compvalue <0, "blue","black")),
        pch = 20, xlab = "empirical abundance", ylab = "difference from mean sim")

plot (sim_CI_RAD$CI_Up, ylim = c(0,  (max (splt_abund25))), type = "l",
      xlab = "rank", ylab = "log abundance", lwd = 3, xlim = c(0, 100))
lines (sim_CI_RAD$CI_Low, lwd = 3)
points (sim_CI_RAD$empir, 
        col = ifelse (sim_CI_RAD$compvalue >0, "red", ifelse(sim_CI_RAD$compvalue <0, "blue","black")),
        pch = 20, xlab = "empirical abundance", ylab = "difference from mean sim")


## Split top 5 ranks
hlf_abund5 <- numeric (length (ATDN_tots)+5)
a <- 1
for (ii in 1:5){
  hlf <- round ((ATDN_tots [ii]*.5), digits = 0)
  hlf_abund5 [a] <- hlf
  hlf_abund5 [a +1] <- ATDN_tots [ii] - hlf
  a <- a + 2
}

hlf_abund5 [11:5032] <- ATDN_tots [6:5027]
hlf_abund5 <- sort(hlf_abund5, decreasing = T)
plot (hlf_abund)

ATDN.mzsm <- fitmzsm (hlf_abund5) 
coef(ATDN.mzsm)
# theta = 693.2148

# Simulate 1000 times
my.theta <- 693.2148
totinds <- sum (hlf_abund)
nsims <- 1000
rand.neut_sim <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim <- append (rand.neut_sim, a)
  cat (ii, "\t")
}

rand.neut_sim2 <- rand.neut_sim
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
rand.neut.TS_2020_DF <- matrix(unlist(rand.neut_sim2),
                               ncol = max.len,
                               nrow = nsims,
                               byrow = T )
rand.neut.TS_2020_DF[is.na(rand.neut.TS_2020_DF) ] <- 0
colnames(rand.neut.TS_2020_DF) <- 1:ncol(rand.neut.TS_2020_DF)

write.csv(rand.neut.TS_2020_DF, file = "rand.neut.ATDN_Half5_DF.csv",row.names = F)

rn_ATDN_half5 <- read.csv("rand.neut.ATDN_Half5_DF.csv", stringsAsFactors = F)
max_abund <- max (rn_ATDN_half5[,1])

CI_025 <- numeric (length (ncol (rn_ATDN_half5)))
CI_975 <- numeric (length (ncol (rn_ATDN_half5)))
mean5 <- numeric (length (ncol (rn_ATDN_half5)))
for (ii in 1:ncol (rn_ATDN_half5)){
  temp_vec <- rn_ATDN_half5 [, ii]
  CI_025 [ii] <- quantile (temp_vec, 0.025)
  CI_975 [ii] <- quantile (temp_vec, 0.975)
  mean5 [ii] <- mean (temp_vec)
}

sim_CI_RAD5 <- data.frame (empir = hlf_abund5, CI_Low = CI_025[1:length (hlf_abund5)],
                           CI_Up = CI_975 [1:length (hlf_abund5)], mean5 = mean5 [1:length (hlf_abund5)],
                           row.names = 1:length (hlf_abund5))

comp.each.rank <- numeric (length (hlf_abund5))
comp.value <- numeric (length (hlf_abund5))
for (ii in 1:length (hlf_abund5)){
  temp_rank <- hlf_abund5[ii]
  if (temp_rank < CI_025[ii]){
    comp.each.rank[ii] <- "below"
    comp.value [ii] <- temp_rank - CI_025[ii]
  }else if (temp_rank > CI_975 [ii]){
    comp.each.rank [ii] <- "above"
    comp.value [ii] <- temp_rank - CI_975[ii]
  }else{
    comp.each.rank [ii] <- "in"
    comp.value [ii] <- 0
  }
}

sim_CI_RAD5$comp <- comp.each.rank
sim_CI_RAD5$compvalue <- comp.value
sim_CI_RAD5$col <- ifelse(sim_CI_RAD5$comp == "above", "#FF0000", 
                          ifelse(sim_CI_RAD5$comp == "below", "#0000FF", "#000000"))

## Split top 20 ranks
hlf_abund20 <- numeric (length (ATDN_tots)+20)
a <- 1
for (ii in 1:20){
  hlf <- round ((ATDN_tots [ii]*.5), digits = 0)
  hlf_abund20 [a] <- hlf
  hlf_abund20 [a +1] <- ATDN_tots [ii] - hlf
  a <- a + 2
}

hlf_abund20 [41:5047] <- ATDN_tots [21:5027]
hlf_abund20 <- sort(hlf_abund20, decreasing = T)
plot (hlf_abund20)

ATDN.mzsm <- fitmzsm (hlf_abund20) 
coef(ATDN.mzsm)
# theta = 695.7198

# Simulate 1000 times
my.theta <- 695.7198
totinds <- sum (hlf_abund)
nsims <- 1000
rand.neut_sim <- numeric ()
for (ii in 1:nsims) {
  a <- 0
  a <- list (rand.neutral(J = totinds, theta = my.theta))
  rand.neut_sim <- append (rand.neut_sim, a)
}


rand.neut_sim2 <- rand.neut_sim
max.len <- max(lengths(rand.neut_sim2))
rand.neut_sim2 <- lapply(rand.neut_sim2, `length<-`, max.len )
rand.neut.TS_2020_DF <- matrix(unlist(rand.neut_sim2),
                               ncol = max.len,
                               nrow = nsims,
                               byrow = T )
rand.neut.TS_2020_DF[is.na(rand.neut.TS_2020_DF) ] <- 0
colnames(rand.neut.TS_2020_DF) <- 1:ncol(rand.neut.TS_2020_DF)

write.csv(rand.neut.TS_2020_DF, file = "rand.neut.ATDN_Half20_DF.csv",
          row.names = F)

rn_ATDN_half20 <- read.csv("rand.neut.ATDN_Half20_DF.csv", stringsAsFactors = F)
max_abund <- max (rn_ATDN_half20[,1])

CI_025 <- numeric (length (ncol (rn_ATDN_half20)))
CI_975 <- numeric (length (ncol (rn_ATDN_half20)))
mean20 <- numeric (length (ncol (rn_ATDN_half20)))
for (ii in 1:ncol (rn_ATDN_half20)){
  temp_vec <- rn_ATDN_half20 [, ii]
  CI_025 [ii] <- quantile (temp_vec, 0.025)
  CI_975 [ii] <- quantile (temp_vec, 0.975)
  mean20 [ii] <- mean (temp_vec)
  
}

sim_CI_RAD20 <- data.frame (empir = hlf_abund20, CI_Low = CI_025[1:length (hlf_abund20)],
                            CI_Up = CI_975 [1:length (hlf_abund20)], mean20 = mean20 [1:length (hlf_abund20)],
                            row.names = 1:length (hlf_abund20))

comp.each.rank <- numeric (length (hlf_abund20))
comp.value <- numeric (length (hlf_abund20))
for (ii in 1:length (hlf_abund20)){
  temp_rank <- hlf_abund20[ii]
  if (temp_rank < CI_025[ii]){
    comp.each.rank[ii] <- "below"
    comp.value [ii] <- temp_rank - CI_025[ii]
  }else if (temp_rank > CI_975 [ii]){
    comp.each.rank [ii] <- "above"
    comp.value [ii] <- temp_rank - CI_975[ii]
  }else{
    comp.each.rank [ii] <- "in"
    comp.value [ii] <- 0
  }
}

sim_CI_RAD20$comp <- comp.each.rank
sim_CI_RAD20$compvalue <- comp.value
sim_CI_RAD20$col <- ifelse(sim_CI_RAD20$comp == "above", "#FF0000", 
                           ifelse(sim_CI_RAD20$comp == "below", "#0000FF", "#000000"))

## Plot 5 and 20
windows (6.8, 3.5)
par (mfrow = c(1, 2))
par (mar = c(3.2, 3.2, 2, 1))
plot (sim_CI_RAD5$mean5, type = "l",
      xlab = " ", ylab = " ", lwd = 3, log = "y", main = "Top 5 ranks split")
polygon(x = c(1:4905, rev(1:4905)),
        y = c(sim_CI_RAD5$CI_Up [1:4905], rev(sim_CI_RAD5$CI_Low [1:4905])),
        col =  adjustcolor("black", alpha.f = 0.20), border = NA)
points (sim_CI_RAD5$empir, 
        col = ifelse (sim_CI_RAD5$compvalue >0, "red", ifelse(sim_CI_RAD5$compvalue <0, "blue","black")),
        pch = 19)
legend ("bottomleft", legend = c("neutral+SE" ,"above neutral", "within neutral", "below neutral"), 
        col = c("grey", "red", "black", "blue"), lwd = c(7, 2, 2, 2), cex = 0.9, bty = "n")
legend ("bottomleft", legend = c(NA, NA, NA, NA), lwd = c(2, NA, NA, NA), cex = 0.9, bty = "n")
mtext (line = 2, side = 1, "Rank")
mtext (line = 2, side = 2, "Abundance")

plot (sim_CI_RAD20$mean20, type = "l",
      xlab = " ", ylab = " ", lwd = 3, log = "y", main = "Top 20 ranks split")
polygon(x = c(1:4905, rev(1:4905)),
        y = c(sim_CI_RAD20$CI_Up [1:4905], rev(sim_CI_RAD20$CI_Low [1:4905])),
        col =  adjustcolor("black", alpha.f = 0.20), border = NA)
points (sim_CI_RAD20$empir, 
        col = ifelse (sim_CI_RAD20$compvalue >0, "red", ifelse(sim_CI_RAD20$compvalue <0, "blue","black")),
        pch = 19)
mtext (line = 2, side = 1, "Rank")


########## Verification of rand.neutral ##########

## Compare SADs for coalescent simulations, rand.neutral simulations, and Alonso & McKane (2004) solution 

# Coalescent sims: JM = 3.9e6;  J = 10000, nu = 0.0001774359 (theta = 691.9998)
coalRad_10000 <- read.csv (file = "coal3.9e6_samp10000.csv")

## Rand.neutral: theta = 692, J = 10000
rnRad_10000 <- read.csv ("rnRad_10000.csv")

## A&M solution: J = 10000, JM = 3.9e6, theta = 692 (nu = 0.0001774359)
AM_J_10000 <- read.csv ("AM2004.J_SAD_10000.csv")


## Coalescent ##
## from RAD to SAD
coal_SAD_mat <- matrix(NA, nrow = 10000, ncol = max(coalRad_10000$X1))

for (ii in 1:nrow(coalRad_10000)){
  temp <- as.vector (t(coalRad_10000[ii,]))
  tempcount <- untb::count (temp)
  tempphi <- phi (tempcount, addnames = F)
  coal_SAD_mat [ii, 1:length(tempphi)] <- tempphi
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}

coal_SAD_mat_0 <- coal_SAD_mat
coal_SAD_mat_0 [is.na(coal_SAD_mat_0)] <- 0

coal_SAD_mean <- numeric (ncol (coal_SAD_mat_0))
coalSAD_CI_025 <- numeric (ncol (coal_SAD_mat_0))
coalSAD_CI_975 <- numeric (ncol (coal_SAD_mat_0))
for (ii in 1:ncol (coal_SAD_mat_0)){
  temp_vec <- coal_SAD_mat_0 [, ii]
  coal_SAD_mean [ii] <- mean (temp_vec, na.rm = T)
  coalSAD_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  coalSAD_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

## Rand.neutral RAD to SAD
rn_SAD_mat <- matrix(NA, nrow = 10000, ncol = max(rnRad_10000$X1))

for (ii in 1:nrow(rnRad_10000)){
  temp <- as.vector (t(rnRad_10000[ii,]))
  tempcount <- untb::count (temp)
  tempphi <- phi (tempcount, addnames = F)
  rn_SAD_mat [ii, 1:length(tempphi)] <- tempphi
  if (ii %% 100 == 0){
    cat(ii, "\t")
  }
}

rn_SAD_mat_0 <- rn_SAD_mat
rn_SAD_mat_0 [is.na(rn_SAD_mat_0)] <- 0

rn_SAD_mean <- numeric (ncol (rn_SAD_mat_0))
rnSAD_CI_025 <- numeric (ncol (rn_SAD_mat_0))
rnSAD_CI_975 <- numeric (ncol (rn_SAD_mat_0))
for (ii in 1:ncol (rn_SAD_mat_0)){
  temp_vec <- rn_SAD_mat_0 [, ii]
  rn_SAD_mean [ii] <- mean (temp_vec, na.rm = T)
  rnSAD_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  rnSAD_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

## Plot coalescent, rand.neutral, A&M
windows (5, 5)
par (mar = c(4, 4, 2, 1))
plot (coal_SAD_mean, type = "l", xlab = "Abundance", main = "SAD comparisons",
      ylab = "N species", log = "xy", lwd = 2, col = "forestgreen", ylim = c(1, 700), xlim = c(1, 30))
polygon(x = c(1:28, rev(1:28)),
        y = c(coalSAD_CI_025 [1:28], rev(coalSAD_CI_975 [1:28])),
        col =  adjustcolor("forestgreen", alpha.f = 0.30), border = NA)
lines (coalSAD_CI_025, col = "forestgreen")
lines (coalSAD_CI_975, col = "forestgreen")
lines (rn_SAD_mean, type = "l", lwd = 2, col = "orange")
polygon(x = c(1:28, rev(1:28)),
        y = c(rnSAD_CI_025 [1:28], rev(rnSAD_CI_975 [1:28])),
        col =  adjustcolor("orange", alpha.f = 0.30), border = NA)
lines (rnSAD_CI_975, col = "orange")
lines (rnSAD_CI_025, col = "orange")
lines (AM_J_10000$n.specs, type = "l", col = "blue", lwd = 2)
legend("bottomleft", c("coalescent", "rand.neutral", "A&M"), lty=1, col=c("forestgreen","orange", "blue"),
       lwd = 3)

## process v. sample noise ##

# Coalescent metacommunity simulations: JM = 3.9e6;  nu = 0.0001774359
coal_data <- read.csv(file  = "coal3_936_all.csv")

#sample noise
to_test <- seq(1:500)
sampList <- vector(mode = "list", length = 1)
sampsize <- 10000
samps <- 500
for (ii in to_test){
  tempdf <- data.frame (matrix (NA, nrow = samps, ncol = sampsize))
  colnames (tempdf) <- 1:ncol(tempdf)
  a <- 1
  tempvec <- as.vector (t(coal_data[ii,-1]))
  tempvec <- tempvec [!is.na(tempvec)]
  names (tempvec) <- paste ("spp", 1:length(tempvec), sep = "")
  tempcen <- census (tempvec)
  for (jj in 1: samps) {
    samp <-  sample (tempcen, sampsize, replace = F) 
    sampcount <- as.count (samp)
    sampcount <- sampcount [sampcount != 0]
    tempdf [a, 1:length(sampcount)] <- sampcount
    a <- a + 1
  }
  tcs <- colSums(tempdf, na.rm = T)  
  wer <- tcs [tcs == 0]
  end <- as.numeric (names (wer[1]))
  tempdf <- tempdf [, 1:end]
  sampList <- append (sampList, list(tempdf))
  cat (ii, "\t")
}

ll <- numeric (501)
for (ii in 2:501){
  ll[ii] <- ncol (sampList[[ii]])
}
max (ll) ## 2049

#mean and SE
samp_means <- data.frame (matrix (NA, nrow = 500, ncol = 2049))
samp_CI_low <- data.frame (matrix (NA, nrow = 500, ncol = 2049))
samp_CI_high <- data.frame (matrix (NA, nrow = 500, ncol = 2049))

for (ii in 2:501){
  tempdf <- sampList [[ii]]
  tempdf [is.na(tempdf)] <- 0
  
  for (jj in 1:ncol (tempdf)){
    temp_vec <- tempdf [, jj]
    samp_means [ii-1, jj] <- mean (temp_vec, na.rm = T)
    samp_CI_low [ii-1, jj] <- quantile (temp_vec, 0.025, na.rm = T)
    samp_CI_high [ii-1, jj] <- quantile (temp_vec, 0.975, na.rm = T)
  }
}

CI_dif <- samp_CI_high - samp_CI_low
CI_dif_prop_samp <- CI_dif / samp_means

### mean of sample noise
meansamp <- numeric (length (ncol (CI_dif)))
for (ii in 1:ncol (CI_dif)){
  meansamp[ii] <- mean (CI_dif [, ii])
}

## process&sample
sims <- 500
sampsize <- 10000
totSamps <- sims * sampsPersim
sampsPersim <- 500
coalRad_sAp <- data.frame (matrix (NA, nrow = totSamps, ncol = sampsize))

a <- 1
for (ii in 1:sims){
  tempvec <- as.vector (t(coal_data[ii,-1]))
  tempvec <- tempvec [!is.na(tempvec)]
  names (tempvec) <- paste ("spp", 1:length(tempvec), sep = "")
  tempcen <- census (tempvec)
  cat (ii, "\t")
  
  for (jj in 1: sampsPersim) {
    samp <-  sample (tempcen, sampsize, replace = F) 
    sampcount <- as.count (samp)
    sampcount <- sampcount [sampcount != 0]
    coalRad_sAp [a, 1:length(sampcount)] <- sampcount
    a <- a + 1
  }
}

coalRad_sAp0 <- coalRad_sAp [, 1:2047]
saveRDS(coalRad_sAp0, file = "coalRAD_samp&proc0RDS")

P_S <- readRDS(file = "coalRAD_samp&proc0RDS")
P_S [is.na(P_S)] <- 0

PS_mean <- numeric (ncol (P_S))
PS_CI_025 <- numeric (ncol (P_S))
PS_CI_975 <- numeric (ncol (P_S))
for (ii in 1:ncol (P_S)){
  temp_vec <- P_S [, ii]
  PS_mean [ii] <- mean (temp_vec, na.rm = T)
  PS_CI_025 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  PS_CI_975 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

var_PS_samp <- PS_CI_975 - PS_CI_025

## Rand.neutral with comparable parameters (same theta), J = 10000
rnRad_10000 <- read.csv ("rnRad_10000.csv")
rnRad_10000 [is.na(rnRad_10000)] <- 0

RN_mean_10000 <- numeric (ncol (rnRad_10000))
RN_CI_025_10000 <- numeric (ncol (rnRad_10000))
RN_CI_975_10000 <- numeric (ncol (rnRad_10000))
for (ii in 1:ncol (rnRad_10000)){
  temp_vec <- rnRad_10000 [, ii]
  RN_mean_10000 [ii] <- mean (temp_vec, na.rm = T)
  RN_CI_025_10000 [ii] <- quantile (temp_vec, 0.025, na.rm = T)
  RN_CI_975_10000 [ii] <- quantile (temp_vec, 0.975, na.rm = T)
}

RN_CIdif <- RN_CI_975_10000 - RN_CI_025_10000

# Plot all three together:
windows (5, 5)
par (mar = c(4, 4, 1, 1))
plot (x = 1:2049, y = CI_dif[1,], log = "xy", ylim = c(1, 60), xlim = c(1, 2050), xlab = "rank", 
      ylab = "SE_high - SE_low", type = "l")
for (ii in 2:nrow (CI_dif)){
  points (x = 1:2049, y = CI_dif[ii,], type = "l")
}
points (x = 1:2047, y = var_PS_samp, col = "red", type = "l", lwd = 2)
points (x = 1:2049, y = meansamp, col = "gray", type = "l", lwd = 2)
points (x = 1:2027, y = RN_CIdif, col = "blue", type = "l", lwd = 2)
legend ("topright", c("sample var", "mean sample var", "process + sample var", "rand.neutral var")
        , col = c("black", "gray", "red", "blue"), lwd = 2)


#### Modified SADISA simulation code ####
  # Modified from: https://github.com/cran/SADISA/blob/master/R/SADISA_sim.R

SADISA_sim2 <- function(parsmc,ii,jj,model = c('pm','dl'),mult = 'single', nsim = 1)
{
  if(!((model[1] == 'pm' | model[1] == 'pr' | model[1] == 'dd') & model[2] == 'dl'))
  {
    stop('Simulations for this metacommunity model is not implemented yet.');
  }
  if(!is.numeric(parsmc)){
    stop('Parameters should be scalars.')
  }
  out <- list();
  
  for(i in 1:nsim){
    out[[i]] <- list();
    if(mult == 'single'){
      ff <- checkiijj(ii,jj,mult); ii <- ff$ii; jj <- ff$jj; rm(ff);
      qq <- jj/(ii + jj);
      out[[i]] <- ms_sim2(parsmc,ii,qq,model);
    } else {
      stop('"mult" should be either "single", "ms", "mg", or "both"')
    }
    cat (i, "\t")
  }
  if(nsim == 1) {
    out <- out[[1]]
  }
  return(out);
}

ms_sim2 <- function(parsmc,ii,qq,model){
  if(model[1] == 'pm') {
    model_es0_int <- pm_estot_int;
  }else if(model[1] == 'pr') {
    model_es0_int <- pr_estot_int;
  }else if(model[1] == 'dd') {
    model_es0_int <- dd_estot_int;
  }
  
  # sample metacommunity abundances
  ff <- function(x)
  {
    return(exp(model_es0_int(x,parsmc,ii,qq)))
  }
  es0 <- stats::integrate(f = ff,lower = 0,upper = 1,rel.tol = 1e-9,abs.tol = 0)$value;
  nx <- stats::rpois(n = length(es0),es0);
  xs <- rep(0,nx);
  for(cnt in 1:nx)
  {
    prob <- stats::runif(1);
    ff2 <- function(x)
    {
      intresult <- rep(0,length(x))
      for(i in 1:length(x))
      {
        intresult[i] <- stats::integrate(ff,0,x[i])$value - prob * es0;
      }
      return(intresult);
    }
    #cat(paste('Compute xs[',cnt,']\n',sep = ''))
    xs[cnt] <- stats::uniroot(f = ff2,interval = c(0,1),tol = .Machine$double.eps)$root;
  }
  if(model[1] == 'dd' && parsmc[2] > 0) # this is alpha
  {
    xs <- xs^(1/(1 - parsmc[2]));
  }
  
  # sample local community abundances
  M <- length(ii);
  fasim <- matrix(0,nrow = nx, ncol = M);
  
  configs <- dec2binmat(M)
  for(ctr in 1:nx)
  {
    nn <- ii * xs[ctr]
    kk <- stats::rnbinom(n = M,size = nn, prob = 1 - qq);
    fasim[ctr,] <- kk;
    if(sum(kk) == 0)
    {
      p0 <- stats::dnbinom(x = 0, size = nn, prob = 1 - qq)
      p <- rep(0,2^M - 1)
      for(i in 2:(2^M))
      {
        p[i - 1] <- prod(p0^(1 - configs[i,]) * (1 - p0)^configs[i,])
      }
      config_possible <- which(c(0,p) != 0)
      if(length(config_possible) == 0)
      {
        totsam <- rowSums(configs)
        config_possible <- which(totsam == 1)
        p <- rep(1,length(config_possible))
      } else
      {
        p <- p[config_possible - 1]
      }
      config_sampled <- configs[DDD::sample2(config_possible,size = 1,prob = p),]
      samples_present <- which(config_sampled == 1)
      for(i in samples_present)
      {
        fele <- stats::runif(n = 1, min = stats::dnbinom(0, size = nn[i], prob = 1 - qq[i]), max = 1)
        fasim[ctr,i] <- stats::qnbinom(fele, size = nn[i], prob = 1 - qq[i])
        if(fasim[ctr,i] <= 0 | fasim[ctr,i] == Inf | is.na(fasim[ctr,i]) | is.nan(fasim[ctr,i]))
        {
          nmax <- 1000000
          probs <- stats::dnbinom(1:nmax, size = nn[i], prob = 1 - qq[i])
          which0 <- which(probs <= 0)
          if(length(which0) > 0)
          {
            nmax <- min(which0) - 1
          }
          fasim[ctr,i] <- DDD::sample2(x = 1:nmax, size = 1, prob = probs[1:nmax])
        }
      }
    }
  }
  fasim <- t(fasim);
  ff <- list();
  for(sam in 1:M)
  {
    ff[[sam]] <- fasim[sam,]
  }
  fasim <- ff;
  rm(ff);
  return(fasim);
}

## Other necessary functions

pm_estot_int <- function(x,parsmc,ii,qq)
{
  return(pmdlms_estot_int(x,c(parsmc,ii),qq));
}

dd_estot_int <- function(x,parsmc,ii,qq)
{
  return(mdd_lestot_int(x,c(parsmc,ii),qq))
}

pr_estot_int <- function(x,parsmc,ii,qq)
{
  th <- parsmc[1];
  ph <- parsmc[2];
  iiqq <- -sum(ii * log(1 - qq));
  y <- rep(0,length(x));
  for(cnt in 1:length(x))
  {
    if(x[cnt] > 0)
    {
      if(ph^2 *x[cnt]/(th + ph) > 100)
      {
        y[cnt] <- -th * ph * x[cnt]/(th + ph);
      } else
      {
        y[cnt] <- -ph * x[cnt] + log(expm1(ph^2 * x[cnt]/(th + ph)));
      }
      y[cnt] <- y[cnt] + log(-expm1(-iiqq * x[cnt])) - log(x[cnt]);
    } else
    {
      y[cnt] <- -Inf;
    }
  }
  y <- y + log(th);
  return(y);
}

#
# UTILS.R #
#' @title Computes integral of a very peaked function
#' @description   # computes the logarithm of the integral of exp(logfun) from 0 to Inf under the following assumptions:
# . exp(logfun) has a single, sharply peaked maximum
# . exp(logfun) is increasing to the left of the peak and decreasing to the right of the peak
# . exp(logfun) can be zero or positive at zero
# . exp(logfun) tends to zero at infinity
#' @param logfun the logarithm of the function to integrate
#' @param xx the initial set of points on which to evaluate the function
#' @param xcutoff when the maximum has been found among the xx, this parameter sets the width of the interval to find the maximum in
#' @param ycutoff set the threshold below which (on a log scale) the function is deemed negligible, i.e. that it does not contribute to the integral)
#' @param ymaxthreshold sets the deviation allowed in finding the maximum among the xx
#' @return the result of the integration
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

integral_peak <- function(logfun, xx = seq(-100,10,2), xcutoff = 2, ycutoff = 40, ymaxthreshold = 1E-12)
{
  # 1/ determine integrand peak
  yy <- xx + logfun(exp(xx));
  yy[which(is.na(yy) | is.nan(yy))] <- -Inf;
  yymax <- max(yy);
  if(yymax == -Inf)
  {
    logQ <- -Inf;
    return(logQ);
  }
  iimax <- which(yy >= (yymax - ymaxthreshold));
  xlft <- xx[iimax[1]] - xcutoff;
  xrgt <- xx[iimax[length(iimax)]] + xcutoff;
  optfun <- function(x) x + logfun(exp(x));
  optres <- stats::optimize(f = optfun, interval = c(xlft,xrgt), maximum = TRUE, tol = 1e-10);
  xmax <- optres$maximum;
  ymax <- optres$objective;
  
  # 2/ determine peak width
  iilft <- which((xx < xmax) & (yy < (ymax - ycutoff)));
  if(length(iilft) == 0)
  {
    xlft <- xx[1] - xcutoff;
  } else
  {
    ilft <- iilft[length(iilft)];
    xlft <- xx[ilft];
  }
  iirgt <- which((xx > xmax) & (yy < (ymax - ycutoff)));
  if(length(iirgt) == 0)
  {
    xrgt <- xx[length(xx)] + xcutoff;
  } else
  {
    irgt <- iirgt[1];
    xrgt <- xx[irgt];
  }
  
  # 3/ compute integral
  intfun <- function(x)
  {
    #if(any(is.nan(logfun(exp(x)))))
    #{
    #   print(exp(x))
    #   print(logfun(exp(x)))
    #   print(logfun)
    #}
    return(exp((x + logfun(exp(x))) - ymax))
  }
  intres <- stats::integrate(f = intfun, lower = xlft, upper = xrgt, rel.tol = 1e-10, abs.tol = 1e-10);
  corrfact <- intres$value;
  logQ <- ymax + log(corrfact);
  return(logQ);
}

#' @title Converts different formats to represent multiple sample data
#' @description Converts the full abundance matrix into species frequencies
#' If S is the number of species and M is the number of samples, then fa is
#' the full abundance matrix of dimension S by M. The  for example
#' fa = [0 1 0;3 2 1;0 1 0] leads to sf = [0 1 0 2;3 2 1 1];
#' @param fa the full abundance matrix with species in rows and samples in columns
#' @return the sample frequency matrix
#' @references Haegeman, B. & R.S. Etienne (2017). A general sampling formula for community structure data. Methods in Ecology & Evolution. In press.
#' @export

convert_fa2sf <- function(fa)
{
  dfa <- dim(fa);
  S <- dfa[1];
  M <- dfa[2];
  aux <- unique(fa);
  da <- dim(aux)[1]
  freq <- rep(0,da);
  for(cnt in 1:da)
  {
    ref <- t(t(rep(1,S))) %*% (aux[cnt,]);
    eqs <- rowSums(fa == ref);
    freq[cnt] = sum(eqs == M);
  }
  sf <- cbind(aux,freq);
  colnames(sf) <- NULL;
  return(sf)
}

convert_sf2fa <- function(sf)
{
  notzero <- which(sf != 0)
  fa <- NULL
  for(i in length(notzero):1)
  {
    fa <- c(fa,rep(notzero[i],sf[notzero[i]]))
  }
  return(fa)
}

difgamln <- function(a,n)
{
  # computes gammaln(a+n)-gammaln(a) as a sum of logarithms
  # for large a the result is numerically more precise than
  # directly evaluating the difference of gammaln
  # . a = real number, or vector of real numbers
  # . n = natural number, or vector of natural numbers
  if(length(a) != length(n))
  {
    stop('Input arguments should have same length.')
  }
  if(min(a) < 0)
  {
    stop('The first argument of difgamln must be non-negative.')
  }
  if(min(n) < 0)
  {
    stop('The second argument of difgamln must be non-negative.')
  }
  cvec <- rep(0,length(a));
  for(ctr in 1:length(a))
  {
    if(a[ctr] < 1e5 && n[ctr] > 100)
    {
      cvec[ctr] <- lgamma(a[ctr] + n[ctr]) - lgamma(a[ctr]);
    } else
    {
      if(n[ctr] == 0)
      {
        cvec[ctr] <- 0;
      } else #n[ctr] > 0
      {
        cvec[ctr] <- sum(log(a[ctr] + (0:(n[ctr] - 1))));
      }
      #if(is.nan(cvec[ctr])) {print(a[ctr]); print(n[ctr])}
      #if(cvec[ctr] == -Inf) {print(a[ctr]); print(n[ctr])}
    }
  }
  return(cvec)
}

checkiijj <- function(ii,jj,mult)
{
  if(is.list(ii))
  {
    if(is.list(ii[[1]]))
    {
      stop('I parameters not correctly specified.')
    }
    ii <- unlist(ii);
  }
  if(is.list(jj))
  {
    if(is.list(jj[[1]]))
    {
      stop('Sample sizes not correctly specified.')
    }
    jj <- unlist(jj);
  }
  if(length(ii) <= 1 && mult != 'single')
  {
    stop('You need to specify more than one I value.')
  }
  if(length(jj) <= 1 && mult != 'single')
  {
    stop('You need to specify more than one sample size.')
  }
  if(length(ii) > 1 && mult == 'single')
  {
    stop('You need to specify only one I value.')
  }
  if(length(jj) > 1 && mult == 'single')
  {
    stop('You need to specify only one sample size.')
  }
  return(list(ii = ii,jj = jj))
}

dec2bin <- function(y,ly)
{
  stopifnot(length(y) == 1, mode(y) == 'numeric')
  q1 <- (y / 2) %/% 1
  r <- y - q1 * 2
  res <- c(r)
  while(q1 >= 1)
  {
    q2 <- (q1 / 2) %/% 1
    r <- q1 - q2 * 2
    q1 <- q2
    res <- c(r, res)
  }
  res <- c(rep(0,ly - length(res)),res)
  return(res)
}

dec2binmat <- function(y)
{
  numrows <- 2^y
  res <- matrix(0,numrows,y)
  for(i in 0:(numrows-1))
  {
    res[i + 1,] <- dec2bin(i,y)
  }
  return(res)
}

is.nonnegativewholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
  ans <- abs(x - DDD::roundn(x)) < tol & x >= 0
  if(any(is.na(ans)))
  {
    ans <- FALSE
  }
  return(all(ans))
}

logminexpm1approx <- function(iiqq,xcnt,be)
{
  if(abs(-iiqq * xcnt^be) < 1E-100)
  {
    ans <- log(iiqq) + be * log(xcnt)
  } else
  {
    ans <- log(-expm1(-iiqq * xcnt^be))
  }
  return(ans)
}


## Other things required 
pmdlms_estot <- function(pars,qq)
{
  logfun <- function(x) pmdlms_estot_int(x,pars,qq);
  lestot <- integral_peak(logfun);
  estot <- exp(lestot);
  return(estot);
}



pmdlms_estot_int <- function(x,pars,qq)
{
  th <- pars[1];
  ii <- pars[2:length(pars)];
  iiqq <- -sum(ii * log(1 - qq));
  y <- rep(0,length(x));
  for(cnt in 1:length(x))
  {
    if(x[cnt] > 0)
    {
      y[cnt] <- log(-expm1(-iiqq * x[cnt])) + log(th) - th * x[cnt] - log(x[cnt]);
    } else
    {
      y[cnt] <- log(iiqq) + log(th);
    }
  }
  return(y);
}

mdd_lesk <- function(pars,qq,k)
{
  logfun <- function(x) mdd_lesk_int(x,pars,qq,k)
  lesk <- integral_peak(logfun);
  return(lesk);
}

mdd_estot <- function(pars,qq)
{
  logfun <- function(x) mdd_lestot_int(x,pars,qq)
  lestot <- integral_peak(logfun);
  estot <- exp(lestot);
  return(estot);
}

mdd_lesk_int <- function(x,pars,qq,k)
{
  th <- pars[1];
  al <- pars[2];
  ii <- pars[3];
  y <- rep(0,length(x));
  if(al > 0)
  {
    be <- 1/(1 - al) #exponent of substitution
    for(cnt in 1:length(x))
    {
      if(k > 1)
      {
        y[cnt] <- sum(log(ii * x[cnt]^be + (1:(k - 1))));
        #if(any(is.nan(y))) { print(k); print(x); print(y)}
      } else
      {
        y[cnt] <- 0;
      }
      y[cnt] <- y[cnt] + ii * x[cnt]^be * log(1 - qq) - th * x[cnt]^be;
    }
    y <- y - lgamma(2 - al);
  } else # al <= 0
  {
    for(cnt in 1:length(x))
    {
      if(k > 1)
      {
        y[cnt] <- sum(log(ii * x[cnt] + (1:(k - 1))));
      } else
      {
        y[cnt] <- 0;
      }
      y[cnt] <- y[cnt] + ii * x[cnt] * log(1 - qq) - al * log(x[cnt]) - th * x[cnt];
    }
    y <- y - lgamma(1 - al);
  }
  y <- y + (1 - al) * log(th) + k * log(qq) - lgamma(k + 1) + log(ii)
  return(y);
}

mdd_lestot_int <- function(x,pars,qq)
{
  th <- pars[1];
  al <- pars[2];
  ii <- pars[3:length(pars)];
  iiqq <- -sum(ii * log(1 - qq));
  y <- rep(0,length(x));
  if(al > 0)
  {
    be <- 1/(1 - al) #exponent of substitution
    for(cnt in 1:length(x))
    {
      if(x[cnt] > 0)
      {
        y[cnt] <- logminexpm1approx(iiqq,x[cnt],be) - be * log(x[cnt]) - th * x[cnt]^be;
        #y[cnt] <- log(-expm1(-iiqq * xcnt^be)) - be * log(x[cnt]) - th * x[cnt]^be;
      } else
      {
        y[cnt] <- log(iiqq);
      }
    }
    y <- y + (1 - al) * log(th) - lgamma(2 - al);
  } else
    if(al < 0)
    {
      for(cnt in 1:length(x))
      {
        if(x[cnt] > 0)
        {
          y[cnt] <- logminexpm1approx(iiqq,x[cnt],1) - (1 + al) * log(x[cnt]) - th * x[cnt];
          #y[cnt] <- log(-expm1(-iiqq * x[cnt])) - (1 + al) * log(x[cnt]) - th * x[cnt];
        } else
        {
          y[cnt] <- -Inf;
        }
      }
      y <- y + (1 - al) * log(th) - lgamma(1 - al);
    } else #al = 0 (pm)
    {
      for(cnt in 1:length(x))
      {
        if(x[cnt] > 0)
        {
          y[cnt] <- logminexpm1approx(iiqq,x[cnt],1) - log(x[cnt]) - th * x[cnt];
          #y[cnt] <- log(-expm1(-iiqq * x[cnt])) - log(x[cnt]) - th * x[cnt];
        } else
        {
          y[cnt] <- log(iiqq);
        }
      }
      y <- y + log(th);
    }
  return(y);
}


