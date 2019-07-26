library(data.table)
library(adegenet)

#cladogenic
### This will pull out each of the unique cladogenetic events at each node from each stochastic map
load("BSM_output.RData")
clado_events_table = clado_events_tables[[2]]
clado_events_table = uniquify_clado_events(clado_events_table)
head(clado_events_table)

types <- unique(clado_events_table$clado_event_type)
col<-funky(4)

list_of_vicariance <- list()
list_of_sympatry <- list()
list_of_subset <- list()
list_of_founder <- list()

for(i in 1:50){
  clado_events_table = clado_events_tables[[i]]
  clado_events_table = uniquify_clado_events(clado_events_table)
  vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "vicariance (v)", ]
  list_of_vicariance[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
  vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "sympatry (y)", ]
  list_of_sympatry[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
  vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "subset (s)", ]
  list_of_subset[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
  vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == "founder (j)", ]
  list_of_founder[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
}

vicariance <- rbindlist((list_of_vicariance))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]
sympatry <- rbindlist((list_of_sympatry))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]
subset <- rbindlist((list_of_subset))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]
founder <- rbindlist((list_of_founder))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]

#anagenic
### This will pull out each of the unique anagenic events at each node from each stochastic map
ana_events_table = ana_events_tables[[20]]
ana_events_table = uniquify_clado_events(ana_events_table)
head(ana_events_table)

types <- unique(ana_events_table$event_type)
col<-funky(3)
list_of_d <- list()
list_of_e <- list()
for(i in 1:50){
  ana_events_table = ana_events_tables[[i]]
  vicariance_clado_events_table = ana_events_table[ana_events_table$event_type == "d", ]
  list_of_d[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
  vicariance_clado_events_table = ana_events_table[ana_events_table$event_type == "e", ]
  list_of_e[[i]] <- as.data.frame(table(round(vicariance_clado_events_table[,c("time_bp")])))
}

dispersal <- rbindlist((list_of_d))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]
list_of_e[[]] <- NULL #apparently no extinction?
extinction <- rbindlist((list_of_e))[,lapply(.SD,mean), list(as.numeric(as.vector(Var1)))]
head(dispersal)
 
vicariance <- as.data.frame(vicariance[order(as.numeric),])
sympatry <- as.data.frame(sympatry[order(sympatry$as.numeric),])
subset <- as.data.frame(subset[order(as.numeric),])
founder <- as.data.frame(founder[order(as.numeric),])
dispersal <- as.data.frame(dispersal[order(as.numeric),]) #anagenic

# add in 0s for events not observed at particular dates - change the "100" to whatever your oldest node is
names(sympatry) <- c("V1","V2")
sympatry <- rbind(sympatry,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which(seq(from = 0, to = 100, by = 1) %in% unique(sympatry$V1)==FALSE)],0)))
sympatry <- sympatry[order(sympatry$V1),]
names(vicariance) <- c("V1","V2")
vicariance <- rbind(vicariance,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which(seq(from = 0, to = 100, by = 1) %in% unique(vicariance$V1)==FALSE)],0)))
vicariance <- vicariance[order(vicariance$V1),]
names(subset) <- c("V1","V2")
subset <- rbind(subset,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which(seq(from = 0, to = 100, by = 1) %in% unique(subset$V1)==FALSE)],0)))
subset <- subset[order(subset$V1),]
names(founder) <- c("V1","V2")
founder <- rbind(founder,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which(seq(from = 0, to = 100, by = 1) %in% unique(founder$V1)==FALSE)],0)))
founder <- founder[order(founder$V1),]
names(dispersal) <- c("V1","V2")
dispersal <- rbind(dispersal,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which(seq(from = 0, to = 100, by = 1) %in% unique(dispersal$V1)==FALSE)],0)))
dispersal <- dispersal[order(dispersal$V1),]

#Calculate the total number of events at every age
all_events <- rbind(vicariance, sympatry, subset, founder, dispersal)
events_per_time <- as.data.frame(matrix(nrow=length(unique(all_events$V1)), ncol = 2))
events_per_time$V1 <- as.vector(unique(all_events$V1))
events_per_time <- events_per_time[order(events_per_time$V1),]
for(i in 1:nrow(events_per_time)){
  events_per_time[i,2] <- as.numeric(sum(all_events[which(all_events$V1==events_per_time[i,1]),"V2"]))
}
events_per_time


# Calculate the weighted frequency of each event type at each year (e.g., #founder/total # events)
vicariance$weighted_freq <- NA
sympatry$weighted_freq <- NA
subset$weighted_freq <- NA
founder$weighted_freq <- NA
dispersal$weighted_freq <- NA
vicariance <- as.data.frame(vicariance)
for(i in 1:nrow(vicariance)){
  vicariance[i,"weighted_freq"] <- as.numeric(vicariance[i,"V2"])/as.numeric(events_per_time[which(events_per_time$V1==vicariance[i,"V1"]),"V2"])
}
sympatry <- as.data.frame(sympatry)
for(i in 1:nrow(sympatry)){
  sympatry[i,"weighted_freq"] <- as.numeric(sympatry[i,"V2"])/as.numeric(events_per_time[which(events_per_time$V1==sympatry[i,"V1"]),"V2"])
}
subset <- as.data.frame(subset)
for(i in 1:nrow(subset)){
  subset[i,"weighted_freq"] <- as.numeric(subset[i,"V2"])/as.numeric(events_per_time[which(events_per_time$V1==subset[i,"V1"]),"V2"])
}
founder <- as.data.frame(founder)
for(i in 1:nrow(founder)){
  founder[i,"weighted_freq"] <- as.numeric(founder[i,"V2"])/as.numeric(events_per_time[which(events_per_time$V1==founder[i,"V1"]),"V2"])
}
dispersal <- as.data.frame(dispersal)
for(i in 1:nrow(dispersal)){
  dispersal[i,"weighted_freq"] <- as.numeric(dispersal[i,"V2"])/as.numeric(events_per_time[which(events_per_time$V1==dispersal[i,"V1"]),"V2"])
}
vicariance[which(is.na(vicariance$weighted_freq)==T), "weighted_freq"] <- 0
sympatry[which(is.na(sympatry$weighted_freq)==T), "weighted_freq"] <- 0
subset[which(is.na(subset$weighted_freq)==T), "weighted_freq"] <- 0
founder[which(is.na(founder$weighted_freq)==T), "weighted_freq"] <- 0
dispersal[which(is.na(dispersal$weighted_freq)==T), "weighted_freq"] <- 0




#### PLOTS ####
library(scales)
pdf("./Timing_of_events.pdf",width = 15, height = 12)
par(mar = c(5,5,2,5))
plot(x = sympatry$V1, y = sympatry$weighted_freq, main = "Timing of Events", type="l", xlab = "Time Before Present (MYA)", ylab = "Relative frequency of events", ylim = c(0, 1), col = alpha("red", 0.4), lwd=3)
par(new=T)
lines(jitter(x=subset$V1, amount = 0.3), y=subset$weighted_freq, col = alpha('blue', 0.4), lwd=3)
par(new=T)
lines(jitter(vicariance$V1,amount = 0.3), vicariance$weighted_freq, col = alpha("dark green", 0.4), lwd=3)
par(new=T)
lines(jitter(dispersal$V1,amount = 0.3), dispersal$weighted_freq, col = alpha("purple", 0.4), lwd = 3)
par(new=T)
lines(jitter(founder$V1,amount = 0.3), founder$weighted_freq, col = alpha("orange", 0.4), lwd = 3)
#text(48,1, "Sympatry", col = "red", adj = c(0, -.1))
#text(48,0.975, "Subset Sympatry", col = "blue", adj = c(0, -.1))
#text(48,0.95, "Vicariance", col = "dark green", adj = c(0, -.1))
#text(48,0.925, "Founder Event", col = "orange", adj = c(0, -.1))
#text(48,0.9, "Anagenic Dispersal", col = "purple", adj = c(0, -.1))
par(new=T)
plot(events_per_time$V1, events_per_time$V2, pch = 16, axes=FALSE, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 3, "Total number of events")
legend("topright", legend = c("Sympatry","Subset Sympatry","Vicariance","Founder Event","Anagenic Dispersal", "Total # Events"), 
       lty = c(1,1,1,1,1,0),
       lwd = c(2,2,2,2,2,0),
       col = c("red","blue","dark green","orange","purple","black"),
       pch = c(NA, NA, NA, NA, NA, 16))
dev.off()

##### Figure with timing of specific dispersal arrival events - i.e., how many times was each area the destination of a disperal at each year.

ana_events_table = ana_events_tables[[2]]
head(ana_events_table)
types <- unique(ana_events_table$ana_dispersal_from)

# Change the area codes based on whatever you used for your analysis.

list_of_A <- list()
list_of_B <- list()
list_of_C <- list()
list_of_D <- list()
list_of_E <- list()
list_of_F <- list()
list_of_G <- list()
list_of_H <- list()
for(i in 1:50){
  ana_events_table = ana_events_tables[[i]]
  temp = ana_events_table[ana_events_table$dispersal_to == "A", ]
  list_of_A[[i]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "B", ]
  list_of_B[[i]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "C", ]
  list_of_C[[i]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "D", ]
  list_of_D[[i]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "E", ]
  list_of_E[[i]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "F", ]
  list_of_F[[i]] <- as.data.frame(table(as.numeric(round(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "G", ]
  list_of_G[[i]] <- as.data.frame(table(as.numeric(round(temp[,c("time_bp")]))))
  temp = ana_events_table[ana_events_table$dispersal_to == "H", ]
  list_of_H[[i]] <- as.data.frame(table(as.numeric(round(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "A", ]
  list_of_A[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "B", ]
  list_of_B[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "C", ]
  list_of_C[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "D", ]
  list_of_D[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "E", ]
  list_of_E[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "F", ]
  list_of_F[[i+50]] <- as.data.frame(table(as.numeric(round(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "G", ]
  list_of_G[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
  temp = clado_events_table[clado_events_table$clado_dispersal_to == "H", ]
  list_of_H[[i+50]] <- as.data.frame(table(round(as.numeric(temp[,c("time_bp")]))))
}

A <- as.data.frame(rbindlist((list_of_A))[,lapply(.SD,mean), list(Var1)])
A$Var1 <- as.numeric(as.vector(A$Var1))
A <- A[order(A$Var1),]
B <- as.data.frame(rbindlist((list_of_B))[,lapply(.SD,mean), list(Var1)])
B$Var1 <- as.numeric(as.vector(B$Var1))
B <- B[order(B$Var1),]
C <- as.data.frame(rbindlist((list_of_C))[,lapply(.SD,mean), list(Var1)])
C$Var1 <- as.numeric(as.vector(C$Var1))
C <- C[order(C$Var1),]
D <- as.data.frame(rbindlist((list_of_D))[,lapply(.SD,mean), list(Var1)])
D$Var1 <- as.numeric(as.vector(D$Var1))
D <- D[order(D$Var1),]
E <- as.data.frame(rbindlist((list_of_E))[,lapply(.SD,mean), list(Var1)])
E$Var1 <- as.numeric(as.vector(E$Var1))
E <- E[order(E$Var1),]
F <- as.data.frame(rbindlist((list_of_F))[,lapply(.SD,mean), list(Var1)])
F$Var1 <- as.numeric(as.vector(F$Var1))
F <- F[order(F$Var1),]
G <- as.data.frame(rbindlist((list_of_G))[,lapply(.SD,mean), list(Var1)])
G$Var1 <- as.numeric(as.vector(G$Var1))
G <- G[order(G$Var1),]
H <- as.data.frame(rbindlist((list_of_H))[,lapply(.SD,mean), list(Var1)])
H$Var1 <- as.numeric(as.vector(H$Var1))
H <- H[order(H$Var1),]

#list_of_F <- list_of_F[-50]
#list_of_F <- list_of_F[-86]

names(A) <- c("V1","V2")
A <- rbind(A,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% A$V1)==FALSE)],0)))
A <- A[order(A$V1),]
names(B) <- c("V1","V2")
B <- rbind(B,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% B$V1)==FALSE)],0)))
B <- B[order(B$V1),]
names(C) <- c("V1","V2")
C <- rbind(C,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% C$V1)==FALSE)],0)))
C <- C[order(C$V1),]
names(D) <- c("V1","V2")
D <- rbind(D,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% D$V1)==FALSE)],0)))
D <- D[order(D$V1),]
names(E) <- c("V1","V2")
E <- rbind(E,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% E$V1)==FALSE)],0)))
E <- E[order(E$V1),]
names(F) <- c("V1","V2")
F <- rbind(F,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% F$V1)==FALSE)],0)))
F <- F[order(F$V1),]
names(G) <- c("V1","V2")
G <- rbind(G,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% G$V1)==FALSE)],0)))
G <- G[order(G$V1),]
names(H) <- c("V1","V2")
H <- rbind(H,as.data.frame(cbind(seq(from = 0, to = 100, by = 1)[which((seq(from = 0, to = 100, by = 1) %in% H$V1)==FALSE)],0)))
H <- H[order(H$V1),]

all_dispersal_events <- rbind(A, B, B, D, E, F, G, H)
all_dispersal_events <- all_dispersal_events[order(all_dispersal_events$V1),]
all_dispersal_events
dispersal_events_per_time <- as.data.frame(matrix(nrow=length(unique(all_dispersal_events$V1)), ncol = 2))
dispersal_events_per_time$V1 <- as.vector(unique(all_dispersal_events$V1))
dispersal_events_per_time <- dispersal_events_per_time[order(dispersal_events_per_time$V1),]
for(i in 1:nrow(dispersal_events_per_time)){
  dispersal_events_per_time[i,2] <- as.numeric(sum(all_dispersal_events[which(all_dispersal_events$V1==dispersal_events_per_time[i,1]),"V2"]))
}
dispersal_events_per_time

A$weighted_freq <- NA
B$weighted_freq <- NA
C$weighted_freq <- NA
D$weighted_freq <- NA
E$weighted_freq <- NA
F$weighted_freq <- NA
G$weighted_freq <- NA
H$weighted_freq <- NA

A <- as.data.frame(A)
for(i in 1:nrow(A)){
  A[i,"weighted_freq"] <- as.numeric(A[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==A[i,"V1"]),"V2"])
}
A[which(is.na(A$weighted_freq)==T), "weighted_freq"] <- 0
B <- as.data.frame(B)
for(i in 1:nrow(B)){
  B[i,"weighted_freq"] <- as.numeric(B[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==B[i,"V1"]),"V2"])
}
B[which(is.na(B$weighted_freq)==T), "weighted_freq"] <- 0
C <- as.data.frame(C)
for(i in 1:nrow(C)){
  C[i,"weighted_freq"] <- as.numeric(C[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==C[i,"V1"]),"V2"])
}
C[which(is.na(C$weighted_freq)==T), "weighted_freq"] <- 0
D <- as.data.frame(D)
for(i in 1:nrow(D)){
  D[i,"weighted_freq"] <- as.numeric(D[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==D[i,"V1"]),"V2"])
}
D[which(is.na(D$weighted_freq)==T), "weighted_freq"] <- 0
E <- as.data.frame(E)
for(i in 1:nrow(E)){
  E[i,"weighted_freq"] <- as.numeric(E[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==E[i,"V1"]),"V2"])
}
E[which(is.na(E$weighted_freq)==T), "weighted_freq"] <- 0
F <- as.data.frame(F)
for(i in 1:nrow(F)){
  F[i,"weighted_freq"] <- as.numeric(F[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==F[i,"V1"]),"V2"])
}
F[which(is.na(F$weighted_freq)==T), "weighted_freq"] <- 0
G <- as.data.frame(G)
for(i in 1:nrow(G)){
  G[i,"weighted_freq"] <- as.numeric(G[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==G[i,"V1"]),"V2"])
}
G[which(is.na(G$weighted_freq)==T), "weighted_freq"] <- 0
H <- as.data.frame(H)
for(i in 1:nrow(H)){
  H[i,"weighted_freq"] <- as.numeric(H[i,"V2"])/as.numeric(dispersal_events_per_time[which(dispersal_events_per_time$V1==H[i,"V1"]),"V2"])
}
H[which(is.na(H$weighted_freq)==T), "weighted_freq"] <- 0


####### PLOT #######
library(scales)
pdf("./Dispersal_destinations.pdf",width = 15, height = 12)
par(mar = c(5,5,2,5))
plot(x = A$V1, y = A$weighted_freq, main = "Dispersal Destinations", type="l", xlab = "Time Before Present (MYA)", ylab = "Relative frequency of dispersal destinations", ylim = c(0, 1), col = alpha("red", 0.4), lwd=3)
par(new=T)
lines(jitter(B$V1, amount = 0.3), B$weighted_freq, col = alpha('orange', 0.4), lwd=3)
par(new=T)
lines(jitter(C$V1,amount = 0.3), C$weighted_freq, col = alpha("purple", 0.4), lwd=3)
par(new=T)
lines(jitter(D$V1,amount = 0.3), D$weighted_freq, col = alpha("dark blue", 0.4), lwd=3)
par(new=T)
lines(jitter(E$V1,amount = 0.3), E$weighted_freq, col = alpha("black", 0.4), lwd = 3)
par(new=T)
lines(jitter(F$V1,amount = 0.3), F$weighted_freq, col = alpha("green", 0.4), lwd = 3)
par(new=T)
lines(jitter(G$V1,amount = 0.3), G$weighted_freq, col = alpha("yellow", 0.4), lwd = 3)
par(new=T)
lines(jitter(H$V1,amount = 0.3), H$weighted_freq, col = alpha("cyan", 0.4), lwd = 3)
#text(48,1, "Sympatry", col = "red", adj = c(0, -.1))
#text(48,0.975, "Subset Sympatry", col = "blue", adj = c(0, -.1))
#text(48,0.95, "Vicariance", col = "dark green", adj = c(0, -.1))
#text(48,0.9, "Anagenic Dispersal", col = "purple", adj = c(0, -.1))
par(new=T)
plot(dispersal_events_per_time$V1, dispersal_events_per_time$V2, pch = 16, axes=FALSE, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 3, "Total number of dispersal events")
legend("top", legend = c("Madagascar region","S Africa","Tropical Africa","N Africa/Eurasia","India", "SE Asia / Australasia", "N America","Neotropics", "Total # Events"), 
       lty = c(1,1,1,1,1,1,1,1,0),
       lwd = c(2,2,2,2,2,2,2,2,0),
       col = c("red","orange","purple","dark blue","black", "green","yellow","cyan","black"),
       pch = c(NA, NA, NA, NA, NA, NA, NA, NA, 16))
dev.off()


  
  
