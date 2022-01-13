#Simulation model of Farrior, Bohlman, Hubbell and Pacala, "Dominance of the Suppressed: Power law size structure in tropical forests."
#Last updated 15 Nov 2015   cfarrior@gmail.com

#To run, create a folder: "c:/usr/Farrior_etal_SleeperDist/"
#enter the following in R without the #'s: 
#setwd("c:/usr/Farrior_etal_SleeperDist/")
#source("Farrior_etal_SleeperSimulation.r")
#runV = c(filestem="PaperParameters",PA=1000,deltaT=1,Fnot=0.02,dnot=.02,G.c=6.05,G.u=.534,mu.c=0.0194/2,mu.u=0.0347-0.0194/2,mu=0.0194/2,cutT=400,landscape_m2=125000,censuscutT = 5)
#main_cohorts(runV=runV,tag=1,newplot=TRUE,plotting=TRUE)
#When plotting = TRUE, the simulation runs much more slowly but you can watch the patch develop and crash. 
#The best fit BCI power law (Fig. 1) is added for comparison. 
#*'s indicate empty bins and blue line is a smoothing of the data to incorporate those empty bins. 
#To finish, the size distribution for all recorded snapshots of the simulation is plotted in the same way. 

library( tidyverse )
library( ggthemes )

phi   <- 0.03615016
theta <- 1.2819275


runV  <- c(filestem="default",
           PA=5000,
           deltaT=1,
           Fnot=0.05,
           dnot=0.02,
           G.c=6.05,
           G.u=0.534,
           mu.c=0.0097,
           mu.u=0.025,
           mu=0.0097,
           cutT=400,
           landscape_m2=125000)

tag       <- 1
newplot   <- FALSE
plotting  <- TRUE


######################################
main_cohorts <- function(runV=c(filestem="default",
                                PA=5000,
                                deltaT=1,
                                Fnot=0.05,
                                dnot=0.02,
                                G.c=6.05,
                                G.u=0.534,
                                mu.c=0.0097,
                                mu.u=0.025,
                                mu=0.0097,
                                cutT=400,
                                landscape_m2=125000),
                         tag=1,
                         newplot=FALSE,
                         plotting=TRUE){
#runs the cohorts model
#runV: contains all of the inputs for the simulation, including:
	#filestem: the name to call the output files
	#PA: size in m2 of the plot
	#deltaT: timestep in years
	#Fnot: individuals produced per m2 of sun-lit crown area at dnot (reproduction to seedling)
	#dnot: the initial size of a seedling
	#G.c: the diameter growth rate, in mm/year of individuals in sun
	#G.u: the diameter growth, in mm/year of individuals in the understory
	#mu.c: the mortality rate in years^-1 of individuals in the sun
	#mu.u: the mortality rate in years^-2 of individuals in the understory
	#mu: the stand clearing disturbance rate
	#cutT: the time between recording snapshots of the forest - these used to build a landscape of independent forests
	#landscape_m2: the desired size of the final assembled landscape (landscape_m2 = PA*[total run time]/cutT)
#tag: an identifier for the specific run
#newplot: logical, whether to create a new plotting window
#plotting: logical, whether to plot the size distribution during the simulation, when TRUE slows the run speed considerably
######################################
  
	filestem <- runV[names(runV)=="filestem"][[1]]

	if(tag==1){
	  write.table(cbind(date(),t(runV)),
	              paste("main_cohorts_LOG.txt",sep=""),
	              append=TRUE,
	              sep="\t",
	              col.names=filestem=="start",
	              row.names=FALSE)
	} 
	
	# Not sure why values assigned as this?
	PA            <- as.numeric(runV[names(runV)=="PA"][[1]])
	deltaT        <- as.numeric(runV[names(runV)=="deltaT"][[1]])
	Fnot          <- as.numeric(runV[names(runV)=="Fnot"][[1]])
	dnot          <- as.numeric(runV[names(runV)=="dnot"][[1]])
	G.c           <- as.numeric(runV[names(runV)=="G.c"][[1]])
	G.u           <- as.numeric(runV[names(runV)=="G.u"][[1]])
	mu.c          <- as.numeric(runV[names(runV)=="mu.c"][[1]])
	mu.u          <- as.numeric(runV[names(runV)=="mu.u"][[1]])
	mu            <- as.numeric(runV[names(runV)=="mu"][[1]])
	cutT          <- as.numeric(runV[names(runV)=="cutT"][[1]])
	landscape_m2  <- as.numeric(runV[names(runV)=="landscape_m2"][[1]])
	
	# PA: size in m2 of the plot
	# deltaT: timestep in years
	# cutT: the time between recording snapshots of the forest - 
	  # these used to build a landscape of independent forests
	#landscape_m2: the desired size of the final assembled landscape 
	  # (landscape_m2 = PA*[total run time]/cutT)
	maxT <- landscape_m2 / PA*cutT

	# used to put a limit on seedlings if needed
	# dnot: the initial size of a seedling
	# (phi*dnot^theta) = CROWN AREA of individual of diameter dnot
	# maxF: maximum number of seedlings across area PA!
	maxF <- PA / (phi*dnot^theta)
	
	# ESTABLISH CORHORT (each row is a cohort)
	#columns: cohort diameter, #of individuals, crown class
	data     <- matrix(1,nrow=1,ncol=4) 
	data[,1] <- dnot; # cohort diameter
	data[,2] <- min(maxF, round(PA * Fnot * deltaT)); # of individuals
	data[,3] <- 1 # crown class
	data[,4] <- 1 # Age
	
	ndata    <- data #when the stand gets wiped out, it will come back to this initial setting
	
	# total number of iterations
	t_n      <- 2000
	  
	# track generation time
	gt_v     <- rep( NA, t_n)
	
	if(newplot) {win.graph(width=4,height=4); par(mfrow=c(1,1),oma=c(1,1,1,1))}
	
	for(t in 1:t_n ){ #seq(0,maxT,by=deltaT)
		# Step 1: Mortality
		# Step 1a: Stand clearing disturbance
	  # If stand-level disturbance, set trees back to ndata
	  print( t )
	  
		if( runif(1) < mu*deltaT ) data = ndata
		
		# Step 1b: background mortality
		for(i in seq(1,dim(data)[1])){
			if(data[i,3]==1) data[i,2] = rbinom(1,data[i,2],1-mu.c*deltaT) # if canopy
			if(data[i,3]==2) data[i,2] = rbinom(1,data[i,2],1-mu.u*deltaT) # if understory
		}
		
		# WHAT DOES THIS DO?!? Drop plots with 0 individuals?
		data                <- data[data[,2]>0,,drop=FALSE]
		
		#Step 2: Growth
		data[data[,3]==1,1] <- data[data[,3]==1,1] + G.c*deltaT # if canopy
		data[data[,3]==2,1] <- data[data[,3]==2,1] + G.u*deltaT # if understory
		
		# Drop time steps with average size == 0? 
		data                <- data[data[,1]>0,,drop=FALSE]

		# update the age of each cohort
		data[,4] <- data[,4] + 1
		
		#Step 3: Reproduce through dispersal from far away
		# ESTABLISH ANOTHER CORHORT
		babies <- round(PA*Fnot*deltaT) 
		
		# how could babies be == 0?
		if( babies > 0 ){
		  
			babyMatrix      <- matrix(1,nrow=1,ncol=4)
			# CA: Canopy Area?
			CA              <- sum( phi*data[,1]^theta * data[,2] )
			babyMatrix[,1]  <- dnot
			babyMatrix[,2]  <- babies
			# I ASSUME that if area of parents exceeds PA,
			# then new recruits are in the understory
			if(CA > PA) babyMatrix[,3] <- 2
			babyMatrix[,4]  <- 1 # age
			
			# Stack one cohort below another!
			data            <- rbind(data,babyMatrix)
			
		}
			
		#Step 4: Assign crown class
		CA                <- sum(phi*data[,1]^theta*data[,2])
		# If CA is < than PA, then all cohorts are in the canopy
		if(CA<=PA)        data[,3]=1
		# Understand how CCassign() works
		if(CA>PA)   data  <- CCassign(data,PA,deltaT)
		
		# 
		can_df <- subset( data, data[,3] == 1 )
		nom    <- sum( phi*can_df[,1]^theta*can_df[,2]* can_df[,4] )
		den    <- sum( phi*can_df[,1]^theta*can_df[,2] )
		gt_v[t]<- nom / den  
		
		#Step 5: Record and plot
		if(floor(t/cutT)==t/cutT){
			if(dim(data)[1]>0) write.table(data,
			                               paste(filestem,tag,".txt",sep=""),
			                               sep="\t",
			                               col.names=t==cutT,
			                               row.names=FALSE,
			                               append=t!=cutT)
		}
		if(plotting) if(max(data[,1])>10) bdata = SizeDistPlot(data,sizem2=PA,xmax=5000,ymin=1e-4)
	}
	
	tdata  <- read.table(paste(runV[[1]],tag,".txt",sep=""),sep="\t",header=TRUE)
	tbdata <- SizeDistPlot(as.matrix(tdata),sizem2=as.numeric(runV[names(runV)=="landscape_m2"][[1]]),main=paste(runV[[1]],tag))

}# end main_cohorts

######################################
CCassign <- function(data, PA, deltaT){
  
# main_cohorts function 
# assigns crown class based on the PPA assumption
# assumes all individuals have the same crown area and height allometries
# assumes that CAtot>PA 
######################################
  
  # Average area of individuals in each cohort 
	CAv             <- phi*data[,1]^theta 
	# Order based on avg. cohort area
	data            <- data[order(CAv,decreasing=TRUE),]
	CAv             <- CAv[order(CAv,decreasing=TRUE)]
	# Total area of each cohort (ordered!)
	cohortCAv       <- CAv*data[,2] 
	# Cumulative area
	cacaV           <- cumsum(cohortCAv)
	# Proportion of indiv. below canopy
	und             <- data[cacaV>PA, , drop=FALSE]
	# Proportion of indiv. in the canopy
	can             <- data[cacaV<PA, , drop=FALSE]
	
	# # when one cohort is so big that cacaV[1] > PA
	# # then put bigger indiv. in the canopy
	# if( nrow(can) == 0 ){
	#   can           <- und[1,]
	#   und           <- und[-1,]
	# }
	
	# total area occupied by the canopy
	canCA           <- max(0, sum( phi*can[,1]^theta*can[,2] ) )
	# assume only largest cohort in understory might belong to canopy
	tosplit         <- und[1,,drop=FALSE]
	opencan         <- PA - canCA # "leftover" canopy area
	# how many individuals in "tosplit" could "graduate" into the canopy?
	splitind_incan  <- floor( opencan/(phi*tosplit[1,1]^theta) )
	# subtract n. of individuals from the understory
	und[1,2]        <- und[1,2] - splitind_incan 
	# specify n. of individuals "graduated" to canopy
	tosplit[,2]     <- splitind_incan
	can             <- rbind(can,tosplit) # update canopy with "graduates"
	can[,3]         <- 1; # all of these belong to the canopy
	und[,3]         <- 2; # all of these belong to understory
	data            <- rbind(can, und)
	#always have a tree in the canopy, even if it's bigger than the plot area
	if(dim(can)[1]==0) data[1,3] <- 1 
	# remove "graduates" if there is 0 of them
	data            <- data[data[,2]>0,,drop=FALSE]
	return(data)
	
}# end CCassign


	
######################################
SizeDistPlot <- function(data,
                        logby=.14,
                        win=3,
                        sizem2=NaN,
                        main="",
                        census=5,
                        ploton=TRUE,
                        newplot=TRUE,
                        col=1,xaxt="s",yaxt="s",
                        ymin=1e-4,xmax=3500,
                        binV=NaN,justbdata=FALSE,
                        powerfit=list(alpha=2.1283,xmin=22,C=19.6211)
                        ){
#Minimal plotting function used in the main_cohorts simulations to plot the size distribution during the simulation run
#this function slows down the speed of the simulation quite a bit.
#data: either a vector of diameters of individuals in the community OR a matrix with the first column the diameter and the second column the number of individuals
#logby: the bin width in logspace for the diameter bins
#win: the number of bins on either side of the target bin to use to smooth over (win=3, makes a window size of 7 bins actually)
#sizem2: the size of the sampled area in m2
#main: header for the plot
#ploton: whether to plot
#newplot: whether to generate a new plotting window
#col: color to use to plot the data
#census: only needs setting if using BCI data with census <=1 where the binning of measurements was different
#xaxt and yaxt: whether to plot an xaxis and yaxis (for use with multiplot figures) 
#ymin: minimum value of y axis (individuals mm^-1 Ha^-1) for plotting
#xmax: maximum value of x (diameter) for plotting
######################################

	if(is.matrix(data)) dV = get_dV(data)
	if(is.vector(data)) dV = data

	dV = dV[order(dV)]
	dV = round(dV)
	dV = dV[dV>10]
	
	if(length(dV)==0) return(NaN)
	
	#make bin categories to group individuals by
	if(is.na(binV[1])){
		binV = exp(seq(log(min(dV)),log(max(dV))+logby*win,by=logby))
		#fix bin categories for censuses BCI censuses that were binned (census 0,1)
		if(census<=1) binV = c(seq(10,55,by=5),exp(seq(log(60),log(max(dV))+logby*win,by=logby)))
	}
	bdata = NULL
	for(i in seq(1,length(binV)-1)){
		d = dV[dV<binV[i+1]]
		d = d[d>=binV[i]]
		bdata = rbind(bdata,c(exp(mean(c(log(binV[i+1]),log(binV[i])))),length(d)/(binV[i+1]-binV[i])))
	}
	bdata[,2] = bdata[,2]/(sizem2)*10000 #generates individuals/mm/Ha
	if(justbdata) return(bdata)
	
	sbin = bdata[,1]
	if(win==0) sdata = bdata
	if(win>0){
		sdata = bdata[seq(1,win),]
		for(i in seq(win+1,length(sbin)-win)){
			n=NULL
			for(j in seq(i-win,i+win)){
				d = dV[dV<sbin[j+1]]
				d = d[d>=sbin[j]]
				d = d[!is.na(d)]
				n = c(n,length(d)/(sbin[j+1]-sbin[j]))
			}
			sdata = rbind(sdata,c(sbin[i],mean(n))) #this is not the mean in log space, because we need the zeros here
		}
	}	
	
	sdata[,2] = sdata[,2]/(sizem2)*10000 #generates individuals/mm/Ha
	
	if(ploton){
		if(newplot) plot(bdata[,1],bdata[,2],main=main,pch=20,xlab="Diameter (mm)",ylab=expression(paste("Individuals ( ",plain(mm^-1)," ", plain(Ha^-1)," )",sep="")),col=col,cex=1,log="xy",xlim=c(10,xmax),ylim=c(ymin,500),xaxt=xaxt,yaxt=yaxt)
	
		if(!newplot) points(bdata[,1],bdata[,2],pch=20,col=col,cex=1)
		zeros = bdata[bdata[,2]==0,,drop=FALSE]
		points(zeros[,1],exp(log(ymin)+log(2.5))+seq(1,dim(zeros)[1])*0,pch=8,cex=1,col=col)
		
		if(col==1) lines(sdata[,1],sdata[,2],col="blue",lwd=2)
		if(col!=1) lines(sdata[,1],sdata[,2],col=col,lwd=2)
		
		dd = seq(powerfit$xmin,max(dV))
		C = 1/(sum(dd^-powerfit$alpha))*length(dV[dV>powerfit$xmin])/sizem2*10000
		lines(dd,C*dd^-powerfit$alpha,lwd=1)
	}
	return(bdata)	
}#end SizeDistPlot


######################################
get_dV = function(data){
# function to turn data frame of cohorts into a list of diameters of individuals
# data is a matrix with columns : (1) diameter (2) number of individuals
######################################
	dV = NULL
	if(dim(data)[1]>0) for(i in seq(1,dim(data)[1])) dV = c(dV,data[i,1][[1]]+seq(1,data[i,2][[1]])*0)
	return(dV)
}#end get_dV


# Simulate average age --------------------------------------------------------

# simulate average age
sim_avg_age <- function( ii, t_n ){
  
  filestem <- runV[names(runV)=="filestem"][[1]]
  
  if(tag==1){
    write.table(cbind(date(),t(runV)),
                paste("main_cohorts_LOG.txt",sep=""),
                append=TRUE,
                sep="\t",
                col.names=filestem=="start",
                row.names=FALSE)
  } 
  
  # Not sure why values assigned as this?
  PA            <- as.numeric(runV[names(runV)=="PA"][[1]])
  deltaT        <- as.numeric(runV[names(runV)=="deltaT"][[1]])
  Fnot          <- as.numeric(runV[names(runV)=="Fnot"][[1]])
  dnot          <- as.numeric(runV[names(runV)=="dnot"][[1]])
  G.c           <- as.numeric(runV[names(runV)=="G.c"][[1]])
  G.u           <- as.numeric(runV[names(runV)=="G.u"][[1]])
  mu.c          <- as.numeric(runV[names(runV)=="mu.c"][[1]])
  mu.u          <- as.numeric(runV[names(runV)=="mu.u"][[1]])
  mu            <- as.numeric(runV[names(runV)=="mu"][[1]])
  cutT          <- as.numeric(runV[names(runV)=="cutT"][[1]])
  landscape_m2  <- as.numeric(runV[names(runV)=="landscape_m2"][[1]])
  
  # PA: size in m2 of the plot
  # deltaT: timestep in years
  # cutT: the time between recording snapshots of the forest - 
  # these used to build a landscape of independent forests
  #landscape_m2: the desired size of the final assembled landscape 
  # (landscape_m2 = PA*[total run time]/cutT)
  maxT <- landscape_m2 / PA*cutT
  
  # used to put a limit on seedlings if needed
  # dnot: the initial size of a seedling
  # (phi*dnot^theta) = CROWN AREA of individual of diameter dnot
  # maxF: maximum number of seedlings across area PA!
  maxF <- PA / (phi*dnot^theta)
  
  # ESTABLISH CORHORT (each row is a cohort)
  #columns: cohort diameter, #of individuals, crown class
  data     <- matrix(1,nrow=1,ncol=4) 
  data[,1] <- dnot; # cohort diameter
  data[,2] <- min(maxF, round(PA * Fnot * deltaT)); # of individuals
  data[,3] <- 1 # crown class
  data[,4] <- 1 # Age
  
  ndata    <- data #when the stand gets wiped out, it will come back to this initial setting
  
  # track generation time
  gt_ca_v  <- rep( NA, t_n)
  gt_ba_v  <- rep( NA, t_n)
  
  if(newplot) {win.graph(width=4,height=4); par(mfrow=c(1,1),oma=c(1,1,1,1))}
  
  for(t in 1:t_n ){ #seq(0,maxT,by=deltaT)
    # Step 1: Mortality
    # Step 1a: Stand clearing disturbance
    # If stand-level disturbance, set trees back to ndata
    print( t )
    
    # if( runif(1) < mu*deltaT ) data = ndata
    
    # Step 1b: background mortality
    for(i in seq(1,dim(data)[1])){
      if(data[i,3]==1) data[i,2] = rbinom(1,data[i,2],1-mu.c*deltaT) # if canopy
      if(data[i,3]==2) data[i,2] = rbinom(1,data[i,2],1-mu.u*deltaT) # if understory
    }
    
    # WHAT DOES THIS DO?!? Drop plots with 0 individuals?
    data                <- data[data[,2]>0,,drop=FALSE]
    
    #Step 2: Growth
    data[data[,3]==1,1] <- data[data[,3]==1,1] + G.c*deltaT # if canopy
    data[data[,3]==2,1] <- data[data[,3]==2,1] + G.u*deltaT # if understory
    
    # Drop time steps with average size == 0? 
    data                <- data[data[,1]>0,,drop=FALSE]
    
    # update the age of each cohort
    data[,4] <- data[,4] + 1
    
    #Step 3: Reproduce through dispersal from far away
    # ESTABLISH ANOTHER CORHORT
    babies <- round(PA*Fnot*deltaT) 
    
    # how could babies be == 0?
    if( babies > 0 ){
      
      babyMatrix      <- matrix(1,nrow=1,ncol=4)
      # CA: Canopy Area?
      CA              <- sum( phi*data[,1]^theta * data[,2] )
      babyMatrix[,1]  <- dnot
      babyMatrix[,2]  <- babies
      # I ASSUME that if area of parents exceeds PA,
      # then new recruits are in the understory
      if(CA > PA) babyMatrix[,3] <- 2
      babyMatrix[,4]  <- 1 # age
      
      # Stack one cohort below another!
      data            <- rbind(data,babyMatrix)
      
    }
    
    #Step 4: Assign crown class
    CA                <- sum(phi*data[,1]^theta*data[,2])
    # If CA is < than PA, then all cohorts are in the canopy
    if(CA<=PA)        data[,3]=1
    # Understand how CCassign() works
    if(CA>PA)   data  <- CCassign(data,PA,deltaT)
    
    # canopy
    can_df      <- subset( data, data[,3] == 1 )
    nom_ca      <- sum( phi*can_df[,1]^theta*can_df[,2]* can_df[,4] )
    den_ca      <- sum( phi*can_df[,1]^theta*can_df[,2] )
    nom_ba      <- sum( can_df[,1]*can_df[,2]* can_df[,4] )
    den_ba      <- sum( can_df[,1]*can_df[,2] )
    gt_ca_v[t]  <- nom_ca / den_ca  
    gt_ba_v[t]  <- nom_ba / den_ba
    
    #Step 5: Record and plot
    if(floor(t/cutT)==t/cutT){
      if(dim(data)[1]>0) write.table(data,
                                     paste(filestem,tag,".txt",sep=""),
                                     sep="\t",
                                     col.names=t==cutT,
                                     row.names=FALSE,
                                     append=t!=cutT)
    }
    if(plotting) if(max(data[,1])>10) bdata = SizeDistPlot(data,sizem2=PA,xmax=5000,ymin=1e-4)
    
  }

  data.frame( ca = gt_ca_v,
              ba = gt_ba_v,
              it = ii )
  
}

time  <- Sys.time()
prova <- sim_avg_age(1, 2000)
Sys.time() - time

# final simulation
sim_l1   <- lapply(1:5, sim_avg_age, 2000)
sim_l2   <- lapply(11:20, sim_avg_age, 2000)
sim_l3   <- lapply(21:30, sim_avg_age, 2000)
sim_l4   <- lapply(31:40, sim_avg_age, 2000)
sim_l5   <- lapply(41:50, sim_avg_age, 2000)
sim_l6   <- lapply(51:60, sim_avg_age, 2000)
sim_l7   <- lapply(61:70, sim_avg_age, 2000)
sim_l8   <- lapply(71:80, sim_avg_age, 2000)
sim_l9   <- lapply(81:90, sim_avg_age, 2000)
sim_l10  <- lapply(91:100, sim_avg_age, 2000)

# put it all together
sim_df   <-   Reduce( function(...) append(...),
                    list( sim_l1, sim_l2, sim_l3,
                          sim_l4, sim_l5, sim_l6,
                          sim_l7, sim_l8, sim_l9,
                          sim_l10 )  ) %>%
              bind_rows %>% 
              mutate( year = rep(1:2000, 100),
                      it   = as.factor(it) )

# store results
write.csv(sim_df, 'sim_gt_nodisturbance.csv', row.names=F)

sim_mean <- sim_df %>% 
              subset( !is.nan(ca) ) %>% 
              subset( !is.nan(ba) ) %>% 
              group_by( year ) %>% 
              summarise( ca_mean   = mean(ca),
                         ca_median = median(ca),
                         ba_mean   = mean(ba),
                         ba_median = median(ba) ) %>% 
              ungroup


ggplot(sim_df) +
  geom_line( aes( year, ca, 
                  color = it),
             alpha = 0.5 ) +
  geom_line( data = sim_mean,
             aes( year, ca_mean),
             lwd = 2, color = 'black',
             ) +
  theme_minimal() +
  theme( legend.position = 'none' ) +
  labs( x = "Time step",
        y = "Parents age (based on canopy area)" ) +
  ggsave( 'generation_time_PPA.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )

ggplot(sim_df) +
  geom_line( aes( year, ba, 
                  color = it),
             alpha = 0.5 ) +
  geom_line( data = sim_mean,
             aes( year, ba_mean),
             lwd = 2, color = 'black',
  ) +
  theme_minimal() +
  theme( legend.position = 'none' ) +
  labs( x = "Time step",
        y = "Parents age (based on basal area)" )


# Simulate maximum diameter ----------------------------------------------------

# simulate average age
sim_diam_distrib <- function( ii, t_n ){
  
  filestem <- runV[names(runV)=="filestem"][[1]]
  
  if(tag==1){
    write.table(cbind(date(),t(runV)),
                paste("main_cohorts_LOG.txt",sep=""),
                append=TRUE,
                sep="\t",
                col.names=filestem=="start",
                row.names=FALSE)
  } 
  
  # Not sure why values assigned as this?
  PA            <- as.numeric(runV[names(runV)=="PA"][[1]])
  deltaT        <- as.numeric(runV[names(runV)=="deltaT"][[1]])
  Fnot          <- as.numeric(runV[names(runV)=="Fnot"][[1]])
  dnot          <- as.numeric(runV[names(runV)=="dnot"][[1]])
  G.c           <- as.numeric(runV[names(runV)=="G.c"][[1]])
  G.u           <- as.numeric(runV[names(runV)=="G.u"][[1]])
  mu.c          <- as.numeric(runV[names(runV)=="mu.c"][[1]])
  mu.u          <- as.numeric(runV[names(runV)=="mu.u"][[1]])
  mu            <- as.numeric(runV[names(runV)=="mu"][[1]])
  cutT          <- as.numeric(runV[names(runV)=="cutT"][[1]])
  landscape_m2  <- as.numeric(runV[names(runV)=="landscape_m2"][[1]])
  
  # PA: size in m2 of the plot
  # deltaT: timestep in years
  # cutT: the time between recording snapshots of the forest - 
  # these used to build a landscape of independent forests
  #landscape_m2: the desired size of the final assembled landscape 
  # (landscape_m2 = PA*[total run time]/cutT)
  maxT <- landscape_m2 / PA*cutT
  
  # used to put a limit on seedlings if needed
  # dnot: the initial size of a seedling
  # (phi*dnot^theta) = CROWN AREA of individual of diameter dnot
  # maxF: maximum number of seedlings across area PA!
  maxF <- PA / (phi*dnot^theta)
  
  # ESTABLISH CORHORT (each row is a cohort)
  #columns: cohort diameter, #of individuals, crown class
  data     <- matrix(1,nrow=1,ncol=4) 
  data[,1] <- dnot; # cohort diameter
  data[,2] <- min(maxF, round(PA * Fnot * deltaT)); # of individuals
  data[,3] <- 1 # crown class
  data[,4] <- 1 # Age
  
  ndata    <- data #when the stand gets wiped out, it will come back to this initial setting
  
  # track generation time
  gt_ca_v  <- rep( NA, t_n)
  gt_ba_v  <- rep( NA, t_n)
  data_l   <- list()
  
  if(newplot) {win.graph(width=4,height=4); par(mfrow=c(1,1),oma=c(1,1,1,1))}
  
  for(t in 1:t_n ){ #seq(0,maxT,by=deltaT)
    # Step 1: Mortality
    # Step 1a: Stand clearing disturbance
    # If stand-level disturbance, set trees back to ndata
    print( t )
    
    # if( runif(1) < mu*deltaT ) data = ndata
    
    # Step 1b: background mortality
    for(i in seq(1,dim(data)[1])){
      if(data[i,3]==1) data[i,2] = rbinom(1,data[i,2],1-mu.c*deltaT) # if canopy
      if(data[i,3]==2) data[i,2] = rbinom(1,data[i,2],1-mu.u*deltaT) # if understory
    }
    
    # WHAT DOES THIS DO?!? Drop plots with 0 individuals?
    data                <- data[data[,2]>0,,drop=FALSE]
    
    #Step 2: Growth
    data[data[,3]==1,1] <- data[data[,3]==1,1] + G.c*deltaT # if canopy
    data[data[,3]==2,1] <- data[data[,3]==2,1] + G.u*deltaT # if understory
    
    # Drop time steps with average size == 0? 
    data                <- data[data[,1]>0,,drop=FALSE]
    
    # update the age of each cohort
    data[,4] <- data[,4] + 1
    
    #Step 3: Reproduce through dispersal from far away
    # ESTABLISH ANOTHER CORHORT
    babies <- round(PA*Fnot*deltaT) 
    
    # how could babies be == 0?
    if( babies > 0 ){
      
      babyMatrix      <- matrix(1,nrow=1,ncol=4)
      # CA: Canopy Area?
      CA              <- sum( phi*data[,1]^theta * data[,2] )
      babyMatrix[,1]  <- dnot
      babyMatrix[,2]  <- babies
      # I ASSUME that if area of parents exceeds PA,
      # then new recruits are in the understory
      if(CA > PA) babyMatrix[,3] <- 2
      babyMatrix[,4]  <- 1 # age
      
      # Stack one cohort below another!
      data            <- rbind(data,babyMatrix)
      
    }
    
    #Step 4: Assign crown class
    CA                <- sum(phi*data[,1]^theta*data[,2])
    # If CA is < than PA, then all cohorts are in the canopy
    if(CA<=PA)        data[,3]=1
    # Understand how CCassign() works
    if(CA>PA)   data  <- CCassign(data,PA,deltaT)
    
    # canopy
    can_df      <- subset( data, data[,3] == 1 )
    nom_ca      <- sum( phi*can_df[,1]^theta*can_df[,2]* can_df[,4] )
    den_ca      <- sum( phi*can_df[,1]^theta*can_df[,2] )
    nom_ba      <- sum( can_df[,1]*can_df[,2]* can_df[,4] )
    den_ba      <- sum( can_df[,1]*can_df[,2] )
    gt_ca_v[t]  <- nom_ca / den_ca  
    gt_ba_v[t]  <- nom_ba / den_ba
    
    #Step 5: Record and plot
    if(floor(t/cutT)==t/cutT){
      if(dim(data)[1]>0) write.table(data,
                                     paste(filestem,tag,".txt",sep=""),
                                     sep="\t",
                                     col.names=t==cutT,
                                     row.names=FALSE,
                                     append=t!=cutT)
    }
    if(plotting) if(max(data[,1])>10) bdata = SizeDistPlot(data,sizem2=PA,xmax=5000,ymin=1e-4)
    
    data_l[[t]] <- data
    
  }
  
  data_l
 
}

# final simulation
sim_l1   <- lapply(1:10,   sim_diam_distrib, 2000)
sim_l2   <- lapply(11:20,  sim_diam_distrib, 2000)
sim_l3   <- lapply(21:30,  sim_diam_distrib, 2000)
sim_l4   <- lapply(31:40,  sim_diam_distrib, 2000)
sim_l5   <- lapply(41:50,  sim_diam_distrib, 2000)
sim_l6   <- lapply(51:60,  sim_diam_distrib, 2000)
sim_l7   <- lapply(61:70,  sim_diam_distrib, 2000)
sim_l8   <- lapply(71:80,  sim_diam_distrib, 2000)
sim_l9   <- lapply(81:90,  sim_diam_distrib, 2000)
sim_l10  <- lapply(91:100, sim_diam_distrib, 2000)
sim_l    <- list( sim_l1, sim_l2, sim_l3, 
                  sim_l4, sim_l5, sim_l6, 
                  sim_l7, sim_l8, sim_l9,
                  sim_l10 ) %>% 
            Reduce( function(...) append(...), .)


# saveRDS( sim_l, 'C:/CODE/gwd_metrics/results/sim_l_nodisturbance.RDS')
sim_l <- readRDS( 'C:/CODE/gwd_metrics/results/farrior/sim_l_nodisturbance.RDS')
sim_l <- readRDS( 'C:/CODE/gwd_metrics/results/farrior/sim_l_disturbance.RDS')

# Format names of dataframes within this list
give_names <- function( x ){
  
  names_to_df <- function( x ){
    x %>% 
      as.data.frame %>% 
      setNames( c('diam_avg','N_indiv',
                  'layer','age') ) 
  }
  
  lapply( x, names_to_df )
  
}

# provide names
sim_l <- lapply(sim_l, give_names)

sim_l[[50]][[1500]] %>% head

  

# Calculate new ----------------------------------------------------------------


# give names by making list a data frame
give_names <- function(x){
  x %>% as.data.frame %>% 
    setNames( c('diam','n','layer','age') )
}

# calculate maximum diameter (as 95th percentile)
get_max_diam <- function(x, burnin){
  
  x <- lapply(x[burnin:2000], give_names)
  
  x %>% 
    bind_rows %>% 
    .$diam %>% 
    # 95th percentile of largest diameter
    quantile(prob=0.95)
  
}

# calculate maximum diameter
max_diam_v <- sapply(sim_l, get_max_diam, 1)
max_diam   <- mean(max_diam_v)


# Recalculate generation time --------------------------------------------------

# Calc. generation time based on max diameter
gen_time_max_diam <- function(ii, burnin, max_diam, sim_l){
  
  # provide names to each stable stage distribution
  x           <- sim_l[[ii]][burnin:2000] %>% 
                  lapply( give_names )
  
  # generation time at time t
  gen_time_at_t <- function(x, max_diam){
    
    demo        <- x %>% subset( diam > (max_diam/2) )
    nom_ba      <- sum( demo[,1] * demo[,2] * demo[,4] )
    den_ba      <- sum( demo[,1] * demo[,2] )
    
    nom_ba / den_ba
    
  }
  
  out <- sapply( x, gen_time_at_t, max_diam)
  out <- replace( out, is.nan(out), 0 )
  
  data.frame( year = c(burnin:2000),
              gt   = out,
              it   = ii )
  
}

# Calc. generation time based on canopy cohorts
gen_time_canopy <- function(ii, burnin, sim_l){
  
  # provide names to each stable stage distribution
  x           <- sim_l[[ii]][burnin:2000] %>% 
                   lapply( give_names )
  
  # generation time at time t
  gen_time_at_t <- function(x, max_diam){
    
    demo        <- x %>% subset( layer == 1 )
    nom_ba      <- sum( demo[,1] * demo[,2] * demo[,4] )
    den_ba      <- sum( demo[,1] * demo[,2] )
    
    nom_ba / den_ba
    
  }
  
  out <- sapply( x, gen_time_at_t, max_diam)
  out <- replace( out, is.nan(out), 0 )
  
  data.frame( year = c(burnin:2000),
              gt   = out,
              it   = ii )
  
}

# 100 time seris of generation times
gt_maxdiam_l  <- lapply(1:100, gen_time_max_diam, 1, max_diam, sim_l )
gt_canopy_l   <- lapply(1:100, gen_time_canopy,   1, sim_l )
gt_maxdiam_df <- gt_maxdiam_l %>% bind_rows 
gt_canopy_df  <- gt_canopy_l %>% bind_rows 

# calculate means across simulations (canopy estimation)
gt_canopy_mean  <- gt_canopy_df %>% 
                    group_by( year ) %>% 
                    summarise( gt = mean(gt) ) %>% 
                    ungroup %>% 
                    as.data.frame %>% 
                    mutate( Calculation = 'Canopy' )

# calculate means across simulations (max diam estimation)
gt_maxdiam_mean <- gt_maxdiam_df %>% 
                    group_by( year ) %>% 
                    summarise( gt = mean(gt) ) %>% 
                    ungroup %>% 
                    as.data.frame %>% 
                    mutate( Calculation = 'Max diam.' )

# plot means
list( gt_canopy_mean, 
      gt_maxdiam_mean ) %>% 
  bind_rows %>% 
  ggplot() +
  geom_line( aes(year, gt,
                 color = Calculation ) ) +
  scale_color_colorblind() +
  theme_minimal() +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 15) ) +
  labs( x = 'Year',
        y = 'Generation time' ) +
  ggsave( 'results/genT_canopy_vs_maxdiam_nodisturbance.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )


ggplot( gt_canopy_df ) +
  geom_line( aes( year, gt, 
                  color = as.factor(it)),
             alpha = 0.5 ) +
  geom_line( data = gt_canopy_mean,
             aes( year, gt ),
             lwd = 2 ) +
  theme_minimal() +
  theme( legend.position = 'none',
         axis.title = element_text( size = 20),
         axis.text  = element_text( size = 15) ) +
  labs( x = "Time step",
        y = "Parents age (based on diameters)" ) 
  ggsave( 'results/genT_canopy_nodisturbance.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )


ggplot( gt_maxdiam_df ) +
  geom_line( aes( year, gt, 
                  color = as.factor(it)),
             alpha = 0.5 ) +
  geom_line( data = gt_maxdiam_mean,
             aes( year, gt ),
             lwd = 2 ) +
  theme_minimal() +
  theme( legend.position = 'none',
         axis.title = element_text( size = 20),
         axis.text  = element_text( size = 15) ) +
  labs( x = "Time step",
        y = "Parents age (based on diameters)" ) 
  ggsave( 'results/genT_maxdiam_nodisturbance.tiff',
          width = 6.3, height = 6.3, compression = 'lzw' )
