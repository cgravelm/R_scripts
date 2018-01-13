library(raster)
library(ggplot2)
library(reshape)

# This parses the model outputs and adds up the multiple paths to the given map.

###########################################
## FIRST PART: GATHER GENERAL STATISTICS ##
###########################################

# Will collect some data that may be of importance.
results <- data.frame(run.n=double(0),
                     switchbacks=double(0),
                     optimization=double(0),
                     viewshed=double(0),
                     dist=double(0),
                     ticks=double(0),
                     max.slope=double(0),
                     mean.slope=double(0),
                     std=double(0),
                     end.x=double(0),
                     end.y=double(0))


# Identifies the folder in which all the output files are located
my_dir = "" ## Provide the link to the folder containing the csv outputs

# Lists all the file names of the csv I need to look into.
all_files = list.files(path = my_dir, all.files = TRUE, full.names = TRUE, pattern = "\\.csv$")

# Identifies the number of files to analyze.
list_size = length(all_files)

##################
## PROGRESS BAR ##
##################

Sys.time()

## This gets the progress bar
total <- list_size
pb <- txtProgressBar(min = 0, max = total, style = 3)

###########################
## UPLOAD THE DEM RASTER ##
###########################

DEM <- raster("") ## Provide the link to the DEM raster ASCII file.

# Get the dimensions of the raster to collate the route data
dem.x <- dim(DEM)[2]
dem.y <- dim(DEM)[1]

# Creates a dataframe that will hold the number of time each cell was walked on.
routes <- data.frame(matrix(0,ncol = dem.y,nrow = dem.x))

for(l in 1:list_size[1]){
  file.name = strsplit(all_files[l],"/")
  file.name = unlist(file.name)
  name.size = length(file.name)
  
  # map name for storing the R value is new.file
  new.file = file.name[name.size]
  
  # Separate out the components of the map name
  filename.split = strsplit(new.file,"_")
  filename.split = unlist(filename.split)
  time.stamp = filename.split [6]
  opt = filename.split [3] 
  switch = filename.split[4]
  viewshed = filename.split[5]
  
  nrow_results = nrow(results)
  new_row = nrow_results + 1
  results[new_row,] <- as.numeric(0)
  
  results$run.n[new_row] = time.stamp
  
  full.name <- paste(my_dir,new.file,sep="")
  
  # This is not pretty but given the weird format of the output, it does the job
  db <- read.table(full.name, fill = TRUE, nrow = 6, stringsAsFactors = FALSE, sep = ",", quote="\"")
  db.2 <- db[6,]
  colnames(db.2) <- db[5,]
  
  results$switchbacks[new_row] <- db.2$switchbacks
  results$optimization[new_row] <- db.2$optimization
  results$viewshed[new_row] <- db.2$viewshed_threshold
  start.x <- db.2$`start-x`
  start.y <- db.2$`start-y`
  end.x <- db.2$`end-x`
  end.y <- db.2$`end-y`
  
  # Import the real data without headers
  ds <- read.table(full.name, fill = TRUE, skip = 20, stringsAsFactors = FALSE, sep = ",")
  
  colnames(ds) <- c("t.1","x","col.1","pen.1","t.2","y","col.2","pen.2","t.3","slope","col.3","pen.3","t.4","dist","col.4","pen.4")
  
  ds <- na.omit(ds)
  
  if(nrow(ds) > 1){
  
  results$ticks[new_row] <- nrow(ds)
  results$dist[new_row] <- max(ds$dist, na.rm = T)
  results$max.slope[new_row] <- max(abs(ds$slope), na.rm = T)
  results$mean.slope[new_row] <- mean(abs(ds$slope), na.rm = T)
  results$std[new_row] <- sd(ds$slope, na.rm = T)
  results$end.x[new_row] <- ds$x[nrow(ds)]
  results$end.y[new_row] <- ds$y[nrow(ds)]
  
  assign('results',results, envir = .GlobalEnv)
  
  # This part creates the vector line between the two sites.
  x.y <- ds[2:nrow(ds),c(2,6)]
  x.y <- unique(x.y)

  i = 1

  for (i in 1:nrow(x.y))
  {
    # This is just for the local routes, to keep it at 10 km. They move to a goal at 30km as-crow-flies, but I record only their movement up to 10km.
    if (ds$dist[i] < 10){
    path.x <- as.numeric(x.y[i,1])
    path.y <- as.numeric(x.y[i,2])

    routes[path.x,path.y] <- routes[path.x,path.y] + 1
    assign('routes',routes, envir = .GlobalEnv)
    }
   }
  }
  Sys.sleep(0.1)
  
  # update progress bar
  setTxtProgressBar(pb, l)
}

################################################
## USING ROUTES TO IDENTIFY MOST POPULAR PATH ##
################################################

df <- routes
df$observation <- 1:nrow(df) 
df.melt <- melt(df, id.vars = "observation")
df.melt.fin <- subset(df.melt, df.melt$value != 0)

dat <- as.data.frame(sapply(df.melt.fin, function(x) 
  gsub(paste("X", collapse = '|'), '', x)))

colnames(dat) <- c("long","lat","value")
dat <- na.omit(dat)
dat$value <- as.numeric(as.character(dat$value))

dat$value <- log(dat$value) # Logging the values allow them to differentiate better in the figure.

dat$value <- dat$value / (max(dat$value))

dat$x <- as.numeric(as.character(dat[,1]))
dat$y <- as.numeric(as.character(dat[,2]))

# To transform into a usable raster by GRASS (get those values from DEM raster)
# For the local DEM at 100m resolution
cellsize <- 99.999532194444
xllcorner <- -1546592.764849999920
yllcorner <- -859940.537545749918

# For the global DEM at 1km resolution
# cellsize <- 999.006046081081
# xllcorner <- -1735185.553420000011
# yllcorner <- -1095974.338819053955

mid <- round(cellsize / 2)

dat$x <- (dat$x * cellsize ) + xllcorner - mid # xmin extent of the original map
dat$y <- (dat$y * cellsize ) + yllcorner - mid # ymin extent of the original map
dat <- dat[,c(4,5,3)]

# Identify the most used routes via stats

s <- summary(dat$value)
dat.sub <- subset(dat, dat$value > as.numeric(s[5])) # This preserves only the entries above the 3rd quartile

r.sub <- rasterFromXYZ(dat.sub)

# Create an ASCII raster from it.

writeRaster(r.sub,"", overwrite = T) ## Provide the link to the raster created (ASCII).

##########################
## SENSITIVITY ANALYSES ##
##########################

########################################
## Calibrating the parameter settings ##
########################################

# Calculate the statistical differences (t-test p-value) in the mean speed of simulations using different parameter settings.
vs.random <- read.csv(".csv", skip = 6, stringsAsFactors=FALSE) ## Provide the link to the BehaviorSpace output file
vs.random <- vs.random[,c(3,10,16,18,19,22,24,25)]

colnames(vs.random) <- c("view","opt","switch","dist","hours","crow","mean_slope","max_slope")
vs.random <- na.omit(vs.random)
vs.random$speed <- vs.random$dist/vs.random$hours

## Plot speed based on optimization, switchbacks, and viewshed

ggplot(results, aes(x = 1, y = as.numeric(speed), group = opt)) + 
  geom_boxplot(notch = T) + 
  facet_grid(opt~view, scales="free_y") +
  ggtitle("Speed - Viewshed\n") +
  xlab("") + ylab("Speed traveled (km/h)\n") +
  theme_bw(base_size = 16) %+replace% theme(strip.background  = element_blank()) +
  guides(fill=FALSE) +
  theme(axis.text.y=element_text(size=10), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

ggplot(results, aes(x = 1, y = as.numeric(speed), group = opt)) + 
  geom_boxplot(notch = T) + 
  facet_grid(opt~switch, scales="free_y") +
  ggtitle("Speed - Switchbacks\n") +
  xlab("") + ylab("Speed traveled (km/h)\n") +
  theme_bw(base_size = 16) %+replace% theme(strip.background  = element_blank()) +
  guides(fill=FALSE) +
  theme(axis.text.y=element_text(size=10), axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

# Test for difference significance

switch.x <- data.frame(optimization=double(0),
                       switch.a=double(0),
                       switch.b=double(0),
                       p.value=double(0))

# These can be changed to evaluate the impact of switchbacks and viewshed-threshold on the results.
ttst.sensitivity <- function(opt.pos, switch.pos)
{
  for (i in opt.pos)
  {
    
    for (j in switch.pos)
    {
      
      for (k in switch.pos)
      {  
        if (j != k)
        {
          
          nrow_res = nrow(switch.x)
          new_row = nrow_res + 1
          switch.x[new_row,] <- as.numeric(0)
          
          switch.x[new_row,]$optimization <- i
          switch.x[new_row,]$switch.a <- j
          switch.x[new_row,]$switch.b <- k
          sub.a <- subset(vs.random$speed, vs.random$opt == i & vs.random$view == j)
          sub.b <- subset(vs.random$speed, vs.random$opt == i & vs.random$view == k)
          
          p <- t.test(sub.a,sub.b)$p.value
          switch.x[new_row,]$p.value <- p
          assign('switch.x',switch.x, envir = .GlobalEnv)
          
        }
      }
    }
  }
}

ttst.sensitivity(c("\"Speed\"","\"Exploration\"","\"Distance\""), c(0, 0.2))

write.csv(switch.x, ".csv") ## Provide the link to the csv recording the t.test p-values.

################################
## Calibrating number of runs ##
################################

ttstplot <- function(data, variable)
{
  # This selects only the variable looked at here.
  b = match(variable, colnames(data))
  full_samp <- data[,b]
  
  # Creates a dataset to record the p-values
  ttst_mean <- nrow(300)
  ttst_sd <- nrow(300)
  
  for (i in 2:300)
  {
    ttst1 <- nrow(100)
    
    for (y in 1:100)
    {
      samp <- sample(full_samp, i, replace = FALSE)
      ttst1[y] <- t.test(full_samp, samp)$p.value 
    }
    ttst_mean[i] <- mean(ttst1, na.rm = TRUE)
    ttst_sd[i] <- sd(ttst1, na.rm = TRUE)
  }
  
  assign('ttst1',ttst1, envir = .GlobalEnv)
  
  ggplot() +
    geom_errorbar(aes(x=c(1:300), ymin=ttst_mean-ttst_sd, ymax=ttst_mean+ttst_sd), width=0.25) +
    geom_point(aes(x=c(1:300), y=ttst_mean)) +
    geom_segment(aes(x = 0, y = 0.05, xend = 300, yend = 0.05, colour = "red")) + 
    xlab("\nSample size") + ylab("t-test p-value\n") +
    theme_bw(base_size = 16) +
    guides(colour = "none")
}

vs.big <- read.csv(".csv", skip = 6, stringsAsFactors=FALSE) # Provide the link to the BehaviorSpace output file where 1000s of runs use the same parameters.
vs.big <- vs.big[,c(3,10,16,18,19,22,24,25)]
colnames(vs.big) <- c("view","opt","switch","dist","hours","crow","mean_slope","max_slope")

sub <- subset(vs.big, opt == "\"Exploration\"")
ttstplot(sub, "mean_slope")
