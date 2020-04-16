## 3D GRP model that runs for a single lake, date and species.
## Also returns summary data file and 3D array data
## Updated to include Whitledge et al 2003 SMB parameters
## Date specific ShadW and shadED (prey) values for Acton and Hoover were
## based on weighted averages from gillnet catch.
## Original model code by G. Steinhart
## Last edited: 4/1/2020 by R. Budnik

##### Load Libraries #####
library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(PBSmapping)
library(maps)
library(data.table)
library(grid)
library(tcltk)

rm(list = ls())                             ## clear environment

#####  Menu and read files code  #####
lake <- tk_select.list(c("Acton", "Hoover"), title = "Pick a lake:")

## Read in shoreline polygons saved in reservoir data files folder
shoreline <- readOGR(paste("C:/Users/richa/OneDrive/Desktop/Reservoir data files/",
                                 lake, "/", lake, "_lake.shp", sep = ""))

## Read in depth data
depthdata.rg <- readOGR(paste("C:/Users/richa/OneDrive/Desktop/Reservoir data files/",
                             lake, sep = ""), paste(lake, "_depths", sep = ""))

## Read in lake acoustic data
allHAD <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir data files/", lake,
                         "/", lake, "_acoustics.csv", sep = ""))
allHAD$Date <- as.Date(allHAD$Date, "%m/%d/%Y")  ## changes date format

##  Read in water quality data
allWQ <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir data files/", lake,
                     "/", lake, "_WQ.csv", sep = ""))
allWQ$Trip_Date <- as.Date(allWQ$Trip_Date, "%m/%d/%Y")  ## changes date format

##  Read in gillnet data
allgillnet <- read.csv(paste("C:/Users/richa/OneDrive/Desktop/Reservoir data files/", lake,
                           "/", lake, "_gillnet.csv", sep = ""))
allgillnet$Trip_Date <- as.Date(allgillnet$Trip_Date, "%m/%d/%Y")  ## changes date format

##### Reads bioenergetic database #####
bioparam <- read.csv("C:/Users/richa/OneDrive/Desktop/Reservoir data files/Bioenergetics.csv")

## Input date for simulation and create 5-d window to include sampling not on exact day
sdates <- unique(allHAD$Date)
dayoichar <- tk_select.list(as.character(sdates), title = "Pick a date:")
dayoi <- as.Date(dayoichar)
day1 <- dayoi - 2
day3 <- dayoi + 2

## Subset data based on the date of interest, omits NA values from water quality data, and selects only correct gillnet dates
dHAD <- subset(allHAD, Date == dayoi)
dHAD$Easting <- as.numeric(as.character(dHAD$Easting))
dWQ <- subset(allWQ, Trip_Date >= day1 & Trip_Date <= day3,
                 select = c(UTMEasting, UTMNorthing, ReadDepth, Temp, DO))
dWQDO <- na.omit(dWQ)
dgillnet <- subset(allgillnet, Trip_Date >= day1 & Trip_Date <= day3)

##Fills in GN data from sampling during same months for dates when no data were available
if(lake == "Hoover" & dayoi == "2012-04-12"| dayoi == "2012-04-26"| dayoi == "2012-05-08"| dayoi == "2012-05-22") dgillnet <- subset(allgillnet, Trip_Date = 5/8/2006) 
if(lake == "Hoover" & dayoi == "2012-06-19") dgillnet <- subset(allgillnet, Trip_Date = 6/5/2006) 
if(lake == "Hoover" & dayoi == "2012-07-03") dgillnet <- subset(allgillnet, Trip_Date = 7/13/2006) 
if(lake == "Hoover" & dayoi == "2011-08-25"| dayoi =="2012-08-03"| dayoi == "2012-08-14"| dayoi == "2012-08-29") dgillnet <- subset(allgillnet, Trip_Date = 8/7/2006) 
if(lake == "Hoover" & dayoi == "2012-09-11"| dayoi == "2012-09-25") dgillnet <- subset(allgillnet, Trip_Date = 9/5/2006) 
if(lake == "Hoover" & dayoi == "2011-10-25" | dayoi == "2012-10-10"| dayoi == "2012-10-23") dgillnet <- subset(allgillnet, Trip_Date = 10/10/2006) 

## Process gillnet data based on gape size selected
gape <- tk_select.list(c("139 (LMB2)", "172 (LMB3)", "139 (SAE1)", "158 (SAE2)", "169 (SAE3)", 
                                  "138 (WAE2)", "149 (WAE3)"), title = "Select gape size (mm):")
if(gape == "139 (LMB2)") yoy <- sum(dgillnet$Length < 139) ## sets gape limit for avg size of Age-2 LMB from Michaletz 1997
if(gape == "172 (LMB3)") yoy <- sum(dgillnet$Length < 172) ## sets gape limit for avg size of Age-3 LMB from Michaletz 1997
if(gape == "139 (SAE1)") yoy <- sum(dgillnet$Length < 139) ## sets gape limit for avg size of Age-1 SAE from Knight et al. 1984 (based on WAE gape)
if(gape == "158 (SAE2)") yoy <- sum(dgillnet$Length < 158) ## sets gape limit for avg size of Age-2 SAE from Knight et al. 1984 (based on WAE gape)
if(gape == "169 (SAE3)") yoy <- sum(dgillnet$Length < 169) ## sets gape limit for avg size of Age-3 SAE from Knight et al. 1984 (based on WAE gape)
if(gape == "138 (WAE2)") yoy <- sum(dgillnet$Length < 138) ## sets gape limit for avg size of Age-2 WAE from Knight et al. 1984
if(gape == "149 (WAE3)") yoy <- sum(dgillnet$Length < 149) ## sets gape limit for avg size of Age-3 WAE from Knight et al. 1984

preyprop <- yoy / nrow(dgillnet)  ## calculates the proportion of GS that are YOY

## Process acoustic data
dHAD$cDensity <- preyprop * dHAD$Density    ## adjusts for proportion gape appropriate
dHAD$cBiomass <- preyprop * dHAD$Biomass

dHAD$cDensity <- dHAD$cDensity / 10000      ## converts #/ha to #/m3
dHAD$cBiomass <- dHAD$cBiomass * 0.1        ## converts kg/ha to g/m3

for(j in 0:max(dHAD$Layer_depth_min)) {     ## fills in vertical profile of abundance, 0-2 m depth filled when side-scan data available
    if(!exists("preydata")) {
        preydata <- subset(dHAD, Layer_depth_min == j, select = c(Layer_depth_min, 
                            Easting, Northing, cDensity, cBiomass))
    } else {
        if(j == 1) {
            temp_dataset <- preydata
            temp_dataset$Layer_depth_min = 1
            preydata <- rbind(preydata, temp_dataset)
            rm(temp_dataset)
        } else {
            temp_dataset <- subset(dHAD, Layer_depth_min == j, select = 
                                       c(Layer_depth_min, Easting, Northing, 
                                         cDensity, cBiomass))
            preydata <- rbind(preydata, temp_dataset)
            rm(temp_dataset)
        }
    }
}

preydata$lnBiomass <- log(preydata$cBiomass + 1)    ## compute ln(biomass)
colnames(preydata)[1] <- "Depth"

##  Define min and max values for plots
maxTEMP = 33
medTEMP = 20
minTEMP = 8

maxDO = 20
medDO = 4
minDO = 0

maxPREY = max(preydata$lnBiomass)
medPREY = median(preydata$lnBiomass)
minPREY = 0

nxsect = max(dHAD$Interval)
maxzWQ = max(dWQ$ReadDepth)
maxzHAD = max(dHAD$Layer_depth_min)
maxz = max(maxzWQ, maxzHAD)  ## Find maximum sample depths from WQ and HAD

depthsTEMP = paste("mT", 0:maxz, sep = "") ## make character names for areal plots (one for every m depth)
depthsDO = paste("mD", 0:maxz, sep = "")
depthsBIOM = paste("mB", 0:maxz, sep = "")


##### Subset and process depth data  #####
depthdata <- as.data.frame(depthdata.rg)
depthd <- as.data.frame(matrix(0, ncol = 3, nrow = nrow(depthdata)))
colnames(depthd) <- c("Easting", "Northing", "Depth")
depthd$Easting <- depthdata$Easting
depthd$Northing <- depthdata$Northing
if(lake != "Findlay1" & lake != "Findlay2" & lake != "Pleasant Hill") depthd$Depth <- -depthdata$Depth_m
if(lake == "Findlay1" | lake == "Findlay2" | lake == "Pleasant Hill") depthd$Depth <- -depthdata$depth_m

## determine if value = 0 and subset rows not equal to all 0's
row_sub = apply(depthd, 1, function(row) all(row !=0 ))
depthd <- depthd[row_sub,]

## subset for really large datasets then make depths positive, speeds up model, but less detail for depth
if(nrow(depthdata) >= 15000 & nrow(depthdata) < 50000) depthd <- depthd[seq(1, nrow(depthd), 10), ]
if(nrow(depthdata) >= 50000 & nrow(depthdata) < 200000) depthd <- depthd[seq(1, nrow(depthd), 20), ]
if(nrow(depthdata) >= 200000) depthd <- depthd[seq(1, nrow(depthd), 100), ]


##### Process Shoreline #####
## Change format to dataframe and define projection, if needed
if(lake == "Findlay1" | lake == "Findlay2" | lake == "Delaware") {
  proj4string(shoreline) <- CRS("+init=epsg:3728")
  shoreline <- spTransform(shoreline, CRS("+proj=utm +zone=17"))
}

if(lake != "Acton") proj4string(shoreline) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(shoreline) <- CRS("+proj=utm +zone=16") 
shoreoutline <-fortify(shoreline)


## Create zero depth data (i.e., shoreline) and add to depth data
zerodepth <- shoreoutline[,1:2]
zerodepth$Depth <- 0
colnames(zerodepth) <- c("Easting", "Northing", "Depth")
depthd <- rbind(depthd, zerodepth)
depthd3D <- depthd


## Set spatial coordinates to create a Spatial object
coordinates(depthd) <- ~Easting + Northing
if(lake != "Acton") proj4string(depthd) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(depthd) <- CRS("+proj=utm +zone=16")

## generates polygons from shoreline data, for making the areal maps and determining which cells are not land
if(lake == "Tappan") {  
    poly <- subset(shoreoutline, shoreoutline$group == 7.1, select = c(long, lat))
} else {
    poly <- subset(shoreoutline, shoreoutline$piece == 1, select = c(long, lat))
}
polyP <- Polygon(cbind(poly$long, poly$lat)) ## Add lat & lon to object
polySP <- SpatialPolygons(list(Polygons(list(polyP), ID=1)))


## Define grid extent for areal graphing and interpolating (grid is larger than reservoir, but land subtracted later)
xmax = round(max(depthd$Easting) + 100, -2)
xmin = round(min(depthd$Easting) - 100, -2)
ymax = round(max(depthd$Northing) + 100, -2)
ymin = round(min(depthd$Northing) - 100, -2)
x.range <- as.numeric(c(xmin, xmax))
y.range <- as.numeric(c(ymin, ymax))

## Create matrix of outside grid bounds and change to SpatialPolygons Dataframe
outside <- matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin),
                  ncol = 2, byrow = TRUE)
outside <- SpatialPolygons(list(Polygons(list(Polygon(outside)), ID = 1)))
outside <- SpatialPolygonsDataFrame(Sr=outside, data=data.frame(polyvar=357))

## Make an overlay grid for range from min to max by units of 50
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 50),
                   y = seq(from = y.range[1], to = y.range[2], by = 50)) 
coordinates(grd) <- ~x + y
if(lake != "Acton") proj4string(grd) <- CRS("+proj=utm +zone=17")
if(lake == "Acton") proj4string(grd) <- CRS("+proj=utm +zone=16")
gridded(grd) <- TRUE


##### Interpolate depth using inverse distance weighted model  #####
idwd <- idw(formula = Depth ~ 1, locations = depthd, newdata = grd)
idwd.output = as.data.frame(idwd)     ## Output defined as dataframe with names
names(idwd.output)[1:3] <- c("Easting", "Northing", "Depth")
idwd.output$Depth <- round(idwd.output$Depth)                                   ## Round depth data

## Generate interpolated (IDW method) temp, DO and biomass data for each depth
for(doi in 0:maxz) {
    ## subsets data based on UTM
    idwTEMP <- subset(dWQ, dWQ$ReadDepth == doi, select = c(UTMEasting, UTMNorthing, Temp))
    idwDO <- subset(dWQDO, dWQDO$ReadDepth == doi, select = c(UTMEasting, UTMNorthing, DO))
    idwHAD <- subset(preydata, preydata$Depth == doi, select = c(Easting, Northing, cDensity))
    idwHAB <- subset(preydata, preydata$Depth == doi, select = c(Easting, Northing, lnBiomass))
    
    if(doi > maxzWQ) {       ## duplicates previous depth layer for deeper depths
        idwTEMP <- subset(dWQ, dWQ$ReadDepth == maxzWQ, select = c(UTMEasting, UTMNorthing, Temp))
        idwDO <- subset(dWQDO, dWQDO$ReadDepth == maxzWQ, select = c(UTMEasting, UTMNorthing, DO))
    }

    if(doi > maxzHAD) {     ## duplicates previous depth layer for deeper depths
        idwHAD <- subset(preydata, preydata$Depth == maxzHAD, select = c(Easting, Northing, cDensity))
        idwHAB <- subset(preydata, preydata$Depth == maxzHAD, select = c(Easting, Northing, lnBiomass))
    }
    
    ## Set spatial coordinates to create Spatial objects and interpolate
    coordinates(idwTEMP) <- ~UTMEasting + UTMNorthing
    if(lake != "Acton") proj4string(idwTEMP) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwTEMP) <- CRS("+proj=utm +zone=16")
    idwt <- idw(formula = Temp ~ 1, locations = idwTEMP, newdata = grd)
    idwt.output = as.data.frame(idwt)   ## Output defined as dataframe with names
    names(idwt.output)[1:3] <- c("Easting", "Northing", "Temp")
    
    coordinates(idwDO) <- ~UTMEasting + UTMNorthing
    if(lake != "Acton") proj4string(idwDO) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwDO) <- CRS("+proj=utm +zone=16")
    idwo <- idw(formula = DO ~ 1, locations = idwDO, newdata = grd)
    idwo.output = as.data.frame(idwo)   ## Output defined as dataframe with names
    names(idwo.output)[1:3] <- c("Easting", "Northing", "DO")
    
    coordinates(idwHAD) <- ~Easting + Northing
    if(lake != "Acton") proj4string(idwHAD) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwHAD) <- CRS("+proj=utm +zone=16")
    idwDEN <- idw(formula = cDensity ~ 1, locations = idwHAD, newdata = grd)
    idwDEN.output = as.data.frame(idwDEN)     ## Output defined as dataframe with names
    names(idwDEN.output)[1:3] <- c("Easting", "Northing", "Density")
    
    coordinates(idwHAB) <- ~Easting + Northing
    if(lake != "Acton") proj4string(idwHAB) <- CRS("+proj=utm +zone=17") 
    if(lake == "Acton") proj4string(idwHAB) <- CRS("+proj=utm +zone=16")
    idwBIO <- idw(formula = lnBiomass ~ 1, locations = idwHAB, newdata = grd)
    idwBIO.output = as.data.frame(idwBIO)     ## Output defined as dataframe with names
    names(idwBIO.output)[1:3] <- c("Easting", "Northing", "lnBiomass")
    
    if(!exists("WQdataset")) {  ## creates dataset with interpolated data
        WQdataset <- idwt.output
        WQdataset$DO <- idwo.output$DO
        WQdataset$Density <- idwDEN.output$Density
        WQdataset$lnBiomass <- idwBIO.output$lnBiomass
        WQdataset$Depth <- rep(doi,nrow(WQdataset))
    } else {                    ## appends data to dataset once original created
        temp_dataset <- idwt.output
        temp_dataset$DO <- idwo.output$DO
        temp_dataset$Density <- idwDEN.output$Density
        temp_dataset$lnBiomass <- idwBIO.output$lnBiomass
        temp_dataset$Depth <- rep(doi,nrow(idwt.output))
        WQdataset <-rbind(WQdataset, temp_dataset)
        rm(temp_dataset)
    }
}

WQdataset <- WQdataset[c(-4)]     ## removes added var1.var column


#####  SHORE EXTRACTION AND PRINTING  #####
## Loops through number of shore polygons and creates map
for (i in 1:nlevels(shoreoutline$piece)) {
    if(lake == "Tappan" && i == 1) {
        poly <- subset(shoreoutline, shoreoutline$group == 7.1, select = c(long, lat))
    } else {
        poly <- subset(shoreoutline, shoreoutline$piece == i, select = c(long, lat))
    }
    polyP <- Polygon(cbind(poly$long, poly$lat)) ## Add lat & lon to object
    polySP <- SpatialPolygons(list(Polygons(list(polyP), ID=i)))
    outside <- gDifference(outside, polySP)
    nam <- paste("polySP", i, sep="")             ## Create sequential names
    assign(nam, polySP)                           ## Assign name
}

## Creates list of ggplot commands to plot all island polygons
islands <- list()
if(nlevels(shoreoutline$piece) >= 2) {
    for (ii in 2:nlevels(shoreoutline$piece)) {
        island <- paste("polySP", ii, sep="")
        island <- get(island)
        islands <- c(islands, geom_polygon(data=island, aes(x=long, y=lat), fill="darkgreen", asp = 1))
    }
}


## Finds only those Easting and Northing coordinates within lake boundary
lakepts <- idwt.output
coordinates(lakepts) <- ~Easting + Northing
if(lake != "Acton") proj4string(lakepts) <- CRS("+proj=utm +zone=17") 
if(lake == "Acton") proj4string(lakepts) <- CRS("+proj=utm +zone=16")
lakepts <- gDifference(lakepts, outside)
lakepts <- as.data.frame(lakepts)

## Goes through each coordinate set in lakepts and finds the depth
for(c in 1:nrow(lakepts)){
    easting <- lakepts[c,1]
    northing <- lakepts[c,2]
    
    ## Finds bottom depth for ploting  (used idwd.oupt, changed to depthd3D)
    dnearest <- sqrt((depthd3D[, 1] - easting)^2 + (depthd3D[, 2] - northing)^2)
    lakepts[c, 3] <- round(max(depthd3D[which(dnearest == min(dnearest)), 3]))
}

#####  BIOENERGETIC MODEL ######
soichar <- unique(bioparam$Species)
soi <- tk_select.list(as.character(soichar), title = "Pick a Species:")
{
  
    ## Read in parameters from bioenergetic database (user defined params at end)
    roi <- which(bioparam$Species == soi)
    CEQ = bioparam$CEQ[roi]
    CA = bioparam$CA[roi]
    CB = bioparam$CB[roi]
    CQ = bioparam$CQ[roi]
    CTO = bioparam$CTO[roi]
    CTM = bioparam$CTM[roi]
    CTL = bioparam$CTL[roi]
    CKone = bioparam$CK1[roi]
    CKfour = bioparam$CK4[roi]

    REQ = bioparam$REQ[roi]
    RA = bioparam$RA[roi]
    RB = bioparam$RB[roi]
    RQ = bioparam$RQ[roi]
    RTO = bioparam$RTO[roi]
    RTM = bioparam$RTM[roi]
    RTL = bioparam$RTL[roi]
    RKone = bioparam$RK1[roi]
    RKfour = bioparam$RK4[roi]
    ACT = bioparam$ACT[roi]
    BACT = bioparam$BACT[roi]
    SDA = bioparam$SDA[roi]

    EGEXEQ = bioparam$EGEXEQ[roi]
    FA = bioparam$FA[roi]
    FB = bioparam$FB[roi]
    FG = bioparam$FG[roi]
    UA = bioparam$UA[roi]
    UB = bioparam$UB[roi]
    UG = bioparam$UG[roi]

    predW = 1000
    if(soi == "LMB2") predW = 204 ##weight based on l/w relationship developed using LMB from Hoover/Acton (avg age-2 length = 243 mm)
    if(soi == "LMB3") predW = 404 ##weight based on l/w relationship developed using LMB from Hoover/Acton (avg age-3 length = 299 mm)
    if(soi == "SAE1") predW = 439 ##weight based on l/w relationship developed using SAE from Hoover/Acton (avg age-1 length = 358 mm)
    if(soi == "SAE2") predW = 1016 ##weight based on l/w relationship developed using SAE from Hoover/Acton (avg age-2 length = 463 mm)
    if(soi == "SAE3") predW = 1700 ##weight based on l/w relationship developed using SAE from Hoover/Acton (avg age-3 length = 542 mm)
    if(soi == "WAE2") predW = 469 ##weight based on l/w relationship developed using WAE from Ohio reservoirs statewide (avg age-2 length = 357 mm)
    if(soi == "WAE3") predW = 757 ##weight based on l/w relationship developed using WAE from Ohio reservoirs statewide (avg age-3 length = 414 mm)

    predED = bioparam$ED[roi]

    shadW = 14   ## Weight of a 90 mm shad, TL-W equation from Denlinger et al 2006
    ## Avg weight of prey below maximum gape size for each species
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-05-03") shadW = 25 
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-06-30") shadW = 30
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-08-09") shadW = 7
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-10-13") shadW = 12
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-04-27") shadW = 12
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-06-22") shadW = 23
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-08-10") shadW = 6
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-10-12") shadW = 7
   
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-06-30") shadW = 36
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-08-09") shadW = 10
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-10-13") shadW = 23
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-04-27") shadW = 17
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-06-22") shadW = 25
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-08-10") shadW = 8
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-10-12") shadW = 8
    
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-06-30") shadW = 30
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-08-09") shadW = 7
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-10-13") shadW = 12
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-04-27") shadW = 12
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-06-22") shadW = 23
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-08-10") shadW = 6
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-10-12") shadW = 7
    
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-06-30") shadW = 35
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-08-09") shadW = 9
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-10-13") shadW = 18
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-04-27") shadW = 13
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-06-22") shadW = 23
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-08-10") shadW = 7
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-10-12") shadW = 7
    
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-06-30") shadW = 36
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-08-09") shadW = 10
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-10-13") shadW = 23
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-04-27") shadW = 16
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-06-22") shadW = 24
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-08-10") shadW = 8
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-10-12") shadW = 8
    
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-06-30") shadW = 30
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-08-09") shadW = 7
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-10-13") shadW = 12
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-04-27") shadW = 12
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-06-22") shadW = 23
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-08-10") shadW = 5
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-10-12") shadW = 7    
    
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-05-03") shadW = 25
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-06-30") shadW = 34
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-08-09") shadW = 8
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-10-13") shadW = 13
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-04-27") shadW = 12
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-06-22") shadW = 23
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-08-10") shadW = 7
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-10-12") shadW = 7

    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 14
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 16
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 19
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 19
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 13
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 18
    
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 17
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 21
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 21
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 22
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 21
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 22   
    
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 14
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 16
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 19
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 19
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 13
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 18       
    
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 15
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 18
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 20
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 21
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 20
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 22        
    
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 16
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 21
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 21
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 22
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 21
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 22    
    
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 14
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 16
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 19
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 18
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 12
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 18
      
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadW = 14
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadW = 17
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadW = 20
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadW = 21
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadW = 18
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadW = 21    
    
      
    shadED = 5063 ##Composite value from various studies 
    ## Weighted avg ED of prey below maximum gape size for each species
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-05-03") shadED = 5103 
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-08-09") shadED = 5099
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2004-10-13") shadED = 5019
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-04-27") shadED = 5073
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-06-22") shadED = 5087
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-08-10") shadED = 5087
    if(lake == "Acton" & soi == "LMB2" & dayoi == "2005-10-12") shadED = 5099
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-08-09") shadED = 5102
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2004-10-13") shadED = 5012
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-04-27") shadED = 4622
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-06-22") shadED = 5077
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-08-10") shadED = 5060
    if(lake == "Acton" & soi == "LMB3" & dayoi == "2005-10-12") shadED = 5077
    
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-08-09") shadED = 5099
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2004-10-13") shadED = 5019
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-04-27") shadED = 5073
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-06-22") shadED = 5087
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-08-10") shadED = 5087
    if(lake == "Acton" & soi == "SAE1" & dayoi == "2005-10-12") shadED = 5099
    
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-08-09") shadED = 5100
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2004-10-13") shadED = 5006
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-04-27") shadED = 5022
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-06-22") shadED = 5081
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-08-10") shadED = 5076
    if(lake == "Acton" & soi == "SAE2" & dayoi == "2005-10-12") shadED = 5091
    
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-08-09") shadED = 5102
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2004-10-13") shadED = 5017
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-04-27") shadED = 5031
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-06-22") shadED = 5079
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-08-10") shadED = 5063
    if(lake == "Acton" & soi == "SAE3" & dayoi == "2005-10-12") shadED = 5079
    
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-08-09") shadED = 5100
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2004-10-13") shadED = 4988
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-04-27") shadED = 5073
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-06-22") shadED = 5094
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-08-10") shadED = 5086
    if(lake == "Acton" & soi == "WAE2" & dayoi == "2005-10-12") shadED = 5099 
    
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-05-03") shadED = 5103
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-06-30") shadED = 5103
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-08-09") shadED = 5099
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2004-10-13") shadED = 4984
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-04-27") shadED = 5073
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-06-22") shadED = 5082
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-08-10") shadED = 5083
    if(lake == "Acton" & soi == "WAE3" & dayoi == "2005-10-12") shadED = 5095
  
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5063
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 5069
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5083
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5064
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5047
    if(lake == "Hoover" & soi == "LMB2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5103
      
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5024
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 4987
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5070
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5030
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5057
    if(lake == "Hoover" & soi == "LMB3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5095   
      
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5063
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 5069
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5083
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5064
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5047
    if(lake == "Hoover" & soi == "SAE1" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5103       
      
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5043
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 5032
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5084
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5059
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5066
    if(lake == "Hoover" & soi == "SAE2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5103        
      
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5022
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 4983
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5075
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5044
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5056
    if(lake == "Hoover" & soi == "SAE3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5096   
      
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5062
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 5067
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5082
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5061
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5043
    if(lake == "Hoover" & soi == "WAE2" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5103
      
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-05-08" | dayoi == "2012-04-12" | dayoi == "2012-04-2012" | dayoi == "2012-05-08" 
      | dayoi == " 2012-05-22") shadED = 5054
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-06-05" | dayoi == "2012-06-19") shadED = 5072
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-07-13" | dayoi == "2012-07-03") shadED = 5084
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-08-07" | dayoi == "2011-08-25" | dayoi == "2012-08-03" | dayoi == "2012-08-14" 
      | dayoi == " 2012-08-29") shadED = 5074
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-09-05" | dayoi == "2012-09-11" | dayoi == "2012-09-25") shadED = 5075
    if(lake == "Hoover" & soi == "WAE3" & dayoi == "2006-10-10" | dayoi == "2011-10-25" | dayoi == "2012-10-10" | dayoi == "2012-10-23") shadED = 5103

    GRPd <- data.frame()
    rownum = 0
    
    ## Generate 3D GRP data for each depth

    for(ccc in 1:nrow(lakepts)) {           ## loop through every cell in 3D array
        easting <- lakepts[ccc,1]
        northing <- lakepts[ccc,2]
        nearest <- sqrt((WQdataset$Easting - easting)^2 + (WQdataset$Northing - northing)^2)
        nearestdata <- WQdataset[which(nearest == min(nearest)), ]
    
        neard <- subset(lakepts, lakepts$x == easting & lakepts$y == northing)
        depthmax <- max(neard[, 3])
    
        for(z in 0:depthmax) {
            rownum = rownum +1
        
            roi <- nearestdata[which(nearestdata$Depth == z), ]
            if(z > maxz)  {
                dc = z - maxz
                roi <- nearestdata[which(nearestdata$Depth == z-dc), ]
            }
        
            wt = roi$Temp
            oxy = roi$DO
            preyd = roi$Density
            preyb = exp(roi$lnBiomass) - 1
            
            if(soi == "SMB") {         ## extra for SMB from Whitledge et al 2003
                if(wt > 22) CQ = 1.95
            }
            
            ##  Consumption Equations
            if(CEQ == 2) {     ## Consumption equation 2
                Y = log(CQ) * (CTM - CTO + 2)
                Z = log(CQ) * (CTM - CTO)
                X = (Z^2 * (1 + (1 + 40 / Y)^0.5)^2) / 400
                V = (CTM - wt) / (CTM - CTO)
                
                ftc = V^X * exp(X * (1 - V))
                if(V^X == "NaN") ftc = 0.00001  ## If water temp>CTM, use low ftc
                
                Cmax = CA * predW^CB        ## this is g prey/g pred
                
                ##  Species-specific oxygen and temperature functions
                if(soi == "LMB2") {  
                    fdo = 1 / (1 + exp(-1 * (oxy - 2.281586) / 0.345078)) ## Type III Stewart et al 1967
                }
                if(soi == "LMB3") {  
                  fdo = 1 / (1 + exp(-1 * (oxy - 2.281586) / 0.345078)) ## Type III Stewart et al 1967
                }
                if(soi == "MUS") {
                    fdo = 1 / (1 + exp(-1 * (oxy - 2.563859) / 0.37121))  ## Type III Adelman and Smith 1970
                }
                if(soi == "WHC") {
                    fdo = 1 / (1 + exp(-1 * (oxy - 3.1485) / 0.345078))  ## Type III
                }
                if(soi == "SAE1") {          
                    fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
                }
                if(soi == "SAE2") {          
                  fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
                }
                if(soi == "SAE3") {          
                  fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
                }
                if(soi == "SMB") {
                    fdo = 1 / (1 + exp(-1 * (oxy - 3.687002) / 0.22578))  ## Type III, Bulkley 1975
                }
                if(soi == "WAE2") {          
                    fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
                }
                if(soi == "WAE3") {          
                  fdo = 1 / (1 + exp(-1 * (oxy - 3) / 0.455))  ## Type III, Brandt et al. 2011
                }
                if(soi == "YEP") {     
                    fdo = 1 / (1 + exp(-1 * (oxy - 2.370402) / 0.330795))  ## Type III, Roberts data
                }
                if(soi == "BCF") {     ## Torrans 2012 & 2008, Green et al 2012
                    fdo = 1 / (1 + exp(-1 * (oxy - 1.19497) / 0.17284)) ## Type III
                }
                
                C = Cmax * ftc * fdo * preyb / (0.865 + preyb)  ##Constantini            
            }
        
            if(CEQ == 3) {     ## Consumption equation 3
                Gone = (1 / (CTO - CQ)) * log((0.98 * (1 - CKone)) / (CKone * 0.02)) 
                Lone = exp(Gone * (wt - CQ))
                Ka = (CKone * Lone) / (1 + CKone * (Lone - 1))
                
                Gtwo = (1 / (CTL - CTM)) * log((0.98 * (1 - CKfour)) / (CKfour * 0.02))
                Ltwo = exp(Gtwo * (CTL - wt))
                Kb = (CKfour * Ltwo) / (1 + CKfour * (Ltwo - 1))
                
                ftc = Ka * Kb
                
                Cmax = CA * predW^CB     ## this is g prey/g pred
                
                if(soi == "STB") {  ## Brandt et al 2009 and Chiba 1998
                    fdo = 1 / (1 + exp(-1 * (oxy - 2.190932) / 0.290913)) ## Type III        
                }
                if(soi == "HSB") {  ## Emily's data
                    fdo = 1 / (1 + exp(-1 * (oxy - 2.190932) / 0.290913)) ## Type III
                }
                
                C = Cmax * ftc * fdo * preyb / (0.865 + preyb)
            } 
            
            if(CEQ == 4) {     ## Consumption equation 4
              ftc = (CQ*wt + CK1*wt^2 + CK4*wt^3)
              
              Cmax = CA * predW^CB     ## this is g prey/g pred
              
              if(soi == "MUS") {  
                fdo = 1 / (1 + exp(-1 * (oxy - 2.563859) / 0.37121))      
              }
              
              C = Cmax * ftc * fdo * preyb / (0.865 + preyb)
            } 
            
            ##  Respiration equations
            if(REQ == 1) {     ## Respiration equation 1
                ftr = exp(RQ * wt)
                VEL = RKone * predW^RKfour  ## no need to use equation set b/c T>RTL
                activity  = exp(RTO * VEL)
                S = SDA * (C - (FA * C))
                R = (RA * predW^RB * ftr) * activity * (13560/shadED)+S #Respiration equation from Zhang et al. 2014
            }
            
            if(REQ == 2) {     ## Respiration equation 2
                V = (RTM - wt) / (RTM - RTO)
                Z = log(RQ) * (RTM - RTO)
                Y = log(RQ) * (RTM - RTO +2)
                X = (Z^2 * (1 + (1 + 40 / Y)^0.5)^2) / 400
                
                ftr = V^X * exp(X * (1 - V))
                activity  = ACT
                
                S = SDA * (C - (FA * C))
                R = (RA * predW^RB * ftr) * activity * (13560/shadED)+S #Respiration equation from Zhang et al. 2014
            }
            
            ##  Egestion and excretion equations
            if(EGEXEQ == 1) {     ## Egestion and excretion equation 1
                FE = FA * C
                U = UA * (C - FE)
            }
            
        
            if(EGEXEQ == 2) {     ## Egestion and excretion equation 2
                p = C/Cmax
                FE = FA * wt^FB * exp(FG * p) * C
                U = UA * wt^UB * exp(UG * p) * (C - FE)
            }
            
             ## Final Equation
            GRP = (shadED / predED)*(C-(R+FE+U))
        
            ## Stores data in master output dataset
            GRPd[rownum,1] <- easting
            GRPd[rownum,2] <- northing
            GRPd[rownum,3] <- z
            GRPd[rownum,4] <- wt
            GRPd[rownum,5] <- oxy
            GRPd[rownum,6] <- preyd
            GRPd[rownum,7] <- preyb
            GRPd[rownum,8] <- GRP
        
        }
    }    

    colnames(GRPd) <- c("Easting", "Northing", "Depth", "Temp", "DO",
                        "Density", "Biomass", "GRP")

    maxGRP = max(GRPd[ , 8:8])
    medGRP = median(GRPd[ , 8])
    minGRP = min(GRPd[ , 8:8])


    ## Generate summary data and write to file
    Summary3D <- data.frame()
    Summary3D[1,1] <- lake
    Summary3D[1,2] <- as.character(dayoi)
    Summary3D[1,3] <- soi
    Summary3D[1,4] <- "Prop hypoxic"
    DOoi <- which(GRPd$DO < 2.0)
    Summary3D[1,5] <- length(DOoi)/nrow(GRPd)
    

    Summary3D[2,2] <- "Mean GRP"
    Summary3D[2,3] <- "Max GRP"
    Summary3D[2,4] <- "Mean GRP pos cells"
    Summary3D[2,5] <- "Proportion pos cells"
    Summary3D[2,6] <- "Mean Temperature"
    Summary3D[2,7] <- "Mean Biomass"
    
    Summary3D[3,1] <- "With hypoxia"
    summaryGRP <- na.omit(GRPd)
    Summary3D[3,2] <- mean(GRPd$GRP)
    Summary3D[3,3] <- max(GRPd$GRP)
    grpoi <- which(GRPd$GRP > 0)
    Summary3D[3,4] <- mean(GRPd$GRP[grpoi])
    Summary3D[3,5] <- length(which(GRPd$GRP > 0)) / nrow(GRPd)
    Summary3D[3,6] <- mean(GRPd$Temp)
    Summary3D[3,7] <- mean(GRPd$Biomass)    

    ## Print results
    write.table(GRPd, paste("C:/Users/richa/OneDrive/Desktop/", dayoi, "_", soi, "_", lake, ".txt", sep = ""), row.names = TRUE)
    write.table(Summary3D, paste("C:/Users/richa/OneDrive/Desktop/", dayoi, "_", soi, "_", lake, "_sum.txt", sep = ""), row.names = TRUE)

}
