mem <- function(expdata, markers=NULL, cluster.id, ref.pop=NULL, IQR.thresh=NULL, statistic="median")
{
  #' @author Julian Spagnuolo
  #' @references Diggins, K.E., et al, (2017) Nature Methods http://dx.doi.org/10.1038/nmeth.4149
  
  if(is.null(markers))
  {
    markers <- names(expdata)
  }
  data <- expdata[,markers]
  data$cluster.id <- cluster.id
  
  #  Initialize variables
  #markers = as.vector(markers)
  n.markers = length(markers)
  n.cells = nrow(expdata)
  n.pops = length(unique(data$cluster.id))
  pop.names = unique(data$cluster.id)
  
  MAGpop = matrix(nrow=n.pops,ncol=n.markers)
  MAGref = matrix(nrow=n.pops,ncol=n.markers)
  IQRpop = matrix(nrow=n.pops,ncol=n.markers)
  IQRref = matrix(nrow=n.pops,ncol=n.markers)
  
  if(statistic == "median")
  {
    # Get population medians and IQRs
    for(i in 1:n.pops)
    {
      pop = pop.names[i]
      MAGpop[i,] = abs(apply(data[which(data$cluster.id == pop),markers],2,FUN=median,na.rm=TRUE))
      IQRpop[i,] = apply(data[which(data$cluster.id == pop),markers],2,FUN=IQR,na.rm=TRUE)
    }
    
    # Get reference population medians and IQRs
    if(!is.null(ref.pop))
    {
      if(length(ref.pop) > 1)
      {
        ref.pop.data <- data[which(data$cluster.id == ref.pop), markers]
        
        MAGref = matrix(rep(abs(apply(ref.pop.data,2,FUN=median,na.rm=TRUE)), each = n.pops),n.pops)
        IQRref = matrix(rep(apply(ref.pop.data,2,FUN=IQR,na.rm=TRUE),each = n.pops),n.pops)
      }
        else
        {
          MAGref = matrix(rep(abs(apply(data[which(data$cluster.id == ref.pop), markers],2,FUN=median,na.rm=TRUE)), each = n.pops),n.pops)
          IQRref = matrix(rep(apply(data[which(data$cluster.id == ref.pop), markers],2,FUN=IQR,na.rm=TRUE),each = n.pops),n.pops)
        }
    }
    if(is.null(ref.pop))
    {
      for(i in 1:n.pops)
      {
        pop = pop.names[i]
        MAGref[i,] = abs(apply(data[which(data$cluster.id != pop), markers],2,FUN=median,na.rm=TRUE))
        IQRref[i,] = apply(data[which(data$cluster.id != pop), markers],2,FUN=IQR,na.rm=TRUE)
      }
    }
  }
  
  if(statistic == "mode")
  {
    #IQRpop <- matrix(nrow = length(unique(data$cluster.id)), ncol = length(markers), dimnames = list(unique(data$cluster.id), markers))
    #IQRref <- matrix(nrow = length(unique(data$cluster.id)), ncol = length(markers), dimnames = list(unique(data$cluster.id), markers))
    for(i in 1:length(unique(data$cluster.id)))
    {
      for(n in 1:length(markers))
      {
        MAGpop[i,n] <- dmode(x=data[which(data$cluster.id == pop.names[i]), markers[n]], gridsize = 14000)
      }
      IQRpop[i,] = apply(data[which(data$cluster.id == pop.names[i]),markers],2,FUN=IQR,na.rm=TRUE)
    }
    
    # Get reference population modes and IQRs
    if(!is.null(ref.pop))
    {
      MAGref = vector(length=length(markers))
      names(MAGref) <- markers
      for(n in markers)
      {
        MAGref[n] <- dmode(data[which(data$cluster.id == ref.pop), n], gridsize = 14000)
      }
      MAGref = matrix(data=MAGref, nrow=n.pops, ncol = length(markers), byrow=TRUE)
      IQRref = matrix(rep(apply(data[which(data$cluster.id == ref.pop), markers],2,FUN=IQR,na.rm=TRUE),each = n.pops),n.pops)
    }
    if(is.null(ref.pop))
    {
      MAGref <- matrix(nrow=n.pops, ncol=length(markers), dimnames=list(pop.names, markers))
      for(i in 1:length(n.pops))
      {
        for(n in markers)
        {
          MAGref[i,n] <- dmode(data[which(data$cluster.id != pop.names[i]),n], gridsize = 14000)
        }
        IQRref[i,] = apply(data[which(data$cluster.id != pop.names[i]), markers],2,FUN=IQR,na.rm=TRUE)
      }
    }
  }
  
  # Set and apply IQR threshold
  if(is.null(IQR.thresh))
  {
    IQR.thresh = 0.5
  }
  
  if(IQR.thresh=="auto")
  {
    #   Use universal IQR threshodling
    IQR.thresh.pop = matrix()
    IQR.thresh.ref = matrix()
    for(i in 1:length(markers))
    {
      MAGpop.belowThresh <- MAGpop[,i] <= quantile(MAGpop[,i])[2]
      IQR.thresh.pop[i] <- min(IQRpop[,i][MAGpop.belowThresh])
      MAGref.belowThresh <- MAGref[,i]<=quantile(MAGref[,i])[2]
      IQR.thresh.ref[i] <- min(IQRref[,i][MAGref.belowThresh])
    }
    IQR.thresh = mean(c(IQR.thresh.pop,IQR.thresh.ref))
    print(IQR.thresh)
  }
  
  for(i in 1:length(markers))
  {
    IQRpop[,i] = pmax(IQRpop[,i],IQR.thresh)
    IQRref[,i] = pmax(IQRref[,i],IQR.thresh)
  }
  
  if(n.pops < 4)
  {
    IQR.thresh = 0.5
  }
  
  # Calculate MEM scores
  MAG.diff = MAGpop - MAGref
  MEM.matrix = abs(MAGpop - MAGref) + (IQRref/IQRpop) - 1
  MEM.matrix[!(MAG.diff >= 0)] <- (-MEM.matrix[!(MAG.diff >= 0)])
  
  # Put MEM values on -10 to +10 scale
  scale.max = max(abs(MEM.matrix[,c(1:ncol(MEM.matrix)-1)]))
  MEM.matrix = cbind((MEM.matrix[,c(1:ncol(MEM.matrix)-1)]/scale.max)*10, MEM.matrix[,ncol(MEM.matrix)])
  
  #Rename rows and columns of all matrices
  rename.table <- function(x)
  {
    colnames(x) = markers[1:length(markers)]
    rownames(x) = pop.names
    return(x)
  }
  
  # Apply rename_table function across matrices
  object.list.labeled <- lapply(list(MAGpop[,1:length(markers)],MAGref[,1:length(markers)],IQRpop[,1:length(markers)],IQRref[,1:length(markers)],MEM.matrix[,1:length(markers)]),rename.table)
  # List all matrices for export
  all.values <- list("MAGpop" = object.list.labeled[1],"MAGref"=object.list.labeled[2],"IQRpop"=object.list.labeled[3],"IQRref"=object.list.labeled[4],"MEM.matrix"=object.list.labeled[5])
  return(all.values)
}