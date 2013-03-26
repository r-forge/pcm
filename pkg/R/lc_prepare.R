lc_prepare<-function(lc.africa.file=lc.africa.file,lc.europe.file=lc.europe.file,out.dir=out.dir,
                     out.proj=out.proj,out.extent=out.extent,out.dim=out.dim,Africa=Africa){
  #sanity check
  if (missing('lc.europe.file')){print("ERROR:lc.europe.file variable not specified. EXIT!");return(1)}
  if (missing('out.dir')){print("ERROR:out.dir variable not specified. EXIT!");return(1)}
  if (missing('out.proj')){print("ERROR:out.proj variable not specified. EXIT!");return(1)}
  if (missing('out.extent')){print("ERROR:out.extent variable not specified. EXIT!");return(1)}
  if (missing('out.dim')){print("ERROR:out.dim variable not specified. EXIT!");return(1)}
  if (missing('Africa')){Africa=T}
  
  if (mode(lc.europe.file)!='character'){print("ERROR:lc.europe.file is not a character string. EXIT!");return(1)}
  if (mode(out.dir)!='character'){print("ERROR:out.dir is not a character string. EXIT!");return(1)}  
  if (mode(out.proj)!='character'){print("ERROR:out.proj is not a character string. EXIT!");return(1)} 
  if (mode(Africa)!='logical'){print("ERROR:keyword Africa accepts only TRUE ot FALSE. EXIT!");return(1)}
  
  out.extent<-as.numeric(out.extent)
  if (length(which(is.na(out.extent)))!=0 | length(out.extent) !=4){print("ERROR:out.extent is not a 4 elements numeric vector in a form c(ULX,ULY,LRX,LRY). EXIT!");return(1)} 
  out.dim<-as.numeric(out.dim)
  if (length(which(is.na(out.dim)))!=0 | length(out.dim) !=2){print("ERROR:out.dim is not a 2 elements numeric vector in a form c(ncol,nrows). EXIT!");return(1)}
  if (! file.exists(lc.europe.file)){print(paste("ERROR:lc.europe.file does not exists:",lc.europe.file));return(1)}
  if (! file.exists(out.dir)){print(paste("ERROR:out.dir does not exists:",out.dir));return(1)}
  
#start processing
  
  out.proj<-try(CRS(out.proj))
  if (class(out.proj)=='try-error'){print("ERROR:Specified proj4 arguments for the out.proj are not valid. EXIT!");return(1)}
  
  if (out.extent[1]>=out.extent[3]){
    print(paste('ERROR:ULX=',out.extent[1],'should be lower than LRX=',out.extent[3],'.EXIT!'))
    return(1)
  }
  if (out.extent[2]<=out.extent[4]){
    print(paste('ERROR:ULY=',out.extent[2],'should be greater than LRY=',out.extent[4],'.EXIT!'))
    return(1)
  }    
  out.extent<-matrix(c(out.extent[1],out.extent[4],out.extent[3],out.extent[2]),nrow=2,dimnames=list(c('x','y'),c('min','max'))) 
  print('out.extent:')
  print(out.extent) 
  pts<-expand.grid(seq(out.extent['y','max'],out.extent['y','min'],length=out.dim[2]),
                   seq(out.extent['x','min'],out.extent['x','max'],length=out.dim[1]))
  pts<-data.frame(x=pts[,2],y=pts[,1])
  coordinates(pts)=~x+y
  proj4string(pts)=out.proj
  pts = spTransform(pts,CRS('+proj=longlat +datum=WGS84 +no_defs'))
  gc()
  
  #Europe
  lc.europe<-try(readGDAL(lc.europe.file,half.cell=c(0,1)))
  if (class(lc.europe)=='try-error'){print("ERROR:GDAL could not open lc.europe.file EXIT!");return(1)}
  
  lc.europe.extent<-lc.europe@bbox
  lc.europe.extent['x','max']<-lc.europe.extent['x','max']-lc.europe@grid@cellsize[1]
  lc.europe.extent['y','min']<-lc.europe.extent['y','min']+lc.europe@grid@cellsize[2]
  print('lc.europe.extent:')
  print(lc.europe.extent)
  grid_europe<-expand.grid(seq(from=lc.europe.extent['y','max'],to=lc.europe.extent['y','min'],length=lc.europe@grid@cells.dim[2]),
                           seq(from=lc.europe.extent['x','min'],to=lc.europe.extent['x','max'],length=lc.europe@grid@cells.dim[1]))
  idx<-which(grid_europe[,1]>=min(pts@coords[,2]) & grid_europe[,1]<=max(pts@coords[,2]) & 
       grid_europe[,2]>=min(pts@coords[,1]) & grid_europe[,2]<=max(pts@coords[,1]))
  gc()
  if (length(idx) != 0){
    neighbor<-get.knnx(grid_europe[idx,],cbind(pts@coords[,2],pts@coords[,1]),k=1,algorithm="kd_tree")
    gc()
    idx_ok <- which(neighbor$nn.dist < mean(lc.europe@grid@cellsize)*1.5)
    neighbor <- neighbor$nn.index[idx_ok]
    data_out<-rep(20,out.dim[1]*out.dim[2])
    data_out[idx_ok]<-c(matrix(lc.europe@data$band1,nrow=lc.europe@grid@cells.dim[2],byrow=T))[idx][neighbor]
    data_out[data_out==16]<-22
    remove(lc.europe,grid_europe,neighbor,idx_ok,idx)
  }else{
    data_out<-rep(255,out.dim[1]*out.dim[2])
    remove(lc.europe,grid_europe,idx)
  }
  gc()
  #Africa
  if (Africa){
      if (missing('lc.africa.file')){print("ERROR:lc.africa.file variable not specified. EXIT!");return(1)}
      if (mode(lc.africa.file)!='character'){print("ERROR:lc.africa.file is not a character string. EXIT!");return(1)}
      if (! file.exists(lc.africa.file)){print(paste("ERROR:lc.africa.file does not exists:",lc.africa.file));return(1)}
      lc.africa<-try(readGDAL(lc.africa.file,half.cell=c(0,1)))
      if (class(lc.africa)=='try-error'){print("ERROR:GDAL could not open lc.africa.file EXIT!");return(1)}
    
      lc.africa.extent<-lc.africa@bbox
      lc.africa.extent['x','max']<-lc.africa.extent['x','max']-lc.africa@grid@cellsize[1]
      lc.africa.extent['y','min']<-lc.africa.extent['y','min']+lc.africa@grid@cellsize[2]
      lc.africa.cellsize<-lc.africa@grid@cellsize
      lc.africa.cells.dim<-lc.africa@grid@cells.dim
      print('lc.africa.extent:')
      print(lc.africa.extent)
      grid_africa<-expand.grid(seq(from=lc.africa.extent['y','max'],to=lc.africa.extent['y','min'],length=lc.africa@grid@cells.dim[2]),
                               seq(from=lc.africa.extent['x','min'],to=lc.africa.extent['x','max'],length=lc.africa@grid@cells.dim[1]))
      idx<-which(grid_africa[,1]>=min(pts@coords[,2]) & grid_africa[,1]<=max(pts@coords[,2]) & 
                   grid_africa[,2]>=min(pts@coords[,1]) & grid_africa[,2]<=max(pts@coords[,1]))
      gc()
      neighbor<-get.knnx(grid_africa[idx,],cbind(pts@coords[,2],pts@coords[,1]),k=1,algorithm="kd_tree")
      gc()
      data_africa<-as.integer(matrix(lc.africa@data$band1,nrow=lc.africa@grid@cells.dim[2],byrow=T))
      remove(lc.africa)
      gc()
      data<-data_africa
      data_africa[data==22]<-25
      data_africa[data==23]<-26
      data_africa[data==24]<-27
      data_africa[data==16]<-24
      data_africa[data==27]<-22
      data_africa[data==25]<-27
      data_africa[data==26]<-24
      idx_ok <- which(neighbor$nn.dist < mean(lc.africa.cellsize)*1.5 & data_africa[idx][neighbor$nn.index] != 0)
      neighbor <- neighbor$nn.index[idx_ok]
      data_out[idx_ok]<-data_africa[idx][neighbor] 
      remove(data,data_africa,grid_africa,neighbor,idx,idx_ok)
      gc()  
  }
  idx_water_image<-matrix(as.numeric(data_out %in% c(20,255)),nrow=out.dim[2])
  shoreline<-as.numeric(imgConvolve(imagedata(idx_water_image,type='grey',ncol=out.dim[2],nrow=out.dim[1]),shoreline_kernel(17),0) ==0)
  
  x.cell.size<-(out.extent['x','max']-out.extent['x','min'])/(out.dim[1]-1)
  y.cell.size<-(out.extent['y','max']-out.extent['y','min'])/(out.dim[2]-1)
  
  check <- try(writeRaster(raster(matrix(data_out*shoreline,ncol=out.dim[1]), xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                              ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                       filename=file.path(out.dir,'lc.tif'), format='GTiff',datatype='INT2U', overwrite=T,
                       options=c("COMPRESS=LZW","PREDICTOR=2")))
  if (class(check)=='try-error'){print("ERROR:Writing file to disk failed.EXIT!");return(1)}
  return(0)
}
