dem_prepare<-function(NWP.dem.file=NWP.dem.file,dem.file=dem.file,out.dir=out.dir,out.proj=out.proj,
                      out.extent=out.extent,out.dim=out.dim,NWP.dem.proj=NWP.dem.proj,
                      NWP.dem.extent=NWP.dem.extent,dem.proj=dem.proj,dem.extent=dem.extent,
                      NWP.dem.half.cell=NWP.dem.half.cell,dem.half.cell=dem.half.cell,
                      NWP.dem.scale=NWP.dem.scale,NWP.dem.offset=NWP.dem.offset,dem.scale=dem.scale,
                      dem.offset=dem.offset,NWP.dem.band=NWP.dem.band){
  #sanity check
  if (missing('NWP.dem.file')){print("ERROR:NWP.dem.file variable not specified. EXIT!");return(1)}
  if (missing('dem.file')){print("ERROR:dem.file variable not specified. EXIT!");return(1)}
  if (missing('out.dir')){print("ERROR:out.dir variable not specified. EXIT!");return(1)}
  if (missing('out.proj')){print("ERROR:out.proj variable not specified. EXIT!");return(1)}
  if (missing('out.extent')){print("ERROR:out.extent variable not specified. EXIT!");return(1)}
  if (missing('out.dim')){print("ERROR:out.dim variable not specified. EXIT!");return(1)}
  
  if (mode(NWP.dem.file)!='character'){print("ERROR:NWP.dem.file is not a character string. EXIT!");return(1)}
  if (mode(dem.file)!='character'){print("ERROR:dem.file is not a character string. EXIT!");return(1)}
  if (mode(out.dir)!='character'){print("ERROR:out.dir is not a character string. EXIT!");return(1)}  
  if (mode(out.proj)!='character'){print("ERROR:out.proj is not a character string. EXIT!");return(1)} 
  
  out.extent<-as.numeric(out.extent)
  if (length(which(is.na(out.extent)))!=0 | length(out.extent) !=4){print("ERROR:out.extent is not a 4 elements numeric vector in a form c(ULX,ULY,LRX,LRY). EXIT!");return(1)} 
  out.dim<-as.numeric(out.dim)
  if (length(which(is.na(out.dim)))!=0 | length(out.dim) !=2){print("ERROR:out.dim is not a 2 elements numeric vector in a form c(ncol,nrows). EXIT!");return(1)}
  if (! file.exists(NWP.dem.file)){print(paste("ERROR:NWP.dem.file does not exists:",NWP.dem.file));return(1)}
  if (! file.exists(dem.file)){print(paste("ERROR:dem.file does not exists:",dem.file));return(1)}
  if (! file.exists(out.dir)){print(paste("ERROR:out.dir does not exists:",out.dir));return(1)}
  
  if (missing('NWP.dem.half.cell')){NWP.dem.half.cell<-c(0,1)}
  if (missing('dem.half.cell')){dem.half.cell<-c(0,1)}
  if (missing('NWP.dem.band')){NWP.dem.band<-1}
  
  NWP.dem<-try(readGDAL(NWP.dem.file,half.cell=NWP.dem.half.cell,band=NWP.dem.band))
  if (class(NWP.dem)=='try-error'){print("ERROR:GDAL could not open NWP.dem.file. EXIT!");return(1)}
  dem<-try(readGDAL(dem.file,half.cell=dem.half.cell))
  if (class(dem)=='try-error'){print("ERROR:GDAL could not open dem.file. EXIT!");return(1)}
  
  if (missing(NWP.dem.proj)){
    NWP.dem.proj<-NWP.dem@proj4string@projargs
  }
  NWP.dem.proj<-try(CRS(NWP.dem.proj))
  if (class(NWP.dem.proj)=='try-error'){print("ERROR:Specified proj4 arguments for the NWP.dem.proj are not valid. EXIT!");return(1)}
  
  if (missing(NWP.dem.extent)){
    NWP.dem.extent<-NWP.dem@bbox
    NWP.dem.extent['x','max']<-NWP.dem.extent['x','max']-NWP.dem@grid@cellsize[1]
    NWP.dem.extent['y','min']<-NWP.dem.extent['y','min']+NWP.dem@grid@cellsize[2]
  }else{
    if (NWP.dem.extent[1]>=NWP.dem.extent[3]){
      print(paste('ERROR:ULX=',NWP.dem.extent[1],'should be lower than LRX=',NWP.dem.extent[3],'.EXIT!'))
      return(1)
    }
    if (NWP.dem.extent[2]<=NWP.dem.extent[4]){
      print(paste('ERROR:ULY=',NWP.dem.extent[2],'should be greater than LRY=',NWP.dem.extent[4],'.EXIT!'))
      return(1)
    }    
    NWP.dem.extent<-matrix(c(NWP.dem.extent[1],NWP.dem.extent[4],NWP.dem.extent[3],NWP.dem.extent[2]),nrow=2,dimnames=list(c('x','y'),c('min','max')))
  }
  print('NWP.dem.extent:')
  print(NWP.dem.extent)
  
  if (missing(dem.proj)){
    dem.proj<-dem@proj4string@projargs
  }
  dem.proj<-try(CRS(dem.proj))
  if (class(dem.proj)=='try-error'){print("ERROR:Specified proj4 arguments for the dem.file are not valid. EXIT!");return(1)}
  
  if (missing(dem.extent)){
    dem.extent<-dem@bbox
    dem.extent['x','max']<-dem.extent['x','max']-dem@grid@cellsize[1]
    dem.extent['y','min']<-dem.extent['y','min']+dem@grid@cellsize[2]
  }else{
    if (dem.extent[1]>=dem.extent[3]){
      print(paste('ERROR:ULX=',dem.extent[1],'should be lower than LRX=',dem.extent[3],'.EXIT!'))
      return(1)
    }
    if (dem.extent[2]<=dem.extent[4]){
      print(paste('ERROR:ULY=',dem.extent[2],'should be greater than LRY=',dem.extent[4],'.EXIT!'))
      return(1)
    }    
    dem.extent<-matrix(c(dem.extent[1],dem.extent[4],dem.extent[3],dem.extent[2]),nrow=2,dimnames=list(c('x','y'),c('min','max')))
  }
  print('dem.extent:')
  print(dem.extent)  
  
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
  
  if(missing('NWP.dem.scale')){NWP.dem.scale=1}
  if(missing('NWP.dem.offset')){NWP.dem.offset=0}
  if(missing('dem.scale')){dem.scale=1}
  if(missing('dem.offset')){dem.offset=0}
#start processing  
  NWP.dem<-list(
    y=seq(NWP.dem.extent['x','min'],NWP.dem.extent['x','max'],length=NWP.dem@grid@cells.dim[1]),
    x=seq(NWP.dem.extent['y','max'],NWP.dem.extent['y','min'],length=NWP.dem@grid@cells.dim[2]),
    z=matrix(NWP.dem@data$band1,ncol=NWP.dem@grid@cells.dim[1],byrow=T)*NWP.dem.scale+NWP.dem.offset
  )
  
  pts<-expand.grid(seq(out.extent['y','max'],out.extent['y','min'],length=out.dim[2]),
                   seq(out.extent['x','min'],out.extent['x','max'],length=out.dim[1]))
  pts<-data.frame(x=pts[,2],y=pts[,1])
  coordinates(pts)=~x+y
  proj4string(pts)=out.proj
  pts = spTransform(pts,NWP.dem.proj)
  idx<-pts@coords[,1]>=NWP.dem.extent['x','min'] & pts@coords[,1]<=NWP.dem.extent['x','max'] & 
       pts@coords[,2]>=NWP.dem.extent['y','min'] & pts@coords[,2]<=NWP.dem.extent['y','max']

  if (length(which(idx))==0){
    print('NWP.dem.extent and out.extent do not overlap.EXIT!')
    return(1)
  }
  
  NWP.dem.out<-rep(-32768,out.dim[1]*out.dim[2])
  NWP.dem.out[idx]<-interp.surface(NWP.dem,cbind(pts@coords[idx,2],pts@coords[idx,1]))
  gc()
  x.cell.size<-(out.extent['x','max']-out.extent['x','min'])/(out.dim[1]-1)
  y.cell.size<-(out.extent['y','max']-out.extent['y','min'])/(out.dim[2]-1)
  check <- writeRaster(raster(matrix(NWP.dem.out,ncol=out.dim[1]),  xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                              ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                       filename=file.path(out.dir,'NWP_dem.tif'), format='GTiff',datatype='INT2S', overwrite=T,
                       options=c("COMPRESS=LZW","PREDICTOR=2"))
  if (class(check)=='try-error'){print("ERROR:Writing file to disk failed.EXIT!");return(1)}
  remove(NWP.dem,idx)
  gc()
  
  dem<-list(
    y=seq(dem.extent['x','min'],dem.extent['x','max'],length=dem@grid@cells.dim[1]),
    x=seq(dem.extent['y','max'],dem.extent['y','min'],length=dem@grid@cells.dim[2]),
    z=matrix(dem@data$band1,ncol=dem@grid@cells.dim[1],byrow=T)*dem.scale+dem.offset
  )
  
  pts<-expand.grid(seq(out.extent['y','max'],out.extent['y','min'],length=out.dim[2]),
                   seq(out.extent['x','min'],out.extent['x','max'],length=out.dim[1]))
  pts<-data.frame(x=pts[,2],y=pts[,1])
  coordinates(pts)=~x+y
  proj4string(pts)=out.proj
  pts = spTransform(pts,dem.proj)
  idx<-pts@coords[,1]>=dem.extent['x','min'] & pts@coords[,1]<=dem.extent['x','max'] & 
    pts@coords[,2]>=dem.extent['y','min'] & pts@coords[,2]<=dem.extent['y','max']
  
  if (length(which(idx))==0){
    print('dem.extent and out.extent do not overlap.EXIT!')
    return(1)
  }
  
  dem.out<-rep(-32768,out.dim[1]*out.dim[2])
  dem.out[idx]<-interp.surface(dem,cbind(pts@coords[idx,2],pts@coords[idx,1]))
  gc()

  check <- writeRaster(raster(matrix(dem.out,ncol=out.dim[1]),xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                              ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                       filename=file.path(out.dir,'dem.tif'), format='GTiff',datatype='INT2S', overwrite=T,
                       options=c("COMPRESS=LZW","PREDICTOR=2"))
  if (class(check)=='try-error'){print("ERROR:Writing file to disk failed.EXIT!");return(1)}
  cor<-6*(NWP.dem.out-dem.out)/100
  cor[NWP.dem.out == -32768 | dem.out == -32768]<- -32768
  check <- writeRaster(raster(matrix(cor,ncol=out.dim[1]),xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                              ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                       filename=file.path(out.dir,'NWP_correction.tif'), format='GTiff',datatype='INT2S', overwrite=T,
                       options=c("COMPRESS=LZW","PREDICTOR=2"))
  if (class(check)=='try-error'){print("ERROR:Writing file to disk failed.EXIT!");return(1)}
  return(0)
}
