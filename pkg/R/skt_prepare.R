skt_prepare<-function(NWP.skt.file=NWP.skt.file,out.dir=out.dir,out.proj=out.proj,
                      out.extent=out.extent,out.dim=out.dim,correction.file=correction.file,NWP.skt.proj=NWP.skt.proj,
                      NWP.skt.extent=NWP.skt.extent,NWP.skt.half.cell=NWP.skt.half.cell,
                      NWP.skt.scale=NWP.skt.scale,NWP.skt.offset=NWP.skt.offset,
                      NWP.skt.band=NWP.skt.band){
  #sanity check
  if (missing('NWP.skt.file')){print("ERROR:NWP.skt.file variable not specified. EXIT!");return(1)}
  if (missing('out.dir')){print("ERROR:out.dir variable not specified. EXIT!");return(1)}
  if (missing('out.proj')){print("ERROR:out.proj variable not specified. EXIT!");return(1)}
  if (missing('out.extent')){print("ERROR:out.extent variable not specified. EXIT!");return(1)}
  if (missing('out.dim')){print("ERROR:out.dim variable not specified. EXIT!");return(1)}
  
  if (mode(NWP.skt.file)!='character'){print("ERROR:NWP.skt.file is not a character string. EXIT!");return(1)}
  if (mode(out.dir)!='character'){print("ERROR:out.dir is not a character string. EXIT!");return(1)}  
  if (mode(out.proj)!='character'){print("ERROR:out.proj is not a character string. EXIT!");return(1)} 
  
  out.extent<-as.numeric(out.extent)
  if (length(which(is.na(out.extent)))!=0 | length(out.extent) !=4){print("ERROR:out.extent is not a 4 elements numeric vector in a form c(ULX,ULY,LRX,LRY). EXIT!");return(1)} 
  out.dim<-as.numeric(out.dim)
  if (length(which(is.na(out.dim)))!=0 | length(out.dim) !=2){print("ERROR:out.dim is not a 2 elements numeric vector in a form c(ncol,nrows). EXIT!");return(1)}
  if (! file.exists(NWP.skt.file)){print(paste("ERROR:NWP.skt.file does not exists:",NWP.skt.file));return(1)}
  if (! file.exists(out.dir)){print(paste("ERROR:out.dir does not exists:",out.dir));return(1)}
  
  if (missing('NWP.skt.half.cell')){NWP.skt.half.cell<-c(0,1)}
  if(missing('NWP.skt.scale')){NWP.skt.scale<-1}
  if(missing('NWP.skt.offset')){NWP.skt.offset<-0}
  if(missing('NWP.skt.band')){NWP.skt.band<-1}
    
  NWP.skt<-try(readGDAL(NWP.skt.file,half.cell=NWP.skt.half.cell,band=NWP.skt.band))
  if (class(NWP.skt)=='try-error'){print("ERROR:GDAL could not open NWP.skt.file. EXIT!");return(1)}

  if (missing(NWP.skt.proj)){
    NWP.skt.proj<-NWP.skt@proj4string@projargs
  }
  NWP.skt.proj<-try(CRS(NWP.skt.proj))
  if (class(NWP.skt.proj)=='try-error'){print("ERROR:Specified proj4 arguments for the NWP.skt.proj are not valid. EXIT!");return(1)}
  
  if (missing(NWP.skt.extent)){
    NWP.skt.extent<-NWP.skt@bbox
    NWP.skt.extent['x','max']<-NWP.skt.extent['x','max']-NWP.skt@grid@cellsize[1]
    NWP.skt.extent['y','min']<-NWP.skt.extent['y','min']+NWP.skt@grid@cellsize[2]
  }else{
    if (NWP.skt.extent[1]>=NWP.skt.extent[3]){
      print(paste('ERROR:ULX=',NWP.skt.extent[1],'should be lower than LRX=',NWP.skt.extent[3],'.EXIT!'))
      return(1)
    }
    if (NWP.skt.extent[2]<=NWP.skt.extent[4]){
      print(paste('ERROR:ULY=',NWP.skt.extent[2],'should be greater than LRY=',NWP.skt.extent[4],'.EXIT!'))
      return(1)
    }    
    NWP.skt.extent<-matrix(c(NWP.skt.extent[1],NWP.skt.extent[4],NWP.skt.extent[3],NWP.skt.extent[2]),nrow=2,dimnames=list(c('x','y'),c('min','max')))
  }
  print('NWP.skt.extent:')
  print(NWP.skt.extent)
  
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
  
  #start processing  
  NWP.skt<-list(
    y=seq(NWP.skt.extent['x','min'],NWP.skt.extent['x','max'],length=NWP.skt@grid@cells.dim[1]),
    x=seq(NWP.skt.extent['y','max'],NWP.skt.extent['y','min'],length=NWP.skt@grid@cells.dim[2]),
    z=matrix(NWP.skt@data$band1,ncol=NWP.skt@grid@cells.dim[1],byrow=T)*NWP.skt.scale+NWP.skt.offset
  )
  
  pts<-expand.grid(seq(out.extent['y','max'],out.extent['y','min'],length=out.dim[2]),
                   seq(out.extent['x','min'],out.extent['x','max'],length=out.dim[1]))
  pts<-data.frame(x=pts[,2],y=pts[,1])
  coordinates(pts)=~x+y
  proj4string(pts)=out.proj
  pts = spTransform(pts,NWP.skt.proj)
  idx<-pts@coords[,1]>=NWP.skt.extent['x','min'] & pts@coords[,1]<=NWP.skt.extent['x','max'] & 
    pts@coords[,2]>=NWP.skt.extent['y','min'] & pts@coords[,2]<=NWP.skt.extent['y','max']
  
  if (length(which(idx))==0){
    print('NWP.skt.extent and out.extent do not overlap.EXIT!')
    return(1)
  }
  NWP.skt.out<-rep(65535,out.dim[1]*out.dim[2])
  NWP.skt.out[idx]<-interp.surface(NWP.skt,cbind(pts@coords[idx,2],pts@coords[idx,1]))*10
  if (! missing('correction.file')){
    if (file.exists(correction.file)){
        correction<-try(readGDAL(correction.file))
        if (class(correction)=='try-error'){print("ERROR:GDAL could not open correction.file. EXIT!");return(1)}
        correction<-c(matrix(correction@data$band1,ncol=correction@grid@cells.dim[1],byrow=T))
        correction[correction==-32768 | is.na(correction)]<-0
        NWP.skt.out<-NWP.skt.out+correction
        NWP.skt.out[! idx]<-65535
    }else{
      print(paste('WARNING:Specified correction file:',correction.file,'does not exist.'))
    }
  }
  NWP.skt.out[NWP.skt.out>3300 | NWP.skt.out<2000]<-65535
  gc()
  x.cell.size<-(out.extent['x','max']-out.extent['x','min'])/(out.dim[1]-1)
  y.cell.size<-(out.extent['y','max']-out.extent['y','min'])/(out.dim[2]-1)
  
  check <- writeRaster(raster(matrix(NWP.skt.out,ncol=out.dim[1]), xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                              ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                       filename=file.path(out.dir,'NWP_skt.tif'), format='GTiff',datatype='INT2U', overwrite=T,
                       options=c("COMPRESS=LZW","PREDICTOR=2"))
  if (class(check)=='try-error'){print("ERROR:Writing file to disk failed.EXIT!");return(1)}
  remove(NWP.skt,idx)
  gc()
  return(0)
}
