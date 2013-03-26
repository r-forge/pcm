sunPosition <- function(year, month, day, hour=hour, min=min, sec=sec,
                        lat=lat, long=long) {
  twopi <- 2 * pi
  deg2rad <- pi / 180
  
  # Get day of the year, e.g. Feb 1 = 32, Mar 1 = 61 on leap years
  month.days <- c(0,31,28,31,30,31,30,31,31,30,31,30)
  day <- day + cumsum(month.days)[month]
  leapdays <- year %% 4 == 0 & (year %% 400 == 0 | year %% 100 != 0) & day >= 60
  day[leapdays] <- day[leapdays] + 1
  
  # Get Julian date - 2400000
  hour <- hour + min / 60 + sec / 3600 # hour plus fraction
  delta <- year - 1949
  leap <- trunc(delta / 4) # former leapyears
  jd <- 32916.5 + delta * 365 + leap + day + hour / 24
  
  # The input to the Atronomer's almanach is the difference between
  # the Julian date and JD 2451545.0 (noon, 1 January 2000)
  time <- jd - 51545.
  
  # Ecliptic coordinates
  
  # Mean longitude
  mnlong <- 280.460 + .9856474 * time
  mnlong <- mnlong %% 360
  mnlong[mnlong < 0] <- mnlong[mnlong < 0] + 360
  
  # Mean anomaly
  mnanom <- 357.528 + .9856003 * time
  mnanom <- mnanom %% 360
  mnanom[mnanom < 0] <- mnanom[mnanom < 0] + 360
  mnanom <- mnanom * deg2rad
  
  # Ecliptic longitude and obliquity of ecliptic
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong[eclong < 0] <- eclong[eclong < 0] + 360
  oblqec <- 23.429 - 0.0000004 * time
  eclong <- eclong * deg2rad
  oblqec <- oblqec * deg2rad
  
  # Celestial coordinates
  # Right ascension and declination
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num / den)
  ra[den < 0] <- ra[den < 0] + pi
  ra[den >= 0 & num < 0] <- ra[den >= 0 & num < 0] + twopi
  dec <- asin(sin(oblqec) * sin(eclong))
  
  # Local coordinates
  # Greenwich mean sidereal time
  gmst <- 6.697375 + .0657098242 * time + hour
  gmst <- gmst %% 24
  gmst[gmst < 0] <- gmst[gmst < 0] + 24.
  
  # Local mean sidereal time
  lmst <- gmst + long / 15.
  lmst <- lmst %% 24.
  lmst[lmst < 0] <- lmst[lmst < 0] + 24.
  lmst <- lmst * 15. * deg2rad
  
  # Hour angle
  ha <- lmst - ra
  ha[ha < -pi] <- ha[ha < -pi] + twopi
  ha[ha > pi] <- ha[ha > pi] - twopi
  
  # Latitude to radians
  lat <- lat * deg2rad
  
  # Azimuth and elevation
  el <- asin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha))
  az <- asin(-cos(dec) * sin(ha) / cos(el))
  elc <- asin(sin(dec) / sin(lat))
  az[el >= elc] <- pi - az[el >= elc]
  az[el <= elc & ha > 0] <- az[el <= elc & ha > 0] + twopi
  
  el <- el / deg2rad
  az <- az / deg2rad
  lat <- lat / deg2rad
  
  return(list(elevation=el, azimuth=az))
}

esdist_correct<-function(year,month,day,img){
  esdist_factor<-function(year,month,day){
    nday = c(1,32,60,91,121,152,182,213,244,274,305,335,366)
    dratio = c(1.034,1.030,1.019,1.001,.985,.972,.967,.971,.982,.998,1.015,1.029,1.034)
    dratio = 1/dratio
    date_text<-paste(formatC(year,width = 4, format = "d", flag = "0"),formatC(month,width = 2, format = "d", flag = "0"),
                     formatC(day,width = 2, format = "d", flag = "0"),sep='-')
    iday<-as.numeric(format(as.Date(date_text),'%j'))
    factor<-approx(nday, dratio, iday, method="linear",rule = 2)
    return(factor$y)
  }
  factor<-esdist_factor(year,month,day)
  img = img * factor
  ilow = which(img < 0)
  if (length(ilow) != 0){img[ilow] = 0}
  ihigh = which(img > 2,arr.ind=T)
  if (length(ihigh) != 0){img[ihigh] = 2}
  return(img)
}

centralwave_avhrr<-function(satname,t,ch){
  #KLM User's Guide (http://www2.ncdc.noaa.gov/docs/klm/html/c7/sec7-1.htm,links to tables, we used the "Moment center wavenumber").
  #http://www.eumetsat.int/Home/Main/Satellites/Metop/Resources/index.htm?l=en, links to callibration files, we used the "Moment center wavenumber"
  cws <- list(noaa7=list(lt225=c(ch3=2670.930,ch4=926.2000),f255t275=c(ch3=2670.300,ch4=926.8000),gt275=c(ch3=2671.90,ch4=927.2200)),
              noaa9=list(lt225=c(ch3=2670.930,ch4=928.5000),f255t275=c(ch3=2674.810,929.0200),gt275=c(ch3=2678.11,ch4=929.4600)),
              noaa11=list(lt225=c(ch3=2663.500,ch4=926.8100),f255t275=c(ch3=2668.150,ch4=927.3600),gt275=c(ch3=2671.40,ch4=927.8300)),
              noaa12=list(lt225=c(ch3=2632.713,ch4=920.0158),f255t275=c(ch3=2636.669,ch4=920.5504),gt275=c(ch3=2639.61,ch4=921.0291)),
              noaa14=list(lt225=c(ch3=2638.652,ch4=928.2603),f255t275=c(ch3=2642.807,ch4=928.8284),gt275=c(ch3=2645.90,ch4=929.3323)),
              noaa15=list(lt225=c(ch3=2694.8241,ch4=925.7260),f255t275=c(ch3=2694.8241,ch4=925.7260),gt275=c(ch3=2694.8241,ch4=925.7260)),
              noaa16=list(lt225=c(ch3=2698.6134,ch4=918.0988),f255t275=c(ch3=2698.6134,ch4=918.0988),gt275=c(ch3=2698.6134,ch4=918.0988)),
              noaa17=list(lt225=c(ch3=2670.7832,ch4=927.0246),f255t275=c(ch3=2670.7832,ch4=927.0246),gt275=c(ch3=2670.7832,ch4=927.0246)),
              noaa18=list(lt225=c(ch3=2663.5665,ch4=927.0246),f255t275=c(ch3=2663.5665,ch4=927.0246),gt275=c(ch3=2663.5665,ch4=927.0246)),
              noaa19=list(lt225=c(ch3=2671.6576,ch4=927.8562),f255t275=c(ch3=2671.6576,ch4=927.8562),gt275=c(ch3=2671.6576,ch4=927.8562)),
              metopA=list(lt225=c(ch3=2689.9123,ch4=926.1104),f255t275=c(ch3=2689.9123,ch4=926.1104),gt275=c(ch3=2689.9123,ch4=926.1104)),
              metopB=list(lt225=c(ch3=2666.9951,ch4=932.6465),f255t275=c(ch3=2666.9951,ch4=932.6465),gt275=c(ch3=2666.9951,ch4=932.6465)))
  
  isat<-which(names(cws) == satname)
  if (length(isat)==0){
    warning("Not supported satellite platform.Returnin central wave length for noaa19!", call. = F)
    isat<-which('noaa19' == names(cws))
  }
  ich <- ch - 2
  if (t < 225.){itemp = 1}else if ((t >= 225) & (t <= 275)){itemp = 2}else if (t >= 275){itemp = 3}
  
  return(cws[[isat]][[itemp]][ich])
}

ch3b_reflectance<-function(ch3,ch4,solzen,satname){
  
  c1 = 1.1910659e-5                 #; First radiation constant, mW/(m^2.sr.cm^-4)
  c2 = 1.438833                     #    ; Second radiation constant, cm K
  L0 = 3.416633306               		#; Weighted solar flux based on mean E-S distance
  
  DEGRAD = pi/180.              	#; To convert degrees to radians
  #;---------------------------------------------------------------------
  
  ch3ref <- ch3
  ch3ref[]<-NA
  # Get the channel 3 central wavenumber for the satellite number passed in,
  # and for the average temperature in the image.  
  
  igood<-which(! (is.na(ch3) | is.na(ch4) | is.na(solzen)),arr.ind=T)
  tave <- mean(ch4[igood])
  cw3 = centralwave_avhrr(satname,tave,3)
  
  # Compute the reflectance.  Radiance units are mW/(m^2.sr.cm^-1).
  
  b3t4 = (c1*cw3^3)/(exp(c2*cw3/ch4[igood])-1.0)    	# Compute Planck radiance B3(T4)
  L3   = (c1*cw3^3)/(exp(c2*cw3/ch3[igood])-1.0)    	# Convert Tb3 to radiance
  ch3ref[igood] = (L3 - b3t4) / (L0*cos(solzen[igood]*DEGRAD) - b3t4)
  
  #  Reset out-of-range values, could happen if solzen is greater than 90 for nigtht.
  
  ilow = which(ch3ref < 0)
  if (length(ilow) != 0){ch3ref[ilow] = 0}
  ihigh = which(ch3ref > 1)
  if (length(ihigh) != 0){ch3ref[ihigh] = 1}
  return(ch3ref)
}

get.knnx.cm<-function(lut_vals,vals){
  idx_fnn<-get.knnx(lut_vals,vals,k=1,algorithm="kd_tree")$nn.index
  idx_change<-(lut_vals[idx_fnn]-vals) > 0
  idx_fnn[idx_change]<-idx_fnn[idx_change]-1
  idx_fnn[idx_fnn == 0]<-1
  return(idx_fnn)
}

PCM<-function(avhrr.file=avhrr.file,avhrr.band.order=avhrr.band.order,
              angles.file=angles.file,angles.band.order=angles.band.order,lc.file=lc.file,
              out.dir=out.dir,satellite=satellite,year=year,month=month,day=day,
              hour=hour,minute=minute,dem.file=dem.file,skt.file=skt.file,
              avhrr.scale=avhrr.scale,avhrr.offset=avhrr.offset,angles.scale=angles.scale,angles.offset=angles.offset,
              skt.scale=skt.scale,skt.offset=skt.offset,dem.scale=dem.scale,dem.offset=dem.offset,
              avhrr.noData=avhrr.noData,angles.noData=angles.noData,lc.noData=lc.noData,skt.noData=skt.noData,
              dem.noData=dem.noData,backingfile=backingfile,ahvrr.normalise=ahvrr.normalise,
              ahvrr.distance=ahvrr.distance,postclassification=postclassification,
              sunazi.satazi=sunazi.satazi,satazi.sunazi=satazi.sunazi){
  
  options(bigmemory.typecast.warning=FALSE)
  options(bigmemory.allow.dimnames=TRUE)
  mis_val<- 65535
  day_angle<-85
  snow_thr<-285
  
  #sanity check
  if (missing('avhrr.file')){print("ERROR:avhrr.file variable not specified. EXIT!");return(1)}
  if (missing('avhrr.band.order')){print("ERROR:avhrr.band.order variable not specified. EXIT!");return(1)}
  if (missing('angles.file')){print("ERROR:angles.file variable not specified. EXIT!");return(1)}
  if (missing('angles.band.order')){print("ERROR:angles.band.order variable not specified. EXIT!");return(1)}
  if (missing('lc.file')){print("ERROR:lc.file variable not specified. EXIT!");return(1)}
  if (missing('out.dir')){print("ERROR:out.dir variable not specified. EXIT!");return(1)}
  if (missing('satellite')){print("ERROR:satellite variable not specified. EXIT!");return(1)}
  if (missing('year')){print("ERROR:year variable not specified. EXIT!");return(1)}
  if (missing('month')){print("ERROR:month variable not specified. EXIT!");return(1)}
  if (missing('day')){print("ERROR:day variable not specified. EXIT!");return(1)}
  if (missing('hour')){print("ERROR:hour variable not specified. EXIT!");return(1)}
  if (missing('minute')){print("ERROR:minute variable not specified. EXIT!");return(1)}
  if (missing('skt.file')){print("WARNING:SKT.file is missing. Results will be of lower quality!");flag.skt<-F}else{flag.skt<-T}
  if (flag.skt & missing('dem.file')){print("ERROR:dem.file variable not specified. EXIT!");return(1)}
  if (mode(avhrr.file)!='character'){print("ERROR:avhrr.file is not a character string. EXIT!");return(1)}
  if (mode(avhrr.band.order)!='character'){print("ERROR:avhrr.band.order is not a character string vector. EXIT!");return(1)}  
  if (mode(angles.file)!='character'){print("ERROR:angles.file is not a character string. EXIT!");return(1)} 
  if (mode(angles.band.order)!='character'){print("ERROR:angles.band.order is not a character string vector. EXIT!");return(1)}  
  if (mode(lc.file)!='character'){print("ERROR:lc.file is not a character string. EXIT!");return(1)}  
  if (flag.skt){if(mode(dem.file)!='character'){print("ERROR:dem.file is not a character string. EXIT!");return(1)}}
  if (flag.skt){if(mode(skt.file)!='character'){print("ERROR:skt.file is not a character string. EXIT!");return(1)}} 
  if (mode(out.dir)!='character'){print("ERROR:out.dir is not a character string. EXIT!");return(1)} 
  
  if (! file.exists(avhrr.file)){print(paste("ERROR:avhrr.file does not exists:",avhrr.file));return(1)}
  if (! file.exists(angles.file)){print(paste("ERROR:angles.file does not exists:",angles.file));return(1)}
  if (! file.exists(lc.file)){print(paste("ERROR:lc.file does not exists:",lc.file));return(1)}
  if (flag.skt){if(! file.exists(dem.file)){print(paste("ERROR:dem.file does not exists:",dem.file));return(1)}}
  if (flag.skt){if(! file.exists(skt.file)){print(paste("ERROR:skt.file does not exists:",skt.file));return(1)}}
  if (! file.exists(out.dir)){print(paste("ERROR:out.dir does not exists:",out.dir));return(1)}
  
  if (length(which(! avhrr.band.order %in% c('ch1','ch2','ch3a','ch3b','ch4','ch5','empty'))) !=0){
    print(paste("ERROR:Specified avhrr.band.order:",avhrr.band.order,"not in a form c('ch1','ch2','ch3a','ch3b','ch4','ch5','empty').EXIT!"))
    return(1)
  }
  if (length(which(! angles.band.order %in% c('sunz','satz','razi','sunazi','satazi'))) !=0){
    print(paste("ERROR:Specified angles.band.order:",angles.band.order,"not in a form c('sunz','satz','razi','sunazi','satazi').EXIT!"))
    return(1)
  }    
  
  idx.avhrr<-which(avhrr.band.order!='empty')
  idx.angles<-which(angles.band.order!='empty')
  nbands.avhrr<-length(avhrr.band.order)
  nbands.angles<-length(angles.band.order)
  
  if (missing('avhrr.scale')){avhrr.scale<-rep(1,nbands.avhrr)}
  if (missing('avhrr.offset')){avhrr.offset<-rep(0,nbands.avhrr)}
  if (missing('angles.scale')){angles.scale<-rep(1,nbands.angles)}
  if (missing('angles.offset')){angles.offset<-rep(0,nbands.angles)}
  if (missing('skt.scale')){skt.scale<-1}
  if (missing('skt.offset')){skt.offset<-0}
  if (missing('dem.scale')){dem.scale<-1}
  if (missing('dem.offset')){dem.offset<-0}
  if (missing('avhrr.noData')){avhrr.noData<-rep(NA,nbands.avhrr)}
  if (missing('angles.noData')){angles.noData<-rep(NA,nbands.angles)}
  if (missing('lc.noData')){lc.noData<-NA}
  if (missing('skt.noData')){skt.noData<-NA}  
  if (missing('dem.noData')){dem.noData<-NA}
  if (missing('ahvrr.normalise')){ahvrr.normalise<-rep(F,nbands.avhrr)}
  if (missing('ahvrr.distance')){ahvrr.distance<-rep(F,nbands.avhrr)}
  if (missing('postclassification')){postclassification<-F}
  if (missing('satazi.sunazi')){satazi.sunazi<-F}
  if (missing('sunazi.satazi')){sunazi.satazi<-F}
  
  avhrr.noData[unlist(attributes(GDALinfo(avhrr.file))$df['NoDataValue']) == avhrr.noData]<-NA
  angles.noData[unlist(attributes(GDALinfo(angles.file))$df['NoDataValue']) == angles.noData]<-NA
  lc.noData[unlist(attributes(GDALinfo(lc.file))$df['NoDataValue']) == lc.noData]<-NA
  if(flag.skt){skt.noData[unlist(attributes(GDALinfo(skt.file))$df['NoDataValue']) == skt.noData]<-NA}
  if(flag.skt){dem.noData[unlist(attributes(GDALinfo(dem.file))$df['NoDataValue']) == dem.noData]<-NA}
  
  if (! missing('backingfile')){
    descriptorfile=paste(basename(backingfile),'desc',sep='.')
    backingpath=dirname(backingfile)
    backingfile=basename(backingfile)
  }

  if (length(idx.avhrr)!=length(avhrr.band.order)){
    nbands.avhrr<-length(idx.avhrr)
    avhrr.band.order<-avhrr.band.order[idx.avhrr]
    avhrr.noData<-avhrr.noData[idx.avhrr]
    avhrr.scale<-avhrr.scale[idx.avhrr]
    avhrr.offset<-avhrr.offset[idx.avhrr]
    ahvrr.normalise<-ahvrr.normalise[idx.avhrr]
    ahvrr.distance<-ahvrr.distance[idx.avhrr]
  }
  
  if (length(idx.angles)!=length(angles.band.order)){
    nbands.angles<-length(idx.angles)
    angles.band.order<-angles.band.order[idx.angles]
    angles.noData<-angles.noData[idx.angles]
    angles.scale<-angles.scale[idx.angles]
    angles.offset<-angles.offset[idx.angles]  
  }

  flag.3a<-length(which(avhrr.band.order=='ch3a'))!=0
  flag.3b<-length(which(avhrr.band.order=='ch3b'))!=0
  
  ncols<-nbands.avhrr+nbands.angles+3  
  info<-GDALinfo(lc.file,silent=T)

  dims<-c(info['rows'],info['columns'])
  x.cell.size<-info['res.x']
  y.cell.size<-info['res.y']
  out.extent<-matrix(c(info['ll.x']+x.cell.size/2,info['ll.x']+x.cell.size/2 +1000*(dims[2]-1),
                       info['ll.y']+y.cell.size/2,info['ll.y']+y.cell.size/2+ 1000*(dims[1]-1)),
                     nrow=2,dimnames=list(c('x','y'),c('min','max')),byrow=T)
  out.proj<-attributes(info)$projection
  print('out.extent:')
  print(out.extent)   
  if (! missing('backingfile')){
    data<-big.matrix(info['rows']*info['columns'],ncols,shared=T,backingfile=backingfile,descriptorfile=descriptorfile,backingpath=backingpath,init=0)
  }else{
    data<-big.matrix(info['rows']*info['columns'],ncols,shared=F,init=0)
  }
  gc()
  data[,seq(nbands.avhrr)]<-as.matrix(readGDAL(avhrr.file,band=idx.avhrr)@data)
  gc()
  data[,seq(nbands.angles)+nbands.avhrr]<-as.matrix(readGDAL(angles.file,band=idx.angles)@data)
  gc()
  data[,nbands.angles+nbands.avhrr+1]<-matrix(readGDAL(lc.file)@data$band1,ncol=1)
  if(flag.skt){
      data[,nbands.angles+nbands.avhrr+2]<-matrix(readGDAL(dem.file)@data$band1,ncol=1)
      data[,nbands.angles+nbands.avhrr+3]<-matrix(readGDAL(skt.file)@data$band1,ncol=1)
  }
  
  if (flag.3a){ch3a.noData<-avhrr.noData[avhrr.band.order == 'ch3a']}
  if (flag.3b){ch3b.noData<-avhrr.noData[avhrr.band.order == 'ch3b']}
  avhrr.noData[avhrr.band.order == 'ch3b' | avhrr.band.order == 'ch3a']<- -100000000 #some big value to test
  
  idx.complete<-mwhich(data,seq(ncols),as.list(c(avhrr.noData,angles.noData,lc.noData,dem.noData,skt.noData)),as.list(rep('neq',ncols)),op='AND')
  gc()
  if(length(idx.complete)!=info['rows']*info['columns']){
      gc()
      if (! missing('backingfile')){
        data<-as.big.matrix(data[idx.complete,],backingfile=backingfile,descriptorfile=descriptorfile,backingpath=backingpath)
      }else{
        data<-as.big.matrix(data[idx.complete,],shared=F)        
      }
  }
  colnames(data)<-c(avhrr.band.order,angles.band.order,'lc','dem','skt')
  
  if (flag.3a){idx.3a<-data[,'ch3a']!=ch3a.noData}
  if (flag.3b){idx.3b<-data[,'ch3b']!=ch3b.noData}
  
  gc()
  data[]<-t(t(data[])*c(avhrr.scale,angles.scale,1,dem.scale,skt.scale)+
                       c(avhrr.offset,angles.offset,0,dem.offset,skt.offset))
  
  if (satazi.sunazi){
    data[,'satazi']<-abs(data[,'satazi']-data[,'sunazi'])
    colnames(d)[colnames(d)=='satazi']<-'razi'
  }
  
  if (sunazi.satazi){
    data[,'satazi']<-abs(data[,'sunazi']-data[,'satazi'])
    colnames(d)[colnames(d)=='satazi']<-'razi'
  }
  
  idx_day <- data[,'sunz'] < 89 
  idx<-which(ahvrr.normalise)
  if (length(idx)!=0){data[idx_day,idx]<-(data[idx_day,idx])/cos(data[idx_day,'sunz']*pi/180)}
  idx<-which(ahvrr.normalise)
  if (length(idx)!=0){
    for (i in idx){data[idx_day,i]<-esdist_correct(year,month,day,data[idx_day,i])}
  }
  idx_night <- ! idx_day
  idx_water<-data[,'lc'] == 20
  idx_land<- ! data[,'lc'] %in% c(20,0)
  
  if (flag.skt){
    data[,'skt']<-data[,'skt']-data[,'ch4']
    colnames(data)<-c(avhrr.band.order,angles.band.order,'lc','dem','sktmch4')
    gc()
    sktmch4_idx<- data[,'ch4']>290 | (idx_land & idx_day & data[,'ch1'] < 0.15) | data[,'dem'] > 2500 | 
       (data[,'sktmch4']<16 & data[,'dem'] >= 1200) | (data[,'sktmch4']<9 & idx_water) | (data[,'lc'] == 0 & idx_night) | 
       (data[,'lc'] == 0 & idx_day & data[,'ch1'] < 0.3)
    gc()
    if (postclassification){
      sktmch4_copy<-as.big.matrix(matrix(data[,'sktmch4'],ncol=1),shared=F)
    }
    data[sktmch4_idx,'sktmch4'] <-0
    remove(sktmch4_idx)
    gc()
  }else{
    colnames(data)<-c(avhrr.band.order,angles.band.order,'lc','dem','sktmch4')
    if (postclassification){
      sktmch4_copy<-as.big.matrix(matrix(rep(20,length(idx.complete)),ncol=1),shared=F)
    }
  }

  idx_sunglint<-acos(sin(data[,'satz']*pi/180)*sin(data[,'sunz']*pi/180)*cos(data[,'razi']*pi/180)+
                       cos(data[,'satz']*pi/180)*cos(data[,'sunz']*pi/180))*180/pi
  gc()
  idx_sunglint<-idx_sunglint > 0 & idx_sunglint < 36 & data[,'sunz'] < 85 
  data[idx_sunglint & (idx_water | data[,'lc'] == 0),'lc']<- -1
  data[idx_sunglint & data[,'lc'] >= 24 & data[,'lc'] <= 27,'lc']<- -2
  remove(idx_sunglint)
  gc()
  
  data[idx_water,'ch1']<-data[idx_water,'ch2']
  texture<-rep(0,dims[1]*dims[2])
  texture[idx.complete]<-data[,'ch1']*100
  gc()
  texture <-c(imgConvolve(imagedata(texture,type='grey',ncol=dims[2],nrow=dims[1]),shoreline_kernel(3),0))[idx.complete]/255
  texture[! idx_water]<-0
  gc()
  #dtn
  dtn<-as.numeric(idx_day)
  dtn[data[,'sunz'] > day_angle]<-2
  dtn[idx_night]<-3
  gc()
  
  #ch4texture
  ch4_texture<-rep(0,dims[1]*dims[2])
  ch4_texture[idx.complete]<-abs(data[,'ch4']-230)
  gc()
  ch4_texture <-c(imgConvolve(imagedata(ch4_texture,type='grey',ncol=dims[2],nrow=dims[1]),shoreline_kernel(3),0))[idx.complete]/50
  texture[data[,'sunz'] >= 88.95 & idx_water]<-ch4_texture[data[,'sunz'] >= 88.95 & idx_water]
  remove(ch4_texture)
  gc()

  prob<-big.matrix(dims[1]*dims[2],1,init=mis_val,type='integer',shared=F)

  if (flag.3a){
    if (length(which(idx.3a))!=0){
        ESF1<-as.big.matrix(matrix(data[idx.3a,'sktmch4']/20+0.1,ncol=1),shared=F)
        if (length(which(idx_day & idx.3a))!=0){
          S1<-rbind(sktmch4=c(156.7073163,0.75528466),
                    ch1=c(0.75528466,0.04312983))
          colnames(S1)<-c('sktmch4','ch1')
          S2<-rbind(sktmch4=c(245.3430733,-0.32366008),
                    ch1=c(-0.32366008,0.04665813))
          colnames(S2)<-c('sktmch4','ch1')
          gc()
          ESF1[idx_day[idx.3a]]<-ics(data.frame(sktmch4=data[idx_day&idx.3a,'sktmch4'],ch1=data[idx_day&idx.3a,'ch1']),S1=S1,S2=S2,stdB="B")@Scores$IC.2
        }
        gc()
        
        ESF2<-big.matrix(length(ESF1),1,init=0,shared=F) 
        if (length(which(idx_day & idx.3a))!=0){
          S1<-rbind(sktmch4=c(156.7073163,0.75528466),
                    ch3=c(0.75528466,0.04312983))
          colnames(S1)<-c('sktmch4','ch3')
          S2<-rbind(sktmch4=c(245.3430733,-0.32366008),
                    ch3=c(-0.32366008,0.04665813))
          colnames(S2)<-c('sktmch4','ch3')     
          ESF2[idx_day[idx.3a]]<-ics(data.frame(sktmch4=data[idx_day & idx.3a,'sktmch4'],ch3=data[idx_day & idx.3a,'ch3a']),S1=S1,S2=S2,stdB="B")@Scores$IC.2
        }
        
        gc()
        .env <- environment()
        LUT3a<-NULL
        data('LUT3a',envir = .env)
        idx_array<-big.matrix(length(ESF1),8,type='integer',shared=F,init=1)
        idx_array[,1]<-get.knnx(LUT3a$lc.vals,data[idx.3a,'lc'],k=1,algorithm="kd_tree")$nn.index
        idx_array[,2]<-dtn[idx.3a]
        #remove(dtn)
        gc()
        idx_array[,3]<-get.knnx.cm(LUT3a$text.vals,texture[idx.3a])
        gc()
        idx_array[,4]<-get.knnx.cm(LUT3a$F1.vals,ESF1[])
        remove(ESF1)
        gc()
        idx_array[,5]<-get.knnx.cm(LUT3a$F2.vals,ESF2[])
        remove(ESF2)
        gc()
        idx_array[,6]<-get.knnx.cm(LUT3a$F3.vals,data[idx.3a,'ch4']-data[idx.3a,'ch5'])  
        gc()
        idx_array[,7]<-get.knnx(LUT3a$view.vals,data[idx.3a,'satz'],k=1,algorithm="kd_tree")$nn.index
        gc()
        if (length(which(idx_day))!=0){idx_array[idx_day[idx.3a],8]<-get.knnx(LUT3a$azi.vals,abs(data[idx_day&idx.3a,'razi']),k=1,algorithm="kd_tree")$nn.index}
        gc()
        prob[idx.complete][idx.3a]<-LUT3a$LUT[idx_array[]]
        prob[idx.complete][idx_night&idx.3a]<-mis_val
        remove(idx_array,LUT3a)
        gc()

    }
  }
  if (flag.3b){
    if (length(which(idx.3b))!=0){
        
        ESF1<-as.big.matrix(matrix(data[idx.3b,'sktmch4']/20+0.1,ncol=1),shared=F)
        if (length(which(idx_day&idx.3b))!=0){
          S1<-rbind(sktmch4=c(156.7073163,0.75528466),
                    ch1=c(0.75528466,0.04312983))
          colnames(S1)<-c('sktmch4','ch1')
          S2<-rbind(sktmch4=c(245.3430733,-0.32366008),
                    ch1=c(-0.32366008,0.04665813))
          colnames(S2)<-c('sktmch4','ch1')
          gc()
          ESF1[idx_day[idx.3b]]<-ics(data.frame(sktmch4=data[idx_day & idx.3b,'sktmch4'],ch1=data[idx_day & idx.3b,'ch1']),S1=S1,S2=S2,stdB="B")@Scores$IC.2
        }
        gc()
        ESF2<-as.big.matrix(matrix((data[idx.3b,'ch4']-data[idx.3b,'ch3b'])/-19+0.2,ncol=1),shared=F)
        idx<-which(data[idx.3b,'sunz']<= day_angle)
        gc()
        if (length(idx)!=0){
          ch3ref<-ch3b_reflectance(data[idx.3b,'ch3b'][idx],data[idx.3b,'ch4'][idx],data[idx.3b,'sunz'][idx],satellite)
          S1<-rbind(sktmch4=c(156.7073163,0.75528466),
                    ch3=c(0.75528466,0.04312983))
          colnames(S1)<-c('sktmch4','ch3')
          S2<-rbind(sktmch4=c(245.3430733,-0.32366008),
                    ch3=c(-0.32366008,0.04665813))
          colnames(S2)<-c('sktmch4','ch3')  
          ESF2[idx]<-ics(data.frame(sktmch4=data[idx.3b,'sktmch4'][idx],ch3=ch3ref),S1=S1,S2=S2,stdB="B")@Scores$IC.2
        }
        gc()
        .env <- environment()
        LUT3b<-NULL
        data('LUT3b',envir=.env)
        idx_array<-big.matrix(length(ESF1),8,type='integer',shared=F,init=1)
        idx_array[,1]<-get.knnx(LUT3b$lc.vals,data[idx.3b,'lc'],k=1,algorithm="kd_tree")$nn.index
        idx_array[,2]<-dtn[idx.3b]
        remove(dtn)
        gc()
        idx_array[,3]<-get.knnx.cm(LUT3b$text.vals,texture[idx.3b])
        remove(texture)
        gc()
        idx_array[,4]<-get.knnx.cm(LUT3b$F1.vals,ESF1[])
        remove(ESF1)
        gc()
        idx_array[,5]<-get.knnx.cm(LUT3b$F2.vals,ESF2[])
        remove(ESF2)
        gc()
        idx_array[,6]<-get.knnx.cm(LUT3b$F3.vals,data[idx.3b,'ch4']-data[idx.3b,'ch5'])  
        gc()
        idx_array[,7]<-get.knnx(LUT3b$view.vals,data[idx.3b,'satz'],k=1,algorithm="kd_tree")$nn.index
        gc()
        if (length(which(idx_day))!=0){idx_array[idx_day[idx.3b],8]<-get.knnx(LUT3b$azi.vals,data[idx_day&idx.3b,'razi'],k=1,algorithm="kd_tree")$nn.index}
        gc()
        prob[idx.complete][idx.3b]<-LUT3b$LUT[idx_array[]]
        remove(idx_array,LUT3b)
        gc()
    }
  }
  if (!! snow_thr){
    idx_warm<-rep(F,dims[1]*dims[2])
    idx_warm[idx.complete]<-data[,'ch4'] >= snow_thr
    prob[idx_warm & prob[] > 50 & prob[] <150]<-200
    remove(idx_warm)
  }
  #remove(idx_day,nc)
  gc()
  
  txt<-unlist(strsplit(basename(avhrr.file),'\\.'))
  txt<-paste(txt[1:length(txt)-1],collapse=".")
    
  file_tif1<-file.path(out.dir,paste(txt,"temp.tif",sep="_"))
  file_tif<-file.path(out.dir,paste(txt,"prob.tif",sep="_"))
  file_vrt<-file.path(out.dir,paste(txt,'temp.vrt',sep="_"))
  
  color_table<-col2rgb(colorRampPalette(c(rgb(red=0,green=170/255,blue=0),rgb(red=0,green=255/255,blue=255/255),
                                          rgb(red=202/255,green=202/255,blue=202/255),rgb(red=0,green=170/255,blue=0)), space = "rgb")(301),alpha=F)
  text<-apply(color_table,2,function(x){paste('<Entry c1="',x['red'],'" c2="',x['green'],'" c3="',x['blue'],'" c4="255"/>',sep='')})
  color_table<-c('<ColorInterp>Palette</ColorInterp>','<ColorTable>',text,'</ColorTable>')
  
  rf <- writeRaster(raster(matrix(prob[],nrow=dims[1],byrow=T), xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                           ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                    filename=file_tif1, format='GTiff',datatype='INT2U', overwrite=TRUE)
  system(paste('gdal_translate -q -of VRT',file_tif1,file_vrt))
  d<-scan(file_vrt, character(0), sep = "\n",strip.white=T)
  d[d=="<ColorInterp>Gray</ColorInterp>"]<-paste(color_table,collapse='\n')
  cat(d,file=file_vrt,sep='\n')
  system(paste('gdal_translate -q ',file_vrt,file_tif,'-co COMPRESS=LZW -co PREDICTOR=2'))
  file.remove(c(file_tif1,file_vrt))
  if (! missing('backingfile')){
    file.remove(file.path(backingpath,backingfile))
    file.remove(file.path(backingpath,descriptorfile))
  }
  if (postclassification){
    #shadow
    pix_size<-c(x.cell.size,y.cell.size)
    clouds<-matrix(as.integer(prob[]>=170 & prob[] <=250),nrow=dims[2])
    gc()
    adjacent<-as.integer(imgConvolve(imagedata(abs(clouds-1),type='grey',ncol=dims[2],nrow=dims[1]),shoreline_kernel(3),0)) > 0
    clouds_edges<-as.integer(imgConvolve(imagedata(clouds,type='grey',ncol=dims[2],nrow=dims[1]),shoreline_kernel(3),0) > 0)[idx.complete]
    clouds_edges[idx_night]<-0
    idx_clouds_edges<-which(clouds_edges>0)
    remove(clouds_edges)
    gc()
    idx_water<-big.matrix(dims[1]*dims[2],1,init=2,shared=F)
    idx_water[idx.complete]<-as.numeric(! idx_land)

    clouds<-as.big.matrix(matrix(c(clouds)*10,ncol=1),type='short',shared=F)
    clouds[clouds[] == 0 & idx_water[] == 1] <- 1
    clouds[clouds[] == 0 & idx_water[] == 0] <- 2 
    clouds[prob[]>300]<-0
    clouds[prob[]>20 & prob[] <170]<-3
    clouds[clouds[] == 1 & adjacent] <-7
    clouds[clouds[] == 2 & adjacent] <-8
    clouds[clouds[] == 3 & adjacent] <-9
    x_out<-seq(from=out.extent['x','min'],to=out.extent['x','max'],length=dims[1]) 
    y_out<-seq(from=out.extent['y','max'],to=out.extent['y','min'],length=dims[2])
    remove(idx_water,idx_land,adjacent)
    yx<-as.big.matrix(cbind(expand.grid(x_out,y_out)[,2],c(t(matrix(expand.grid(y_out,x_out)[,2],nrow=dims[1],byrow=F))))[idx.complete,],type='integer',shared=F)
    gc()
    pts<-data.frame(x=yx[idx_clouds_edges,2],y=yx[idx_clouds_edges,1])
    coordinates(pts)=~x+y
    proj4string(pts)=CRS(out.proj)
    pts = spTransform(pts,CRS("+proj=longlat +datum=WGS84 +no_defs"))
    gc()
    
    i<-length(idx_clouds_edges)
    sun_azi<-sunPosition(year=rep(as.numeric(year),i),month=rep(as.numeric(month),i),day=rep(as.numeric(day),i),hour=rep(as.numeric(hour),i),
                         min=rep(as.numeric(minute),i),sec=rep(0,i),lat=pts@coords[,2],long=pts@coords[,1])$azimuth 
    idx_x<-sun_azi>180
    x_scale <-as.integer(as.numeric(idx_x)*2-1)
    sun_azi[idx_x]<-360-sun_azi[idx_x]
    idx_y<-sun_azi>90
    y_scale <- as.integer(as.numeric(!idx_y)*-2+1)
    sun_azi[idx_y]<-180-sun_azi[idx_y]
    remove(idx_x,idx_y,pts,x_out,y_out)
    gc()
    
    sktmch4<-big.matrix(dims[1]*dims[2],1,init=0,shared=F)
    sktmch4[idx.complete]<-sktmch4_copy[]
    sktmch4[clouds[] != 10 | sktmch4[] <0 ]<-0

    cloud_height<-(as.numeric(imgMaximumFilter(imagedata(sktmch4[],type='grey',ncol=dims[2],nrow=dims[1]),5))[idx.complete])*166.7
    res_array<-cbind(height=cloud_height[idx_clouds_edges],sunz=data[idx_clouds_edges,'sunz']*pi/180,sun_azi=sun_azi*pi/180,
                     x=yx[idx_clouds_edges,2],y=yx[idx_clouds_edges,1],x_scale=x_scale,y_scale=y_scale)
    remove(cloud_height,sktmch4,sun_azi,idx_clouds_edges,sktmch4_copy)
    gc()
    shadow_length<-round(res_array[,'height']/tan(90*pi/180-res_array[,'sunz'])+mean(pix_size),-3)
    shadow_length[shadow_length > 30*mean(pix_size)]<-30*mean(pix_size)
    gc()
    res_array<-as.big.matrix(cbind(res_array,shadow_length=shadow_length),shared=F)
    gc()
    
    array_expanded<-c()
    lengths<-unique(sort(res_array[,'shadow_length']))
    lengths<-lengths[lengths > mean(pix_size)]
    for (j in lengths){
      idx<-which(res_array[,'shadow_length']==j)
      intervals<-seq(mean(pix_size),j,by=mean(pix_size)/2)
      array_expanded<-as.big.matrix(rbind(array_expanded[],cbind(matrix(unlist(res_array[idx,rep(c('height','sunz','sun_azi','x','y','x_scale','y_scale'),each=length(intervals))]),ncol=7,byrow=F),rep(intervals,each=length(idx)))),shared=F)
    }
    res_array<-as.big.matrix(rbind(res_array[],array_expanded[]),shared=F)
    res_array<-cbind(y=res_array[,'y']+res_array[,'y_scale']*cos(res_array[,'sun_azi'])*res_array[,'shadow_length'],
                     x=res_array[,'x']+res_array[,'x_scale']*sin(res_array[,'sun_azi'])*res_array[,'shadow_length'])
    idx<-which(res_array[,'x']>=out.extent['x','min'] & res_array[,'x']<=out.extent['x','max'] & res_array[,'y']>=out.extent['y','min'] & 
                 res_array[,'y']<=out.extent['y','max'])
    remove(array_expanded)
    gc()
    res_array<-as.big.matrix(res_array[idx,],shared=F)
    remove(shadow_length,idx_night,y_scale,x_scale,idx)
    gc()
    shadows<-get.knnx(yx[],res_array[],k=2,algorithm="kd_tree")$nn.index
    gc()
    shad<-rep(0,dims[1]*dims[2])
    shad[idx.complete][shadows]<-1
    shad[idx.complete][data[,'ch2']>0.3]<-0
    
    clouds[clouds[] == 1 & shad==1] <-4
    clouds[clouds[] == 2 & shad==1] <-5
    clouds[clouds[] == 3 & shad==1] <-6
    
    file_tif1<-file.path(out.dir,paste(txt,"temp.tif",sep="_"))
    file_tif<-file.path(out.dir,paste(txt,"cmask.tif",sep="_"))
    file_vrt<-file.path(out.dir,paste(txt,'temp.vrt',sep="_"))
    rf <- writeRaster(raster(matrix(clouds[],nrow=dims[1],byrow=T), xmn=out.extent['x','min']-x.cell.size/2, xmx=out.extent['x','max']+x.cell.size/2, 
                             ymn=out.extent['y','min']-y.cell.size/2, ymx=out.extent['y','max']+y.cell.size/2, out.proj), 
                             filename=file_tif1, format='GTiff',datatype='INT1U', overwrite=TRUE)
    system(paste('gdal_translate -q -of VRT',file_tif1,file_vrt))
    d<-scan(file_vrt, character(0), sep = "\n",strip.white=T)
    color_table<-'<ColorInterp>Palette</ColorInterp>
<ColorTable>
<Entry c1="0" c2="0" c3="0" c4="255"/>
<Entry c1="0" c2="85" c3="255" c4="255"/>
<Entry c1="0" c2="170" c3="0" c4="255"/>
<Entry c1="0" c2="255" c3="255"  c4="255"/>
<Entry c1="0" c2="0" c3="211" c4="255"/>
<Entry c1="0" c2="85" c3="0" c4="255"/>
<Entry c1="85" c2="170" c3="255" c4="255"/>
<Entry c1="85" c2="0" c3="255" c4="255"/>
<Entry c1="85" c2="85" c3="0" c4="255"/>
<Entry c1="170" c2="170" c3="255" c4="255"/>
<Entry c1="202" c2="202" c3="202" c4="255"/>
</ColorTable>'
    d[d=="<ColorInterp>Gray</ColorInterp>"]<-color_table
    cat(d,file=file_vrt,sep='\n')
    system(paste('gdal_translate -q',file_vrt,file_tif,'-co COMPRESS=LZW -co PREDICTOR=2'))
    file.remove(c(file_tif1,file_vrt))
    remove(prob,clouds)
  }
  return(0)
}
