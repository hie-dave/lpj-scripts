library(RCurl)

sites        <- c("AdelaideRiver","AliceSpringsMulga","Litchfield","Tumbarumba")
level        <- 6  # processing level
processing   <- 'default'  # 'default' or 'site_pi'
version      <- '2022_v1' 
destdir      <- getwd()

for (site in sites){

  filelist <- getURL(url=paste0('https://dap.tern.org.au/thredds/catalog/ecosystem_process/ozflux/',site,'/',
                               version,'/L',level,'/',processing,'/','catalog.html'),
                    dirlistonly=TRUE)
  date <- substr(sub(paste0(".*",site),'',filelist),5,21)
  
  #cat('Downloading site:',site,'(Dates:',date,')',fill=F)
  
  download.file(url=paste0('http://dap.tern.org.au/thredds/fileServer/ecosystem_process/ozflux/',
                          site,'/',version,'/L',level,'/',processing,'/',site,'_L',level,'_',date,'.nc'),
                destfile=paste0(destdir,'/',site,'_',processing,'_',version,'_L',level,'_',date,'.nc'))

}
