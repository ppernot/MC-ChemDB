# 1/ Gather samples
# 2/ Generate html index for reactions summaries
# 3/ Generate PDF doc of database

version   = '1.0'
rootDir   = '/home/pernot/Bureau/Titan-APSIS/MC-ChemDB/'
subDB     = 'Photo/'
sourceDir = paste0(rootDir,subDB,'Source/v_',version,'/')
tmpDir    = paste0(rootDir,subDB,'Tmp/v_',version,'/')
publicDir = paste0(rootDir,subDB,'Public/v_',version,'/')
docDir    = paste0(sourceDir,'Doc/')
webDir    = paste0(publicDir,'WWW/')

# Load libraries #####
library(stringr)
library(xtable)
library(bibtex)
library('CHNOSZ')
data('thermo')

# Load functions #####
setwd(sourceDir)
source('Scripts/stoechiometry.R') # Stoechiometry
printBib = function(keys,bibList) {
  if(length(keys) != 0) {
    cat('<H2>References</H2><DL>\n')
    for (key in keys) {
      cat(paste0('<DT>[',key,']</DT><DD>'))
      print(bibList[key],style="html")
      cat('</DD>')
    }
    cat('</DL>')
  }
}
printBibKeys = function(keys) {
  if(any(!is.na(keys))) {
    refs=paste0(sort(keys),collapse=',')
    paste0('[',refs,']')
  } else {
    NA
  }
}
printRQ = function(comments) {
  if(!is.na(comments)) {
    cat('<H2>Comments</H2>\n')
    cat(paste0('<font color="red">',comments,'</font>\n'))
  }
}

# Parameters
sampleSize=100 # Number max of random samples to gather
randomSamples=TRUE # Gather random databases (vs. nominal only)

setwd(paste0(tmpDir,'Reactions/')) # Contains generated samples
listReacs = list.dirs(full.names=FALSE, recursive=FALSE)
listReacs = gsub("./","",listReacs)

# Reorder by increasing mass of (1) ion & (2) other reactants
tabSpecies  = t(data.frame(sapply(listReacs, 
                                  function(x) strsplit(x,' + ',fixed=TRUE))))
mass=matrix(0,ncol=ncol(tabSpecies),nrow=nrow(tabSpecies))
for (i in 1:ncol(tabSpecies)) 
  mass[,i]=as.numeric(unlist(sapply(tabSpecies[,i],getMassList)))
tabSpecies = tabSpecies[order(mass[,1],mass[,2],na.last=FALSE),]
listReacs = apply(tabSpecies,1,function(x) paste0(x,collapse=' + '))
# listReacs=c('N2 + HV','H2 + HV')

# Clean target files to be appended to
indexFile=paste0(webDir,'index.html')
file.remove(indexFile)

spIndexFile=paste0(webDir,'spIndex.html')
file.remove(spIndexFile)

fileList=list.files(path=paste0(publicDir,'Databases'),full.name=TRUE)
file.remove(fileList)

dataTableFile=paste0(publicDir,'dataTable.html')
sink(file=dataTableFile,append=FALSE)
cat('<TABLE BORDER=0>')
sink(file=NULL)


allSpecies = allBibKeys = c()

for (reac in listReacs) {
  sink(file=dataTableFile,append=TRUE)
  cat('<TR><TD COLSPAN=5><HR size=1></TD></TR>\n')
  sink(file=NULL)
  file.append(
    file1=dataTableFile,
    file2=paste0(tmpDir,'Reactions/',reac,'/dataTable.html')) 
  
  # Generate Html index to summary files
  cat(paste0(reac,'\n'))
  sink(file=indexFile,append=TRUE)
  cat(paste0('<BR><A HREF="',webDir,reac,'/summary.html">',reac,'</A>\n')) 
  sink(file=NULL)
   
  # Generate collated Monte Carlo samples
  maxNum=ifelse(randomSamples,sampleSize,0)
  for (i in 0:maxNum) {
    # Append reactions list
    file.append(
      file1=paste0(publicDir,'Databases/run_',sprintf('%04i',i),'.csv'),
      file2=paste0(tmpDir,'Reactions/',reac,'/Samples/run_',
                   sprintf('%04i',i),'.csv')) 
    
    # Copy cross section and BR
    for(pat in c(paste0(sprintf('%04i',i),"_se*.dat$"),
                 paste0(sprintf('%04i',i),"_qy*.dat$")) )
    file.copy(
      from=  list.files(paste0(tmpDir,'Reactions/',reac,'/Samples'),
                        pattern = glob2rx(pat), 
                        full.names = TRUE),
      to  = paste0(publicDir,'Databases/')
    )
  }
  
  # Collate full species list
  species = read.csv(file=paste0(tmpDir,'Reactions/',reac,'/species.txt'),
                     sep=' ',header=FALSE, stringsAsFactors = FALSE)
  allSpecies=c(allSpecies,unlist(species))                   

  # Collate full biblio
  file=paste0(tmpDir,'Reactions/',reac,'/bibKeys.txt')
  if( file.info(file)$size != 0) {
    bibKeys = read.csv(file,sep=' ',header=FALSE,stringsAsFactors = FALSE)
    allBibKeys=c(allBibKeys,unlist(bibKeys))                   
  }
  
}
allSpecies = unique(allSpecies)
masses = sapply(allSpecies, getMassList)

sink(file=dataTableFile,append=TRUE)
cat('</TABLE>')
sink(file=NULL)

# Generate prod-loss file 
sink(file=spIndexFile,append=FALSE)

cat('<H2>Neutrals</H2>')
selIons=grepl('\\+$',allSpecies)
spec = allSpecies[!selIons]
mass  = masses[!selIons]
mord=order(mass)
for (sp in spec[mord]) {
  specDir=paste0(tmpDir,'Species/',sp)
  
  cat(paste0('<BR><B>',sp,'</B> ')) 
  pFile=paste0(specDir,'/prod.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Productions</A>'))
  
  pFile=paste0(specDir,'/loss.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Losses</A>'))
}

cat('<H2>Ions</H2>')
spec = allSpecies[selIons]
mass  = masses[selIons]
mord=order(mass)
for (sp in spec[mord]) {
  specDir=paste0(tmpDir,'Species/',sp)

  cat(paste0('<BR><B>',sp,'</B> ')) 
  pFile=paste0(specDir,'/prod.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Productions</A>'))

  pFile=paste0(specDir,'/loss.html')
  if(file.exists(pFile)) 
    cat(paste0(' <A HREF="',pFile,'">Losses</A>'))
}
sink(file=NULL)

targetHtml = paste0(sourceDir,'speciesList.html')
sink(file=targetHtml, append=FALSE)
cat('<H1>Species List</H1>')

# Dummies 
selAux = is.na(masses) | masses < 1
cat('<H2>Auxiliary species</H2>')
cat(paste0(allSpecies[selAux],collapse='<BR>'))

# Neutrals
trueSpecies=allSpecies[!selAux]
trueMasses=masses[!selAux]
selIons=grepl('\\+$',trueSpecies)
species = trueSpecies[!selIons]
spMass  = trueMasses[!selIons]
mord=order(spMass)
listSp = c()
for (i in seq_along(mord)) 
  listSp[i] = paste0('<font color="blue">',species[mord[i]],
                     '</font> (',signif(spMass[mord[i]],4),')')
cat('<H2>Neutrals</H2>')
cat(paste0(listSp,collapse='<BR>'))

# Cations
species = trueSpecies[selIons]
spMass  = trueMasses[selIons]
mord=order(spMass)
listSp = c()
for (i in seq_along(mord)) 
  listSp[i] = paste0('<font color="blue">',species[mord[i]],
                     '</font> (',signif(spMass[mord[i]],4),')')
cat('<H2>Cations</H2>')
cat(paste0(listSp,collapse='<BR>'))

sink(file = NULL)

setwd(sourceDir)
listHtml=paste0('"',webDir,listReacs,'/summary.html"')
listHtml=paste(listHtml,collapse=' ')
listHtml=paste0('Doc/ReleaseNotes.html ','speciesList.html ',listHtml)
command=paste0('htmldoc --book --toclevels 1 --size A4 ',
               '--compression=5 --fontsize 10 --linkcolor purple ',
               '-f summary.pdf ',listHtml, collapse=' ')
system(command,intern=FALSE)
file.rename(from='summary.pdf',to=paste0(publicDir,'summary.pdf'))

# Full biblio
setwd(docDir)
bibFile='refsDR.bib'
if(!file.exists('bib.Rdata')) {
  cat('*** Processing .bib file\n')
  bib = read.bib(file=bibFile)
  save(bib, file='bib.Rdata')
} else {
  sourceTime = file.info(bibFile)["mtime"]
  bibTime    = file.info('bib.Rdata')["mtime"]
  if(sourceTime > bibTime) {
    cat('*** Processing .bib file\n')
    bib = read.bib(file=bibFile)
    save(bib, file='bib.Rdata')    
  } else {
    cat('*** Loading  processed .bib file\n')
    load('bib.Rdata')
  }
}
targetHtml = paste0(tmpDir,'bibliography.html')
sink(file=targetHtml, append=FALSE)
allBibKeys=sort(unique(allBibKeys))
printBib(allBibKeys,bib)
sink(file = NULL)
publicPDF = paste0(publicDir,'bibliography.pdf')
command=paste0('pandoc -V geometry:margin=2cm ',targetHtml,' -o ',publicPDF)
system(command,intern=FALSE)



