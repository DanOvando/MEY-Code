rm(list=ls())
tic <- proc.time()
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
#RAMdata = read.csv('../Data/RAM Stock Parameters.csv',stringsAsFactors=F)


# Prep Data ---------------------------------------------------------------

load('Data/Policies and Params for Chris.Rdata')

DANdata = for_chris

rm(for_chris)
# DANdata = read.csv('Data/Policies and Params for Chris.csv',stringsAsFactors=F)

# ggKobe(DANdata)

Fisheries = unique(DANdata$IdOrig)
MEYdata = data.frame(matrix(NA,nrow=length(Fisheries),ncol=7))
current_b = vector()
current_f = vector()
current_b_mey = vector()
current_f_mey = vector()
b_mey = vector()
f_mey = vector()
grow_disc = vector()

# Dynamic Optim Functions -------------------------------------------------

DiffF=function(b,bgrid,f0,f1)
{
  #To be zeroed out.  Want to find b* such that f0(b*) = f1(b*)
  fval0 = spline(bgrid,f0,xout=b,method="natural")
  fval1 = spline(bgrid,f1,xout=b,method="natural")
  difference = fval0$y - fval1$y
  return(difference)
}


# Run Dynamic Optim -------------------------------------------------------


for (i in seq(1,length(Fisheries),1))
{
  print(i)
  DANsub=subset(DANdata,DANdata$IdOrig==Fisheries[i])
  bvec = DANsub$b
  fvec = DANsub$Opt

  phi = DANsub$phi[1]
  g = DANsub$g[1]

  if (is.na(phi))
  {phi=.188}

  fpt = ((phi+1)/phi)*(1 - bvec^phi/(phi+1))

  bmey = fzero(DiffF,1.2,bgrid=bvec,f0=fvec,f1=fpt)
  fmey = ((phi+1)/phi)*(1 - bmey$x^phi/(phi+1))

  current_b[i] = DANsub$BvBmsy[1] #current B/Bmsy
  current_f[i] = DANsub$FvFmsy[1] #current F/Fmsy
  b_mey[i] = bmey$x #Bmey/Bmsy
  f_mey[i] = fmey #Fmey/Fmsy
  current_b_mey[i] = DANsub$BvBmsy[1]/bmey$x #B/Bmey
  current_f_mey[i] = DANsub$FvFmsy[1]/fmey #F/Fmey
  grow_disc[i] = ((phi+1)/phi)*g
}
proc.time()  - tic

MEYdata = data.frame(IdOrig = Fisheries,current_b,current_f,
                     b_mey,f_mey,current_b_mey, current_f_mey)
write.csv(MEYdata,file='MEY_results.csv')

kobe_dat <- DANdata %>%
  select(IdOrig, Dbase, SciName, CommName, MSY, SpeciesCat, SpeciesCatName, BvBmsy, FvFmsy) %>%
  unique() %>%
  left_join(MEYdata,by = c('IdOrig'))

head(kobe_dat)

# Make Kobe Plot ----------------------------------------------------------
source("MEY_Functions/ggKobe.R")

kobe_mey <- ggKobe(kobe_dat, xvar = 'current_b_mey', yvar = 'current_f_mey' ) +
  labs(x = 'B/Bmey', y = 'F/Fmey')

ggsave('MEY Kobe.pdf', kobe_mey)

kobe_msy <- ggKobe(kobe_dat, xvar = 'current_b', yvar = 'current_f' ) +
  labs(x = 'B/Bmsy', y = 'F/Fmsy')

ggsave('MSY Kobe.pdf', kobe_msy)

kobes <- arrangeGrob(kobe_msy, kobe_mey)

grid.draw(kobes)

ggsave('Kobe Comparison.pdf', kobes)


# Match to NEIs



