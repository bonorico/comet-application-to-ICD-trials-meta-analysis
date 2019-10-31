## comet_ICD_analysis.R contains R commands to execute 'comet application to ICD trials meta-analysis' (from homonymous repository). 
## Copyright (C) 2019 Federico Bonofiglio

    ## This Program is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## This Program is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with This Program.  If not, see <https://www.gnu.org/licenses/>. 



#### meta-analysis on implantable cardioverter defibrillator clinical trials
#### see: Bonofiglio, Federico, et al. "Meta‐analysis for aggregated survival data with competing risks: a parametric approach using cumulative incidence functions." Research synthesis methods 7.3 (2016): 282-293.
###################################################
### code chunk number 1: dataload
###################################################

if (!require("devtools")) {
    install.packages("devtools")
    library(devtools)
}

install_github("bonorico/comet")

library(comet)

library(xtable)


icd.export <- read.table("dataextraction_icd_final_14032007_export.csv",
                         sep="\t", header=TRUE)


#########################################
 

icd.export <- icd.export[-6,]  # remove Connolly ipd study (no ref.)

fup <- c(27,18,32,57,35,20,36,29,39)/12    # I use fup by koller paper

icd.export$meanfoup <- fup


## dataprep (new code)

data1 <- prep.metacr( event=list(TR1= data.frame( 
                                     E1=d_arr_icd , E2=d_nonarr_icd ), 
     TR2=data.frame(E1=d_arr_cont, E2=d_nonarr_cont) ),
                     
ptime= list(TR1= data.frame( E1=py_icd ,
         E2=py_icd ), 
     TR2=data.frame(E1= py_cont,
         E2= py_cont) ),
                     
 N=list(TR1= data.frame( E1=No..ICD , E2=No..ICD ), 
     TR2=data.frame(E1=No..cont, E2=No..cont) ), studlab=Trial, 
                     fup=meanfoup,                  
data=icd.export )


## alternative data with inflated hazards

cf <- 10  # inflating factor

data1b <- data.frame( ptime= data1$ptime / cf, data1[ , -2] )
     


###################################################
### code chunk number 2: calculations
###################################################

## COMPUTATIONS OF CIFs

args(metacinc)

cifs <- metacinc(event, N, 
                ptime, cause,
                     time, group, study, control="TR2",
                      data=data1)

cifsINFL <- update(cifs, data=data1b) # inflated CIFs

# plot CIFs

# plot prep
cifs.plot <- update(cifs, binn=100)
cifsINFL.plot <- update(cifsINFL, binn=100)


args(plot.metacif)

ell <- c("ICD","CONTROL")
ev <- c("SC death", "Other death")



postscript("FIG-2.ps", width=11)
plot(cifs.plot, x_lab="Time (Years)", ev_lab=ev, arm_lab=ell, 
     title=" ")

dev.off()





### META-ANALYSIS

meta <- update(cifs, predt=6, pool=TRUE)

meta.land <- update(meta, predt=NULL, landmark=TRUE)

# inflated MA
metaINFL <- update(cifsINFL,predt=6, pool=TRUE) 

metaINFL.land <- update(metaINFL, predt=NULL, landmark=TRUE)


## PLOT MA 
#plot prep 
meta.plot <- update(meta, plot=TRUE)
metaINFL.plot <- update(metaINFL, plot=TRUE)

meta.land.plot <- update(meta.land, plot=TRUE)

metaINFL.land.plot <- update(metaINFL.land, plot=TRUE)


#

x <- plot(meta.plot, x_lab="Time (Years)", ev_lab=ev, arm_lab=ell, 
     title="A", fixed=FALSE,y_extra=0.05)$plot 


y <- plot(meta.land.plot, x_lab="Time (Years)", ev_lab=ev, arm_lab=ell, 
     title="B ", text=TRUE, fixed=FALSE, y_extra=0.5)$plot  




 g <- arrangeGrob(x$out, y$out)

 ggsave("FIG-3.ps", g, device=cairo_ps, dpi=800,  
      width = 10.0, height = 15.0, units="in")



  #inflation



z <- plot(metaINFL.plot, x_lab="Time (Years)", ev_lab=ev, arm_lab=ell, 
     title="A ", ci=FALSE, fixed=FALSE, y_extra=0.05 )$plot


r <- plot(metaINFL.land.plot, x_lab="Time (Years)", ev_lab=ev, arm_lab=ell, 
     title=" B",text=TRUE, ci=FALSE, fixed=FALSE, y_extra=0.5 )$plot




postscript("FIG-4.ps", width = 10.0, height = 15,
                horizontal = FALSE, 
                onefile = FALSE, paper = "special")
grid.arrange(z$out, r$out, nrow=2)

dev.off()


# use pdf(file.pdf)  or postscript(file.eps) to print


### print RR results in percent change 

 prr <- function( model, type, cause,func="max", data){  # 

     
 data <- subset(data, MODEL==model & TYPE==type 
             & CAUSE==cause, 
                select=c(LOGPOOL,LOGLOW,LOGUP,TIME))

 data[, 1:3] <- round((exp(data[,1:3])-1)*100,0)
    

  time <- with(data,# max(TIME, na.rm=TRUE)
               switch(func, mean=mean(TIME, na.rm=TRUE),
                            max=max(TIME, na.rm=TRUE),
                             min=min(TIME, na.rm=TRUE) ) 
    )
  
 data <- subset(data, TIME==time)
 

     
pci <- function(x){
      if(x>0)
             noquote(paste("+",x, sep=""))
                   else
            noquote( paste(x, sep="") )
               }

     x1 <- pci(data[,1])
      x2 <- pci(data[,2])
    x3 <- pci(data[,3])
     
noquote(
    paste( x1, " (", x2, "; ", x3, ") ", sep="" )
    )
     
 
  }


#prr( "RANDOM", "RR", "E1", data=meta.plot$MA)





###################################################
### code chunk number 3: tableone
###################################################

trial <- matrix(noquote(
paste(icd.export$Trial," (",icd.export$Year,")", sep="")

),ncol=1)   # first column


tab1 <- with(icd.export,
     
     
     data.frame("Trial(year)"=trial,"Follow-up(years)"=fup, "SC deathsICD"= d_arr_icd, "other deathsICD"= d_nonarr_icd,
                
                "Person-timeICD"= py_icd, "SC deathsCONT"= d_arr_cont, "other deathsCONT"= d_nonarr_cont, 
                
                "Person-timeCONT"= py_cont
                                
                
                )
    
     
     )



colnames(tab1) <- c("","","","","","","","")




###################################################
### code chunk number 4: printtableone
###################################################


print.xtable(
             xtable(tab1, digits=c(0, 0, 2, 0, 0, 1, 0, 0, 1),
                    align=c("l","l","c","r","r","r","r","r","r"), # first term is rownames
                    label="tab1", caption="Data on nine randomized clinical trials collected by \\cite{Koll:Stij:Stey:Lub:2008}.
Experimental = implantable cardioverter defibrillator (ICD); SC = sudden cardiac; Follow-up (mean) and Person-time are in years." ), floating.environment="sidewaystable",
             include.rownames=F,include.colnames=F, booktabs=T,
             , size=5, caption.placement="top", 
             add.to.row= list(pos=list(-1,0),command=c("\\toprule & &  \\multicolumn{3}{c}{\\bfseries \\itshape Experimental}& 
\\multicolumn{3}{c}{\\bfseries \\itshape Control} \\\\ [0.2cm] % ",
                                      
                                               "\\itshape Trial (Year) & \\itshape Follow-up & \\itshape SC deaths & \\itshape Other deaths & \\itshape Person-time & \\itshape SC deaths & \\itshape Other deaths & \\itshape Person-time \\\\"))
             )



