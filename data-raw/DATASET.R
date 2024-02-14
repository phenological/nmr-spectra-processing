#Generator of test data
#avo3: first three methanol samples in avocado NMR dataset from ANPC, 2022
#ppm: corresponding ppm scale
load("avo3.rda")
avo3tsp <- avo3[,crop(ppm,-0.02,0.02)]
avo3Ala <- avo3[,crop(ppm,1.47,1.52)]
avo3StartToAla <- avo3[,crop(ppm,end=1.52)]
avo3SucrToEnd <- avo3[,crop(ppm,start=5.4)]
avo3CalAla <- rbind(shiftSeries(avo3[1,],-60,padding="zeroes")
                    ,shiftSeries(avo3[2,],-67,padding="zeroes")
                    ,shiftSeries(avo3[3,],-53,padding="zeroes"))
sucrose <- new("NMRSignal1D"
               ,id="Sucrobation"
               ,chemicalShift=5.416
               ,shape=list(name="pseudoVoigt"
                           ,params = list(mu=0.85,fwhm=0.0015,base=0)
               )
               ,peaks = list(new("NMRPeak1D",x=5.4113, y=1)
                             ,new("NMRPeak1D",x=5.4207, y=1)
               )
)

avo3CalSucro <- rbind(shiftSeries(avo3[1,],1,padding="zeroes")
                      ,avo3[2,]
                      ,shiftSeries(avo3[3,],2,padding="zeroes"))

y1 <- c(rep(0,7),avo3[1,1:40799])
y2 <- c(avo3[2,11:40806],rep(0,10))
avo3RandShift <- rbind(y1,y2,avo3[3,])
dimnames(avo3RandShift) <- NULL

usethis::use_data(ppm, avo3
                  ,avo3tsp, avo3Ala, avo3StartToAla, avo3SucrToEnd
                  ,avo3RandShift, avo3CalAla, avo3CalSucro
                  ,overwrite = TRUE, internal=TRUE)
