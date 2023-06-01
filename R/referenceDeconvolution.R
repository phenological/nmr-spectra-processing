#Candidate to utilities.R
crop <- function(x,roi){
  x >= roi[1] & x <= roi[2]
}

#' Reference deconvolution using TSP peak as reference
#' @importFrom gsignal hilbert
#' @export
#TBD: handle references other than TSP
#TBD: improved frequency-based criterion to determine the inner extremes of the
#region where satellites' maxxima are sought (higher frequency should mean
#that we can move further because there is less overlap with the main peak)
#For that reason the current threshold is expected to give problems, specially
#since it was calculated for 600 Hz but the default frequency of the function
#is 400 Hz...
referenceDeconvolution <- function(x,y,rOref=c(-0.02,0.02),signals,frequency=400
                                   ,padding="sampling",using=c(9.5,10),from=y
                                   ,optimizationOptions=list(
                                     optimization=list(kind="lm")
                                     ,shape=list(kind="lorentzian"
                                                 ))
                                   ,...
                                   ){
  #inverse fourier transform with normalization
  ifft <- function(x){
  fft(x, inverse = TRUE) / length(x)
  }
  
  #Due to overlap with the main peak, special measures are required to model
  #the satellites. This function takes care of that.
  #TBD: handle the case where the symetrization elongates the signal beyond the
  #roi, or even beyond the range of x (is that even possible on a sensible roi?)
  #TBD: handle the case where neither end of the signal's delta reaches an extreme
  #of the roi
  optimizeOverlappedSignal <- function(x,y,signal,...){
    l_x <- length(x)
    x_signal <- x[crop(x,signal$delta)]
    #print(x_signal)
    y_signal <- y[crop(x,signal$delta)]
    #print(y_signal)
    center <- which.max(y_signal)
    #print(center)
    signal$delta <- x_signal[center]
    l_s <- length(y_signal)
    if (y_signal[1] > y_signal[length(x_signal)]){
      #plot(y_signal)
      y_signal <- c(y_signal[l_s:(center-1)]
                        ,y_signal[center:l_s])
      #plot(y_signal)
      trail_l <- y[1:(l_x - length(y_signal))] #l_x - 2*(l_s - center) - 1
      y_signal <- c(trail_l,y_signal)
      #plot(y_signal)
    }
    else{
      y_signal <- c(y_signal[1:(center-1)],y_signal[center:1])
      trail_r <- y[(length(y_signal) + 1):l_x]
      y_signal <- c(y_signal,trail_r)
     # plot(x,y_signal,type="l")
    }
    
    optimizeSignals(x,y_signal,signals=list(signal),...)
  }
  
  #parse using as chem. shift. Any other form will be left for pad to parse.
  if (is.numeric(using) & length(using) == 2)
    using <- crop(x,using)
  
  #qc input lengths
  if(length(x)!=length(y)){
    cat(crayon::red("referenceDeconvolution >>", "lengths of x and y don't match\n"))
    stop()
  }
  
  #Build default TSP model, if necessary
  if (missing(signals)){
    #Compute the inner extreme of the region where satellites' maxima are sought
    satellites <- 2.7 / frequency
    signals <- list(clean=list(list(delta=0))
                    ,overlap=list(
                      list(delta=c(rOref[1],-satellites))
                      ,list(delta=c(satellites,rOref[2]))
                    ))
  }
  
  #avoid zeros in y
  #y[y==0] <- zero
   
  #parse and construct rOref
  if (is.numeric(rOref) & length(rOref==2)){
    rOref[1] <- max(rOref[1],min(x))
    rOref[2] <- min(rOref[2],max(x))
    rOref <- x >= rOref[1] & x <= rOref[2]
  }
  else{
    cat(crayon::red("referenceDeconvolution >>", "invalid rOref\n"))
    stop()
  }
  
  #calculate trailing zeros outside of the rOref
  #trail_l <- rep(0,which.max(rOref) - 1)
  #trail_r <- rep(0,which.max(rev(rOref)) - 1)
  
  #get reference, pad it
  ref <- y[rOref]
  #plot(ref,type="l")
  ref <- pad(ref,n=which.max(rOref)-1,side=-1,method=padding
             ,using=using,from=from)
  #plot(ref,type="l")
  ref <- pad(ref,n=which.max(rev(rOref))-1,side=1,method=padding
             ,using=using,from=from)
  #plot(x,ref,type="l")
  #ref <- c(trail_l, y[rOref], trail_r)
  
  #temporary workaround: if there are ppm<0 in the rOref, shift the spectrum 
  #to start at zero, fit the signal and then shift back
  hasNegative <- sum(x<0) > 0
  if (hasNegative){
    shift <- x[1]
    x <- x - shift
    for (i in 1:length(signals)){
      for (j in 1:length(signals[[i]])){
        signals[[i]][[j]]$delta <- signals[[i]][[j]]$delta - shift
      }
    }
  }
  #print(signals)
  #fit model
  main <- optimizeSignals(x[rOref],y[rOref],signals=signals$clean
                               ,frequency=frequency
                               ,options=optimizationOptions)
  #print(main$signals)
  satellite1 <- optimizeOverlappedSignal(x[rOref],y[rOref],signals$overlap[[1]]
                                         ,frequency=frequency
                                         ,options=optimizationOptions)
  #print(satellite1$signals)
  satellite2 <- optimizeOverlappedSignal(x[rOref],y[rOref],signals$overlap[[2]]
                                         ,frequency=frequency
                                         ,options=optimizationOptions)
  #print(satellite2$signals)
  if (is.null(main) | is.null(satellite1) | is.null(satellite2)){
      cat(crayon::yellow("referenceDeconvolution>>"
                         ,"signal fitting with fitSignals in the rOref returned NULL\n"
                         ,"deconvolution aborted\n"
                         ,"returning the input spectrum\n"))
      return(y)
  }
  #refModel <- rbind(main$signals,satellite1$signals,satellite2$signals)
  #print(refModel)
  refModel <- rbind(signalsToY(x,main$signals,frequency=frequency)
                  ,signalsToY(x,satellite1$signals,frequency=frequency)
                  ,signalsToY(x,satellite2$signals,frequency=frequency)
                             )
  refModel <- colSums(refModel)
  
  #shift the scale back to its original position, if necessary
  #Warning: only doing it for diagnostics; not needed once stable
  if (hasNegative) x <- x + shift
  
  plot(x[crop(x,c(0.015,0.03))],ref[crop(x,c(0.015,0.03))],type="l")
  plot(x[rOref],y[rOref],type="l")
  lines(x[rOref],refModel[rOref],col="red")
  
  #refModel <- signalsToY(x,refModel,frequency=frequency)
  #lines(x[rOref],refModel[rOref],col="red")
  #refModel <- fitSignals(x,y,signals=signals,roi=rOref,frequency=frequency,...)
  # if (is.null(refModel)){
  #   cat(crayon::yellow("referenceDeconvolution>>"
  #                      ,"signal fitting with fitSignals in the rOref returned NULL\n"
  #                      ,"deconvolution aborted\n"
  #                      ,"returning the input spectrum\n"))
  #   return(y)
  # }
  
  #refModel <- signalsToY(x,refModel$signals,frequency=frequency)
  
  #move to time domain
  ref <- ifft(hilbert(ref))
  refModel <- ifft(hilbert(refModel))
  y <- ifft(hilbert(y))
  
  #deconv., correct and go back to freq. domain
  res <- Re(fft(y * refModel / ref))
  lines(x[rOref],res[rOref],type="l",col="blue")
  res
}
