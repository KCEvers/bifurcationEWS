
#' Recurrence Quantification Analysis
#'
#' @param x Data, can be a vector or matrix or dataframe
#' @param emDim Embedding dimension
#' @param emLag Embedding lag
#' @param theiler Size of Theiler window
#' @param distNorm Distance norm used for constructing recurrence matrix
#' @param targetValue Target value for thresholding recurrence matrix
#' @param rescaleDist Rescale distance matrix
#' @param fix_emRad_or_RR Set radius or recurrence rate as fixed
#' @param ... Optional arguments
#'
#' @return RQA measures
#' @export
#'
#' @examples
runRQA <- function(x, emDim = 1, emLag = 1, theiler = 1, distNorm = "max", targetValue = .05,
                   rescaleDist = "none",
                   fix_emRad_or_RR = c( "emRad-space-diameter", "emRad-constant", "RR")[2],
                   ...){
  RM <- getRM(x,  emDim = emDim,
              emLag= emLag,
              theiler = theiler,
              distNorm = distNorm,
              rescaleDist = rescaleDist,
              targetValue = targetValue,
              fix_emRad_or_RR = fix_emRad_or_RR, ...)
  return(rp_measures_KCE(RM, theiler = theiler, targetValue = targetValue, ...))
}


#' Thresholded Recurrence Matrix
#'
#' @inheritParams runRQA
#'
#' @return Recurrence Matrix
#' @export
#'
#' @examples
getRM <- function(x,
                    emDim,
                    emLag,
                    theiler,
                    distNorm,
                    rescaleDist,
                    targetValue,
                    fix_emRad_or_RR = c( "emRad-space-diameter", "emRad-constant", "RR")[2] ){
  # # If emDim and/or emLag is NA, estimate
  # if (is.na(pars$emLag) & is.na(pars$emDim)) {
  #   optimPar <- casnet::est_parameters(
  #     x,
  #     lagMethods = pars$methodLag,
  #     #c("first.minimum","global.minimum","max.lag")[1],
  #     estimateDimensions = "preferSmallestInLargestHood",
  #     # maxDim   = pars$max.emb.dim,
  #     emLag    = NULL,
  #     #  maxLag   = pars$lag.max,
  #     minVecLength = 20,
  #     nnSize  = NA,
  #     nnRadius = NA,
  #     nnThres  = 10,
  #     theiler  = pars$theiler,
  #     doPlot   = TRUE,
  #     silent   = TRUE
  #   )
  #   pars$emLag <- optimPar$optimLag
  #   pars$emDim <- optimPar$optimDim
  #
  #   print(sprintf("emLag: %d", pars$emLag))
  #   print(sprintf("emDim: %d", pars$emDim))
  # }
#
  # casnet::rp(
  #   x,
  #   emDim   = emDim,
  #   emLag   = emLag,
  #   emRad   = NA, #pars$emRad,
  #   theiler = theiler,
  #   method = distNorm,
  #   targetValue    = targetValue
  # )

  if (fix_emRad_or_RR == "emRad-space-diameter"){

    RM_unthresh <- casnet::rp(
      x,
      emDim   = emDim,
      emLag   = emLag,
      theiler = theiler,
      method = distNorm,
      rescaleDist = rescaleDist
    ) %>% as.matrix()

    # To define a radius in terms of the maximum phase space diameter, look at unthresholded matrix and exclude the diagonal and theiler window
    RM_unthresh_incl = RM_unthresh[col( RM_unthresh) - row(RM_unthresh) > theiler]
    emRad = targetValue * diff(range(RM_unthresh_incl)) + min(RM_unthresh_incl)

    RM <- casnet::rp(
      x,
      emDim   = emDim,
      emLag   = emLag,
      emRad   = emRad, # Percentage of maximum phase space diameter
      theiler = theiler,
      method = distNorm,
      rescaleDist = rescaleDist
    )
  } else if (fix_emRad_or_RR == "emRad-constant"){

    RM <- casnet::rp(
      x,
      emDim   = emDim,
      emLag   = emLag,
      emRad   = targetValue, # Percentage of maximum phase space diameter
      theiler = theiler,
      method = distNorm,
      rescaleDist = rescaleDist
    )
  } else if (fix_emRad_or_RR == "RR"){
    RM <- casnet::rp(
      x,
      emDim   = emDim,
      emLag   = emLag,
      emRad   = NA,
      theiler = theiler,
      method = distNorm,
      targetValue    = targetValue
    )
  }

  # # Binarize recurrence matrix
  # RM <- casnet::mat_di2bi(RM,emRad)

  return(RM)
}


#' Get RQA measures based on a binary matrix
#'
#' A zoo of measures based on singular recurrent points, diagonal, vertical, horizontal, and white line structures (anisotropic) will be caluclated.
#'
#' @param RM A binarized distance matrix
#' @param DLmin Minimal diagonal line length (default = `2`)
#' @param VLmin Minimal vertical line length (default = `2`)
#' @param HLmin Minimal horizontal line length (default = `2`)
#' @param WLmin Minimal white line length (default = `1`) # KCE
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param WLmax Maximal white line length (default = length of diagonal -1) # KCE
#' @param includeDiagonal Should the diagonal be included in the calculation of the Recurrence Rate? If `NA` the value will be decided by the symmetry of the matrix, the diagonal will be removed for Auto RQA (`AUTO = TRUE`) but not for Cross RQA (`AUTO = FALSE`)  (default = `NA`)
#' @param matrices Return matrices? (default = `FALSE`)
#' @param targetValue A value passed to `est_radius(...,type="fixed", targetMeasure="RR", tol = .2)` if `is.na(emRad)==TRUE`, it will estimate a radius (default = `.05`).
#' @param theiler Theiler window
#' @param AUTO Auto-recurrence?
#' @param silent Do not display any messages (default = `TRUE`)'
rp_measures_KCE <- function(RM,
                            DLmin = 2,
                            VLmin = 2,
                            HLmin = 2,
                            WLmin = 1, # KCE
                            DLmax = length(Matrix::diag(RM)),
                            VLmax = length(Matrix::diag(RM)),
                            HLmax = length(Matrix::diag(RM)),
                            WLmax = length(Matrix::diag(RM)), # KCE
                            AUTO      = NULL,
                            theiler   = NA,
                            includeDiagonal = NA,
                            # recurrenceTimes = FALSE,
                            matrices  = FALSE,
                            targetValue = .05,
                            silent     = TRUE){

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM <- casnet::rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE, checkTSPARSE = TRUE, fixTSPARSE = TRUE)

  if(is.null(AUTO)){
    AUTO <- attr(RM,"AUTO")
  }

  # if(is.na(emRad)){
  #   if(!is.null(attributes(RM)$emRad)){
  #     emRad <- attributes(RM)$emRad
  #   } else {
  #     # Check for attributes
  #     if(is.na(targetValue)){
  #       targetValue <- .05
  #     }
  #     if(!is.null(attributes(RM)$emDim)){
  #       emDim <- attributes(RM)$emDim
  #     } else {
  #       emdim <- 1
  #     }
  #     if(!is.null(attributes(RM)$emLag)){
  #       emLag <- attributes(RM)$emLag
  #     } else {
  #       emLag <- 1
  #     }
  #
  #     emRad <- est_radius(RM = RM, emDim = emDim, emLag = emLag, targetValue = targetValue, tol = .2, radiusOnFail = "percentile", silent = silent)
  #     if(emRad$Converged){
  #       emRad <- emRad$Radius
  #     } else {
  #       emRad <- stats::sd(RM,na.rm = TRUE)
  #     }
  #   }
  # }

  if (all(as.vector(RM)==0)){
    warning("Nothing is recurring!")
    return(RM_measures_ = data.frame(
                                     RP_max = Matrix::nnzero(RM, na.counted = FALSE),
                                     RR = 0,
                                     SING_N = 0,
                                     SING_rate = 0,
                                     RP_N = Matrix::nnzero(RM, na.counted = FALSE),
                                     DIV_dl = 0,
                                     REP_av = 0,
                                     N_dl = 0,
                                     N_dlp = 0,
                                     DET = 0,
                                     MEAN_dl = 0,
                                     MAX_dl = 0,
                                     ENT_dl = 0,
                                     ENTrel_dl = 0,
                                     CoV_dl = 0,
                                     N_vl = 0,
                                     N_vlp = 0,
                                     LAM_vl = 1,
                                     TT_vl = 0,
                                     MAX_vl = 0,
                                     ENT_vl = 0,
                                     ENTrel_vl = 0,
                                     CoV_vl = 0,
                                     REP_vl = 1,
                                     N_hlp = 0,
                                     N_hl = 0,
                                     LAM_hl = 1,
                                     TT_hl = 0,
                                     MAX_hl = 0,
                                     ENT_hl = 0,
                                     ENTrel_hl = 0,
                                     CoV_hl = 0,
                                     REP_hl = NA,
                                     N_hvp = NA,
                                     N_hv = NA,
                                     LAM_hv = NA,
                                     TT_hv = NA,
                                     MAX_hv = NA,
                                     ENT_hv = 0,
                                     ENTrel_hv = 0,
                                     CoV_hv = 0,
                                     N_wl= 0,
                                     N_wlp = 0,
                                     LAM_wl = 0,
                                     ENT_wl = 0,
                                     ENTrel_wl = 0,
                                     MEAN_wl = 0,
                                     MAX_wl = 0,
                                     CoV_wl = 0,
                                     DIV_wl = 0,
                                     LAM_DET = 0,
                                     DET_RR = 0,
                                     DIV_vl = 0))
  } else if (all(as.vector(RM)==1)){
    warning("Everything is recurring!")
    return(data.frame(
                      RP_max = Matrix::nnzero(RM, na.counted = FALSE),
                      RR = 1,
                      SING_N = 0,
                      SING_rate = 0,
                      RP_N = Matrix::nnzero(RM, na.counted = FALSE),
                      DIV_dl = 1 / (dim(RM)[1] - theiler),
                      REP_av = 1,
                      N_dl = (dim(RM)[1] - theiler) * 2,
                      N_dlp = Matrix::nnzero(RM, na.counted = FALSE),
                      DET = 1,
                      MEAN_dl = mean(seq(DLmin, dim(RM)[1] - theiler)),
                      MAX_dl = dim(RM)[1] - theiler,
                      ENT_dl = 0,
                      ENTrel_dl = 0,
                      CoV_dl = 0,
                      N_vl = (dim(RM)[1] - theiler) * 2,
                      N_vlp = Matrix::nnzero(RM, na.counted = FALSE),
                      LAM_vl = 1,
                      TT_vl = mean(seq(VLmin, dim(RM)[1] - theiler)),
                      MAX_vl =  dim(RM)[1] - theiler,
                      ENT_vl = 0,
                      ENTrel_vl = 0,
                      CoV_vl = 0,
                      REP_vl = 1,
                      N_hlp = Matrix::nnzero(RM, na.counted = FALSE),
                      N_hl = (dim(RM)[1] - theiler) * 2,
                      LAM_hl = 1,
                      TT_hl = mean(seq(VLmin, dim(RM)[1] - theiler)),
                      MAX_hl =  dim(RM)[1] - theiler,
                      ENT_hl = 0,
                      ENTrel_hl = 0,
                      CoV_hl = 0,
                      REP_hl = NA,
                      N_hvp = NA,
                      N_hv = NA,
                      LAM_hv = NA,
                      TT_hv = NA,
                      MAX_hv = NA,
                      ENT_hv = 0,
                      ENTrel_hv = 0,
                      CoV_hv = 0,
                      N_wl= 0,
                      N_wlp = 0,
                      LAM_wl = 0,
                      ENT_wl = 0,
                      ENTrel_wl = 0,
                      MEAN_wl = 0,
                      MAX_wl = 0,
                      CoV_wl = 0,
                      DIV_wl = 0,
                      LAM_DET = 1,
                      DET_RR = 1,
                      DIV_vl = 1 / dim(RM)[1] - theiler))
  }
  # Compute RQA measures from binarized recurrence matrix
  suppressMessages(out <- rp_calc_KCE(RM,
                                      DLmin = DLmin,
                                      VLmin = VLmin,
                                      HLmin = HLmin,
                                      WLmin = WLmin, # KCE
                                      DLmax = DLmax,
                                      VLmax = VLmax,
                                      HLmax = HLmax,
                                      WLmax = WLmax, # KCE
                                      theiler = theiler,
                                      AUTO  = AUTO,
                                      includeDiagonal = includeDiagonal,
                                      # recurrenceTimes = recurrenceTimes,
                                      matrices  = matrices))

  if(matrices){
    tab <- out$crqaMeasures
  } else {
    tab <- out
  }
  return(invisible(out))
}



# KCE: Changes to rp_calc
# - Added arguments WLmin, WLmax
# - Added WLmin, WLmax, recurrenceTimes to rp_calc_lineMeasures calls
#


#' rp_calc
#'
#' @inheritParams rp_measures_KCE
#'
#' @return CRQA measures and matrices of line distributions (if requested)
#' @export
#'
rp_calc_KCE <- function(RM,
                        DLmin = 2,
                        VLmin = 2,
                        HLmin = 2,
                        WLmin = 1, # KCE
                        DLmax = length(Matrix::diag(RM)),
                        VLmax = length(Matrix::diag(RM)),
                        HLmax = length(Matrix::diag(RM)),
                        WLmax = length(Matrix::diag(RM)), # KCE
                        theiler   = NA,
                        AUTO      = NULL,
                        includeDiagonal = NA,
                        # recurrenceTimes = FALSE,
                        matrices  = FALSE){

  RM <- casnet::rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE)

  if(is.null(AUTO)){
    AUTO <- attributes(RM)$AUTO
  }

  if(is.na(theiler)){
    if(is.na(attributes(RM)$theiler)){
      RM <- casnet::setTheiler(RM, theiler = theiler)
    }
    theiler <- attributes(RM)$theiler
  }

  if(!is.na(includeDiagonal)){
    AUTO <- includeDiagonal
  }

  lineMeasures_global <- rp_calc_lineMeasures_KCE(RM = RM,
                                                  DLmin = DLmin, DLmax = DLmax,
                                                  VLmin = VLmin, VLmax = VLmax,
                                                  HLmin = HLmin, HLmax = HLmax,
                                                  WLmin = WLmin, WLmax = WLmax, # KCE
                                                  # recurrenceTimes = recurrenceTimes, # KCE
                                                  theiler = theiler,
                                                  AUTO = AUTO,
                                                  matrices = matrices)

  if(matrices){
    lineMatrices_global <- lineMeasures_global$crqaMatrices
    lineMeasures_global <- lineMeasures_global$crqaMeasures
  }

  # # Singularities
  # SING <- rp_lineDist_KCE(RM,
  #                         DLmin = 1, DLmax = DLmax,
  #                         VLmin = 1, VLmax = VLmax,
  #                         HLmin = 1, HLmax = HLmax,
  #                         theiler = theiler, AUTO = AUTO)
  #
  # # table(SING_N$verticals.dist)
  # # table(SING_N$horizontals.dist)
  # SING_N <- table(SING$diagonals.dist)[1]
  # SING_rate <- SING_N / RP_N

  #Output
  out <- as.data.frame(list(lineMeasures_global))

  if(matrices){
    return(list(crqaMeasures = out,
                crqaMatrices = lineMatrices_global)
    )
  } else {
    return(out)
  }
}


# KCE:
# Changes to rp_calc_lineMeasures()
# - Added argument to function WLmin = 2
# - Added argument to function WLmax = length(Matrix::diag(RM))-1
# - Added white lines to most calculations
# - Added argument recurrenceTimes

#' Compute RQA measures
#'
#' @inheritParams rp_measures_KCE
#' @param d Vector of diagonals to be extracted from matrix `RP` before line length distributions are calculated. A one element vector will be interpreted as a windowsize, e.g., `d = 50` will extract the diagonal band `-50:50`. A two element vector will be interpreted as a band, e.g. `d = c(-50,100)` will extract diagonals `-50:100`. If `length(d) > 2`, the numbers will be interpreted to refer to individual diagonals, `d = c(-50,50,100)` will extract diagonals `-50,50,100`. If `length(d)` is `NULL`, 1 or 2, the theiler window is applied before diagonals are extracted. The theiler window is ignored if `length(d)>2`, or if it is larger than the matrix or band indicated by parameter `d`. A warning will be given is a theiler window was already applied to the matrix.
#'
#' @return Dataframe with RQA line measures
#' @importFrom invctr  `%00%`
#' @export
#'
#' @examples
rp_calc_lineMeasures_KCE <- function(RM,
                                     DLmin = 2,
                                     VLmin = 2,
                                     HLmin = 2,
                                     WLmin = 1, # KCE
                                     DLmax = length(Matrix::diag(RM)),
                                     VLmax = length(Matrix::diag(RM)),
                                     HLmax = length(Matrix::diag(RM)),
                                     WLmax = length(Matrix::diag(RM)), # KCE
                                     d         = NULL,
                                     theiler   = NA,
                                     # invert    = FALSE,
                                     # recurrenceTimes = FALSE, # KCE
                                     AUTO      = TRUE,
                                     matrices  = FALSE){

  # Total nr. recurrent points
  RP_N <- Matrix::nnzero(RM, na.counted = FALSE)
  recmatsize <- casnet::rp_size(RM, AUTO = AUTO, theiler = theiler)

  #Proportion recurrence / Recurrence Rate
  RR <- RP_N/recmatsize$rp_size_theiler

  if(length(RR)==0){RR<-0}

  if(RR==1){
    warning("Everything is recurring!\nReturning empty vector")
    return(casnet::rp_empty())
  }

  RP_max     = recmatsize$rp_size_theiler

  #Get line segments
  # if(Matrix::nnzero(RM)>0)
  lineSegments <- rp_lineDist_KCE(RM,
                                  DLmin = DLmin, DLmax = DLmax,
                                  VLmin = VLmin, VLmax = VLmax,
                                  HLmin = HLmin, HLmax = HLmax,
                                  WLmin = WLmin, WLmax = WLmax, # KCE
                                  d = d, theiler = theiler,
                                  # invert = FALSE, # KCE
                                  # recurrenceTimes = TRUE, #KCE
                                  # AUTO = AUTO,
                                  matrices = matrices
                                  )

  dlines <- lineSegments$diagonals.dist%00%0
  vlines <- lineSegments$verticals.dist%00%0
  hlines <- lineSegments$horizontals.dist%00%0
  wlines <- lineSegments$whites.dist%00%0 # KCE

  #Frequency tables of line lengths
  freq_dl <- table(dlines)
  freq_vl <- table(vlines)
  freq_hl <- table(hlines)
  freq_hv <- table(c(hlines,vlines))
  freq_wl <- table(wlines) # KCE

  freqvec_dl <- as.numeric(names(freq_dl))
  freqvec_vl <- as.numeric(names(freq_vl))
  freqvec_hl <- as.numeric(names(freq_hl))
  freqvec_hv <- as.numeric(names(freq_hv))
  freqvec_wl <- as.numeric(names(freq_wl)) # KCE


  # Number of lines
  N_dl <- sum(freq_dl, na.rm = TRUE)%00%0
  N_vl <- sum(freq_vl, na.rm = TRUE)%00%0
  N_hl <- sum(freq_hl, na.rm = TRUE)%00%0
  N_hv <- sum(freq_hv, na.rm = TRUE)%00%0
  N_wl <- sum(freq_wl, na.rm = TRUE)%00%0 # KCE

  #Number of recurrent points on diagonal, vertical and horizontal lines
  N_dlp <- sum(freqvec_dl*freq_dl, na.rm = TRUE)
  N_vlp <- sum(freqvec_vl*freq_vl, na.rm = TRUE)
  N_hlp <- sum(freqvec_hl*freq_hl, na.rm = TRUE)
  N_hvp <- sum(freqvec_hv*freq_hv, na.rm = TRUE)
  N_wlp <- sum(freqvec_wl*freq_wl, na.rm = TRUE) # KCE

  #Determinism / Horizontal and Vertical Laminarity
  DET    <- N_dlp/RP_N
  LAM_vl <- N_vlp/RP_N
  LAM_hl <- N_hlp/RP_N
  LAM_hv <- N_hvp/(RP_N*2)
  LAM_wl <- N_wlp/RP_N

  #Array of probabilities that a certain line length will occur (all >1)
  P_dl <- freq_dl/N_dl
  P_vl <- freq_vl/N_vl
  P_hl <- freq_hl/N_hl
  P_hv <- freq_hv/N_hv
  P_wl <- freq_wl/N_wl # KCE

  #Entropy of line length distributions
  ENT_dl <- -1 * sum(P_dl * log(P_dl))
  ENT_vl <- -1 * sum(P_vl * log(P_vl))
  ENT_hl <- -1 * sum(P_hl * log(P_hl))
  ENT_hv <- -1 * sum(P_hv * log(P_hv))
  ENT_wl <- -1 * sum(P_wl * log(P_wl)) # KCE

  #Relative Entropy (Entropy / Max entropy)
  ENTrel_dl = ENT_dl/(-1 * log(1/DLmax))
  ENTrel_vl = ENT_vl/(-1 * log(1/VLmax))
  ENTrel_hl = ENT_hl/(-1 * log(1/HLmax))
  ENTrel_hv = ENT_hv/(-1 * log(1/max(c(HLmax,VLmax), na.rm = TRUE)))
  ENTrel_wl = ENT_wl/(-1 * log(1/WLmax)) # KCE

  #Meanline
  MEAN_dl = mean(dlines, na.rm = TRUE)%00%0
  MEAN_vl = mean(vlines, na.rm = TRUE)%00%0
  MEAN_hl = mean(hlines, na.rm = TRUE)%00%0
  MEAN_hv = mean(c(hlines,vlines), na.rm = TRUE)%00%0
  MEAN_wl = mean(wlines, na.rm = TRUE)%00%0 # KCE

  #Maxline
  MAX_dl = max(freqvec_dl, na.rm = TRUE)%00%0
  MAX_vl = max(freqvec_vl, na.rm = TRUE)%00%0
  MAX_hl = max(freqvec_hl, na.rm = TRUE)%00%0
  MAX_hv = max(freqvec_hv, na.rm = TRUE)%00%0
  MAX_wl = max(freqvec_wl, na.rm = TRUE)%00%0 # KCE

  # REPetetiveness
  REP_av  <- (N_hlp+N_vlp) / N_dlp
  REP_hl  <-  N_hlp/N_dlp
  REP_vl  <-  N_vlp/N_dlp

  #Coefficient of determination
  CoV_dl = stats::sd(dlines)/mean(dlines)
  CoV_vl = stats::sd(vlines)/mean(vlines)
  CoV_hl = stats::sd(hlines)/mean(hlines)
  CoV_hv = stats::sd(c(hlines,vlines))/mean(c(hlines,vlines))
  CoV_wl = stats::sd(wlines)/mean(wlines) # KCE

  #Divergence
  DIV_dl = 1/MAX_dl
  DIV_vl = 1/MAX_vl
  DIV_hl = 1/MAX_hl
  DIV_hv = 1/MAX_hv
  DIV_wl = 1/MAX_wl # KCE

  # KCE
  # Find singular points by finding which points have no neighbours (left, right, upper, lower)
  N <- ncol(RM); E <- rep(0, N) # Empty vector
  RM_neighbours = rbind(E, RM[1:(N-1),]) + rbind(RM[2:N,], E) + cbind(E, RM[,1:(N-1)]) + cbind(RM[,2:N], E)
  SING_N = sum((RM_neighbours == 0) * RM)
  SING_RATE =  SING_N / RP_N

  out <- data.frame(
    RP_N      = RP_N,
    RR = RR,
    RP_max = RP_max,
    DIV_dl    = DIV_dl,
    REP_av    = REP_av,
    N_dl      = N_dl,
    N_dlp     = N_dlp,
    DET       = DET,
    MEAN_dl   = MEAN_dl,
    MAX_dl    = MAX_dl,
    ENT_dl    = ENT_dl,
    ENTrel_dl = ENTrel_dl,
    CoV_dl    = CoV_dl,
    N_vl      = N_vl,
    N_vlp     = N_vlp,
    LAM_vl    = LAM_vl,
    TT_vl     = MEAN_vl,
    MAX_vl    = MAX_vl,
    ENT_vl    = ENT_vl,
    ENTrel_vl = ENTrel_vl,
    CoV_vl    = CoV_vl,
    REP_vl    = REP_vl,
    N_hlp     = N_hlp,
    N_hl      = N_hl,
    LAM_hl    = LAM_hl,
    TT_hl     = MEAN_hl,
    MAX_hl    = MAX_hl,
    ENT_hl    = ENT_hl,
    ENTrel_hl = ENTrel_hl,
    CoV_hl    = CoV_hl,
    REP_hl    = REP_hl,
    N_hvp     = N_hvp,
    N_hv      = N_hv,
    LAM_hv    = LAM_hv,
    TT_hv     = MEAN_hv,
    MAX_hv    = MAX_hv,
    ENT_hv    = ENT_hv,
    ENTrel_hv = ENTrel_hv,
    CoV_hv    = CoV_hv,
    N_wl = N_wl, # KCE
    N_wlp = N_wlp, # KCE
    LAM_wl = LAM_wl, # KCE
    ENT_wl = ENT_wl, # KCE
    ENTrel_wl = ENTrel_wl, # KCE
    MEAN_wl = MEAN_wl, # KCE
    MAX_wl = MAX_wl, # KCE
    CoV_wl = CoV_wl, # KCE
    DIV_wl = DIV_wl, # KCE,
    LAM_DET = LAM_vl / DET, # KCE,
    DET_RR = DET / RR, # KCE,
    DIV_vl = 1 / MAX_vl, # KCE,
    SING_N = SING_N,
    SING_RATE = SING_RATE
  )

  if(matrices){
    return(list(
      crqaMeasures = out,
      crqaMatrices = list(RM     = RM,
                          dlines = dlines,
                          vlines = vlines,
                          hlines = hlines,
                          wlines = wlines, # KCE
                          freq_dl = freq_dl,
                          freq_vl = freq_vl,
                          freq_hl = freq_hl,
                          freq_wl = freq_wl)) # KCE
    )
  } else {
    return(out)
  }
}






#' Extract line lengths from recurrence matrix
#'
#' @inheritParams rp_measures_KCE
#' @inheritParams rp_calc_lineMeasures_KCE
#'
#' @description Extract lengths of diagonal, vertical and horizontal line segments from a recurrence matrix.
#'
#' @details Based on the Matlab function `linedists` by Stefan Schinkel, Copyright (C) 2009 Stefan Schinkel, University of Potsdam, http://www.agnld.uni-potsdam.de
#'
#' References:
#' S. Schinkel, N. Marwan, O. Dimigen & J. Kurths (2009):
#' "Confidence Bounds of recurrence-based complexity measures
#' Physics Letters A,  373(26), pp. 2245-2250
#'
#' Copyright (C) 2009 Stefan Schinkel, University of Potsdam
#' <http://www.agnld.uni-potsdam.de>
#'
#' @author Fred Hasselman
#' @return A list object with distributions of line lengths. If `matrices = TRUE` dataframes are returned whose columns represent the nonzero diagonals, verticals, or, horizontals.
#' @importFrom invctr  `%[]%`
#' @export
rp_lineDist_KCE <- function(RM,
                            DLmin = 2,
                            VLmin = 2,
                            HLmin = 2,
                            WLmin = 1, # KCE
                            DLmax = length(Matrix::diag(RM))-1,
                            VLmax = length(Matrix::diag(RM))-1,
                            HLmax = length(Matrix::diag(RM))-1,
                            WLmax = length(Matrix::diag(RM))-1, # KCE
                            d         = NULL,
                            theiler   = NA,
                            # recurrenceTimes = FALSE, # KCE
                            AUTO      = NULL,
                            matrices  = FALSE){

  if(!all(as.vector(RM)==0|as.vector(RM)==1)){stop("Matrix should be a binary (0,1) matrix!")}

  if(length(d)<=2){
    suppressMessages(RM <- casnet::setTheiler(RM, theiler))
  } else {
    if(!is.null(attributes(RM)$theiler)){
      message(paste0("Value found in attribute 'theiler'... assuming a theiler window was already applied to the matrix."))
    }
  }

  if(!is.null(d)){
    if(length(d)==1){d <- seq(-d,d)}
    if(length(d)==2){d <- seq(min(d),max(d))}
  }


  B <- casnet::rp_nzdiags(RM)
  V <- Matrix::as.matrix(RM)[,colSums(Matrix::as.matrix(RM))>0]

  RPt <- Matrix::t(RM)

  H <- Matrix::as.matrix(RPt)[,colSums(Matrix::as.matrix(RPt))>0]
  rm(RPt)

  # Get white line distribution (recurrence time of the first type, see Figure 1.3b in Marwan and Webber, 2015, Chapter 1) # KCE
  RM_inv <- Matrix::Matrix(1-RM,sparse = TRUE) %>% casnet::setTheiler(theiler = theiler) # Invert recurrence matrix; make sure diagonal and Theiler window is still set to 0

  W <- Matrix::as.matrix(RM_inv)[,colSums(Matrix::as.matrix(RM_inv))>0]
  rm(RM_inv)

  if(NCOL(B)==0|is.null(dim(B))){

    B <- matrix(0)
  }
  if(NCOL(V)==0|is.null(dim(V))){
    V <- matrix(0)
  }
  if(NCOL(H)==0|is.null(dim(H))){
    H <- matrix(0)
  }
  if(NCOL(W)==0|is.null(dim(W))){ # KCE
    W <- matrix(0)
  }

  # Get diagonal lines & pad with zeros
  diagonals   <- rbind.data.frame(rep(0,dim(B)[2]),
                                  B,
                                  rep(0,dim(B)[2])
  )

  # get nonzero vertical Lines & pad with zeros
  verticals <- rbind.data.frame(rep(x = 0, times = dim(V)[2]),
                                V,
                                rep(x= 0, times = dim(V)[2])
  )

  # get nonzero horizontal Lines & pad with zeros
  horizontals <- rbind.data.frame(rep(0,dim(H)[2]),
                                  H,
                                  rep(0,dim(H)[2])
  )

  # get nonzero white Lines & pad with zeros # KCE
  whites <- rbind.data.frame(rep(x = 0, times = dim(W)[2]),
                             W,
                             rep(x= 0, times = dim(W)[2])
  )

  # Get indices of line lengths
  diagonals.ind   <- tidyr::gather(diagonals,   key = "diagonal",   value = "segment")
  verticals.ind   <- tidyr::gather(verticals,   key = "vertical",   value = "segment")
  horizontals.ind <- tidyr::gather(horizontals, key = "horizontal", value = "segment")
  whites.ind   <- tidyr::gather(whites,   key = "whites",   value = "segment") # KCE

  invert = FALSE
  D <- diagonals.ind$segment
  names(D) <- paste0(diagonals.ind$diagonal,ifelse(invert,"DT","D"))
  V <- verticals.ind$segment
  names(V) <- paste0(verticals.ind$vertical,ifelse(invert,"VT","V"))
  H <- horizontals.ind$segment
  names(H) <- paste0(horizontals.ind$horizontal,ifelse(invert,"HT","H"))
  W <- whites.ind$segment # KCE
  names(W) <- paste0(whites.ind$whites,ifelse(invert,"WT","W")) # KCE

  # Get consecutive nonzero segments from indices, their difference is the segment length
  # We added a row of 0s so we'll get sequences of -1, 1, -1
  diagonals.dist   <- sort(which(diff(D)==-1)-which(diff(D)==1))
  verticals.dist   <- sort(which(diff(V)==-1)-which(diff(V)==1))
  horizontals.dist <- sort(which(diff(H)==-1)-which(diff(H)==1))
  whites.dist   <- sort(which(diff(W)==-1)-which(diff(W)==1)) # KCE

  diagonals.dist   <- diagonals.dist[diagonals.dist%[]%c(DLmin,DLmax)]
  verticals.dist   <- verticals.dist[verticals.dist%[]%c(VLmin,VLmax)]
  horizontals.dist <- horizontals.dist[horizontals.dist%[]%c(HLmin,HLmax)]
  whites.dist   <- whites.dist[whites.dist%[]%c(WLmin,WLmax)] # KCE

  if(length(diagonals.dist)==0){diagonals.dist <- NA}
  if(length(verticals.dist)==0){verticals.dist <- NA}
  if(length(horizontals.dist)==0){horizontals.dist <- NA}
  if(length(whites.dist)==0){whites.dist <- NA} # KCE

  return(list(diagonals.dist   = diagonals.dist,
              verticals.dist   = verticals.dist,
              horizontals.dist = horizontals.dist,
              whites.dist = whites.dist, # KCE
              diagonals.mat = diagonals[-c(1,NROW(diagonals)),],
              verticals.mat = verticals[-c(1,NROW(verticals)),],
              horizontals.mat = horizontals[-c(1,NROW(horizontals)),])[c(TRUE,TRUE,TRUE,TRUE,matrices,matrices,matrices,matrices)] # KCE
  )
}

