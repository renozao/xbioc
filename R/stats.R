# Common statistical utilities
# 
# Author: Renaud Gaujoux
###############################################################################


#' Compute Permutation-based FDR as in SAM
#' 
#' @param obs observed statistics
#' @param perms matrix of statistics obtained for each permutation (column).
#' @param alternative alternative hypothesis to test
#' @param digit FDR precision in number of digits
#' 
#' @export 
samFDR <- function(obs, perms, alternative = c('two.sided', 'greater', 'less'), digit = 2) {
  
  tt.sam <- obs
  ttstar.sam <- perms
  numgene <- length(tt.sam)
  nperms <- nrow(perms)
  N <- 10^digit
  cutp.sam <- seq(0, max(abs(tt.sam), na.rm = TRUE), length = N)
  sigGene.sam <- setNames(rep(NA_real_, numgene), names(obs))
  fdr.sam <- ncall.sam <- rep(NA_real_, 100)
  alternative <- match.arg(alternative)
  
  lapply(1:N, function(i){
        th <- cutp.sam[i]
        # compute number of hits in original and permutations
        if(alternative == 'two.sided') {
          i_call <- which(abs(tt.sam) >= th)
          n_perms <- sum(abs(ttstar.sam) >= cutp.sam[i], na.rm = TRUE)
          
        } else if(alternative == 'greater') {
          i_call <- which(tt.sam >= th)
          n_perms <- sum(ttstar.sam >= cutp.sam[i], na.rm = TRUE)
          
        } else if(alternative == 'less') {
          i_call <- which(tt.sam <= -th)
          n_perms <- sum(ttstar.sam <= -th, na.rm = TRUE)
        }
        
        # compute FDR
        n_call <- length(i_call)
        fdr <- n_perms / nperms / n_call
        fdr.sam[i] <<- fdr
        ncall.sam[i] <<- n_call
        sigGene.sam[i_call] <<- i
      })
  
  fdr.sam <- pmin(fdr.sam,1)
  fdr.sam <- make.monotone(fdr.sam)
  
  return (list(fdr.sam = fdr.sam, ncall.sam = ncall.sam, sigGene.sam = setNames(fdr.sam[sigGene.sam], names(sigGene.sam))))
}

make.monotone <-
    function(x) {
  n=length(x)
  for(j in 2:n){ x[j]=min(x[j-1],x[j])}
  return(x)
}


