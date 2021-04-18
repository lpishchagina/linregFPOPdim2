#' @title lrdata_gen2D
#'  
#' @description Generation of data (x, y) when x~N(xmean,1) and y = kCoef*x + aCoef+e, e~N(0,1) with a given value of xmean and vectors of kCoef, aCoef and changepoints  
#' 
#' @param n number of data point.
#' @param chpts a vector of increasing changepoint indices. Last index is always less than 'n'.
#' By default, 'chpts = NULL' (the data without changepoints). 
#' @param kcoef vector of successive kCoef for y, by default 'kCoef = seq(from = 1, to = length(chps)+1, by = 1)'.
#' @param aCoef vector of successive aCoef for y, by default 'aCoef = seq(from = 1, to = length(chps)+1, by = 1)'.
#' @param xmean value of mean for x, by default 'xmean = 0'.
#' @param noise standard deviation of an additional normal noise, by default 'noise = 1'.
#'  
#' @return matrix of data of dimension 2 x n with a given values of mean for x and kCoef, aCoef for linear regression  by the segmentation.
#'  
#' @examples
#' Data1 <- lrdata_gen2D(n = 10, chpts = NULL, xmean = 0, kCoef = 2, aCoef = 1, noise = 1)
#' Data2 <- lrdata_gen2D(n = 10, chpts = c(5), kCoef = c(10,20), aCoef = c(1,10), noise = 1)

lrdata_gen2D <- function(n, chpts = NULL, xmean = 0, kCoef = seq(from = 1, to = length(chps)+1, by = 1), aCoef = seq(from = 1, to = length(chps)+1, by = 1), noise = 1)
{
  #---stop---#
  if (!is.null(chpts) && n <= chpts[length(chpts)]){stop('last element of changepoints is always less than n')}
  if(!is.null(chpts) && !is.numeric(chpts)){stop('changepoints are not all numeric')}
  if(is.unsorted(chpts)){stop('changepoints should be an increasing vector')}
  
  if(!is.numeric(xmean)){stop('xmean is not all numeric')}
  if(!is.numeric(kCoef)){stop('kCoef are not all numeric')}
  if(!is.numeric(aCoef)){stop('aCoef are not all numeric')}
  
  if ((length(chpts)+1) !=  length(aCoef)){stop('The length of the aCoef is always equal to the number of changepoints plus one')}
  if ((length(chpts)+1) !=  length(kCoef)){stop('The length of the kCoef is always equal to the number of changepoints plus one')}  
  
  if(length(kCoef) != (length(aCoef))){stop('kCoef and aCoef vectors are of different size')}
  
  if(!is.double(noise)){stop('noise is not a double')}
  if(noise < 0){stop('noise must be non-negative')}
  #---function---#	
  
  data <- matrix(0,2,n)
  InttT<- diff(c(0,chpts,n))
  # rnorm(mu,noise) = mu + rnorm(0,noise)
  data[1,] <- rnorm(n, xmean, noise)
  k <-  rep(kCoef, InttT)
  a <-  rep(aCoef, InttT)
  e <- rnorm(n, 0, noise)
  data[2,] <- k* data[1,]+a
  return(data)
}






