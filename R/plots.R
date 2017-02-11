# Plotting Utilities
# 
# Author: Renaud Gaujoux
###############################################################################


#' Open a File Graphic Device
#'
#' Opens a graphic device depending on the file extension.
#' 
#' @param filename path to the image file to create.
#' @param width output width
#' @param height output height
#' @param ... other arguments passed to the relevant device function
#' such as \code{\link{png}} or \code{\link{pdf}}.
#' 
#' importFrom grDevices bmp jpeg pdf png svg tiff
#' @export
gfile <- function(filename, width, height, ...){ 
  # Get file type
  r = regexpr("\\.[a-zA-Z]*$", filename)
  if(r == -1) stop("Improper filename")
  ending = substr(filename, r + 1, r + attr(r, "match.length"))
  
  f = switch(ending,
      pdf = function(x, ...) pdf(x, ...),
      svg = function(x, ...) svg(x, ...),			
      png = function(x, ...) png(x, ...),
      jpeg = function(x, ...) jpeg(x, ...),
      jpg = function(x, ...) jpeg(x, ...),
      tiff = function(x, ...) tiff(x, compression = "lzw", ...),
      bmp = function(x, ...) bmp(x, ...),
      stop("File type should be: pdf, svg, png, bmp, jpg, tiff")
  )
  
  args <- c(list(filename), list(...))	
  if( !missing(width) ){
    args$width <- as.numeric(width)
    args$height <- as.numeric(height)
    if( !ending %in% c('pdf','svg') && is.null(args[['res']]) ){
      args$units <- "in"
      args$res <- 300
    }
  }
  do.call('f', args)	
}

