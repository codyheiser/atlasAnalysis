#' Tissue Clustering Utility Functions
#' @author: C Heiser
#' adapted from S Vandekar `tissue_clustering.Rmd`
library(raster)
set.seed(18)


downsample <- function(img, outimg, mask=NULL, fact=10, fun='mean'){
  #' Downsamples image
  #'
  #' @param img image to downsample as filepath or raster
  #' @param outimg path to image file to write downsampled output
  #' @param mask path to tissue mask corresponding to `img`. if NULL, don't mask `img`
  #' prior to downsampling
  #' @param fact factor to downsample by in x and y pixel directions
  #' @param fun function to pass to `aggregate` to downsample pixels
  #'
  #' @return ans writes downsampled `img` to `outimg`
  if(is.character(img)){
    message(basename(img))
    img = raster(img)
  }
  if(!is.null(mask)){
    if(is.character(mask)) mask = raster(mask)
    if(all(dim(img)==ceiling(dim(mask)/2))){
      mask = aggregate(mask, fact=2, fun=fun)
      fact = fact/2
      extent(img) = extent(mask)
    }
  }
  ans = aggregate(img, fact=fact, fun=fun, filename=outimg, overwrite=TRUE)
  ans = outimg
}


gaussSmooth <- function(
  img,
  outimg,
  mask,
  maskThr,
  sigma=5,
  transform=log10,
  invtransform=function(x) 10^x,
  offset=1,
  divmean = function(img){ cellStats(img, stat='mean')}
){
  #' Smooths image with simple Gaussian kernel on a scale defined by `transform`
  #'
  #' @param img image to smooth as filepath or raster
  #' @param outimg path to image file to write smoothed output
  #' @param mask path to tissue mask corresponding to `img`. if NULL, don't mask `img`
  #' prior to smoothing
  #' @param maskThr threshold below which to exclude pixels from smoothing based on
  #' `mask` values
  #' @param sigma size of Gaussian kernel in x and y pixel directions
  #' @param transform function to transform pixel values with before smoothing
  #' @param invtransform function to undo transformation via `transform` to get pixel
  #' values back to original space
  #' @param offset pseudovalue to add prior to transform to avoid `transform`(0)
  #' @param divmean function to calculate mean to divide pixel values by when
  #' transforming
  #'
  #' @return ans writes smoothed `img` to `outimg`
  message(basename(img))
  ans = raster(img)
  maskimg = raster(mask)
  values(ans)[values(maskimg)<=maskThr] = NA
  values(ans) = transform(values(ans)/divmean(ans) + offset)
  if(sigma>0){
    # not a Gaussian smooth, but can handle boundaries better
    ans = focal(
      ans,
      focalWeight(ans, sigma, type='circle'),
      filename=outimg,
      overwrite=TRUE,
      fun='sum',
      na.rm=TRUE
    )
    #fun=function(x, na.rm=TRUE){ prod(x+1, na.rm=na.rm)^(1/sum(!is.na(x)))-1}
    # set NAs to zero
    values(ans) = invtransform(values(ans))-offset
    values(ans)[is.na(values(ans))] = 0
    writeRaster(ans, filename=outimg, overwrite=TRUE)
  }
  ans = outimg
}


subsampImgs <- function(
  imgs,
  masks,
  maskThr=0,
  subsamp=4,
  transform=log10,
  offset=1
){
  #' Gets subset of data from the images
  #'
  #' @param imgs character vector, matrix, or data frame of images to load. rows are
  #' slides/regions/subjects; columns are images.
  #' @param masks character vector of masks defining where there is tissue in the image
  #' @param maskThr numeric noninclusive lower threshold for mask
  #' @param subsamp integer specifying the factor to subsample each dimension of the
  #' image.
  #' @param transform function to transform pixel values with
  #' @param offset pseudovalue to add prior to transform to avoid `transform`(0)
  #'
  #' @return locData dataframe with the subsample of data for all of the images.
  #' images are identified by the rownames of `imgs`
  imgs = as.data.frame(imgs)
  ss = 1:nrow(imgs)
  imgnames = names(imgs)
  # samples every 4th downsampled voxel from the DAPI
  subsamp = rep(c(FALSE, TRUE), c(subsamp-1, 1))
  locData = do.call(rbind, lapply(ss, function(rind){
    message(rownames(imgs)[rind])
    # get DAPI for this slideRegion
    img = raster(masks[rind])
    arrIndsImg = img
    values(arrIndsImg) = outer(
      rep(subsamp, length.out=nrow(img)),
      rep(subsamp, length.out=ncol(img)),
      FUN="&"
    )
    # indsImg selects the pixels that will be used to build the model
    indsImg = (img>dapiThr & arrIndsImg)
    rm(arrIndsImg, img)
    # for channel
    names(imgnames) = imgnames
    res = as.data.frame(do.call(cbind, mclapply(imgnames, function(marker){
      # gets row for this subject
      img = raster(imgs[rind, marker])
      # select target voxels
      res = values(img)[ values(indsImg) ]
      res = transform(res + offset)
    }, mc.cores = 2 )))
    # set NAs
    res[ , imgnames[which(!imgnames %in% names(res)) ] ] = NA
    # reorder them to be the original
    res = res[ ,imgnames]
    # return result for this slideRegion
    res[,'ID'] = rownames(imgs)[rind]
    res
  }
  ) ) # combine them all into one data frame
}


labelImage <- function(
  imgnames,
  centers,
  weights,
  mask,
  maskThr=0,
  outname,
  transform=log10,
  offset=1
){
  #' Takes k-means centers and assigns each pixel to one of the centroids
  #' Always divides the image by the mean, adds offset, and applies transformation
  #'
  #' @param imgnames tif image names for markers used in the clustering
  #' @param centers kmeans centroid coordinates, with columns in the same order as
  #' imgnames
  #' @param weights Multiplicative weights applied in model fitting
  #' @param mask path to tissue mask corresponding to `imgnames`
  #' @param maskThr threshold below which to exclude pixels from labeling based on
  #' `mask` values
  #' @param outname filename for output image
  #' @param transform function to transform pixel values with before assigning labels
  #' @param offset pseudovalue to add prior to transform to avoid `transform`(0)
  #'
  #' @return outname writes image with cluster labels to `outname`
  imgs = stack(imgnames)
  mask = raster(mask)
  values(imgs)[values(mask)<=maskThr] = NA
  imgout = calc(imgs, fun=function(vec){
    nas = apply(is.na(vec), 1, any)
    ans = rep(0, nrow(vec))
    vec = transform(vec[!nas,]+offset)
    if(nrow(vec)!=0){
      ans[!nas] = as.numeric(
        as.character(
          c(knn1(centers, test=sweep(vec, 2, weights, '*'), cl=1:nrow(centers)))
        )
      )
    }
    return(ans)
  }, filename=outname, overwrite=TRUE)
  return(outname)
}


moran <- function (x, y=NULL, atlas, w = matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), 3, 3)){
  #' Computes Moran Index
  #' TODO: check that this does the same thing as Moran function when NAs are present
  #' reference:
  #' https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0126158
  #'
  #' @param x
  #' @param y
  #' @param atlas
  #' @param w
  #'
  #' @return mI Moran Index
  if(is.null(y)) y = x
  if(is.character(x)) x = raster(x)
  if(is.character(y)) y = raster(y)
  if(is.character(atlas)) atlas = raster(atlas)
  z <- x - cellStats(x, mean)
  z2 <- y - cellStats(y, mean)
  wZiZj <- focal(z2, w = w, fun = "sum", na.rm = TRUE, pad = TRUE)
  wZiZj <- wZiZj * z
  wZiZj <- c(by(values(wZiZj), values(atlas), sum))
  zsq <- sqrt(c(by(values(z^2), values(atlas), sum)))
  # ineffecient if x=y
  z2sq <- sqrt(c(by(values(z2^2), values(atlas), sum)))
  n <- table(values(atlas)) - c(by(values(is.na(z) | is.na(z2)), values(atlas), sum))
  values(z) = as.numeric(!is.na(values(z)) & !is.na(values(z2)))
  W <- c(
    by(
      values(
        focal(z, w = w, fun = "sum", na.rm = TRUE, pad = TRUE)
      ),
      values(atlas),
      sum
    )
  )
  NS0 <- n/W
  mI <- NS0 * wZiZj/zsq/z2sq
  return(mI)
}


snr <- function(imgname, mask, maskThr=0, transform=log10, offset=1){
  #' Computes signal-to-noise ratio
  #'
  #' @param imgname image to calculate snr as filepath
  #' @param mask filepath to tissue mask corresponding to `imgname`
  #' @param maskThr numeric noninclusive lower threshold for mask
  #' @param transform function to transform pixel values with
  #' @param offset pseudovalue to add prior to transform to avoid `transform`(0)
  #'
  #' @return snr list with `mu` and `sigma` values from snr
  imgs = raster(imgname)
  mask = raster(mask)
  values(imgs)[ values(mask)<=maskThr ] = NA
  mu = cellStats(imgs, 'mean')
  sigma = cellStats(imgs, stat='sd')
  #values(imgs) = transform(values(imgs)/mu+offset)
  #snrTransform = cellStats(imgs, 'mean')/cellStats(imgs, 'sd')
  c(mu=mu, sigma=sigma)
}
