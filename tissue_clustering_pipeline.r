#' Tissue Clustering Pipeline
#' @author: C Heiser
#' adapted from S Vandekar `tissue_clustering.Rmd`
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(class))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(raster))


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
    indsImg = (img>maskThr & arrIndsImg)
    rm(arrIndsImg, img)
    # for channel
    names(imgnames) = imgnames
    res = as.data.frame(do.call(cbind, mclapply(imgnames, function(marker){
      # gets row for this subject
      img = raster(imgs[rind, marker])
      # select target voxels
      res = values(img)[values(indsImg)]
      res = transform(res + offset)
    }, mc.cores = 2 )))
    # set NAs
    res[, imgnames[which(!imgnames %in% names(res))]] = NA
    # reorder them to be the original
    res = res[, imgnames]
    # return result for this slideRegion
    res[, 'ID'] = rownames(imgs)[rind]
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


setup_data <- function(
  wd,
  slides,
  markers=c("COLLAGEN","DAPI","MUC2","NAKATPASE","OLFM4","PANCK","SOX9","VIMENTIN")
){
  #' Sets up data structure to loop through when clustering images
  #'
  #' @param wd path to base directory containing image data
  #' @param slides list of slide names (subdirectories of `wd`) to cluster with
  #' @param markers list of markers to use for tissue clustering
  #'
  #' @return sr
  sr <- list()
  for(i in 1:length(slides)){
    file <- file.path(wd, sprintf('%s/%s_thresholds.csv', slides[i], slides[i]))
    sr[[i]] <- read.csv(file)
  }
  sr <- do.call(rbind,sr)
  sr$slideRegion <- sprintf('%s_region%s',sr$Slide,sr$Region)
  sr <- sr[!duplicated(sr$slideRegion), c('Slide', 'Region', 'slideRegion')]
  # set up expected file paths
  sr[,markers] <- sapply(
    markers,
    function(m) file.path(
      wd,
      sr$Slide,
      sprintf('%s_ADJ_region_%03d.tif', m, sr$Region)
    )
  )
  sr[,'mask'] <- file.path(
    wd,
    sr$Slide,
    sprintf('%s_%02d_TISSUE_MASK.tif', sr$Slide, sr$Region)
  )
  sr[,'tumorMask'] <- file.path(
    wd,
    sr$Slide,
    sprintf('Tumor_mask_region_%03d.png', sr$Region)
  )
  sr[,'epiMask'] <- file.path(
    wd,
    sr$Slide,
    sprintf('%s_region_%03d_epi_mask.png', sr$Slide, sr$Region)
  )
  sr[,'strMask'] <- file.path(
    wd,
    sr$Slide,
    sprintf('%s_region_%03d_stroma_mask.png', sr$Slide, sr$Region)
  )
  # check that all the images exist
  test <- sapply(sr[,markers], file.exists)
  message(
    paste0(
      sum(test) / (dim(test)[1] * dim(test)[2]) * 100,
      " % of ",
      (dim(test)[1] * dim(test)[2]),
      " expected images were found in ",
      wd
    )
  )
  # only use regions with tissue mask
  sr = sr[file.exists(sr$mask),]
  return(sr)
}


downsample_data <- function(
  sr,
  fact=8,
  markers=c("COLLAGEN","DAPI","MUC2","NAKATPASE","OLFM4","PANCK","SOX9","VIMENTIN"),
  masks=c("mask","epiMask","strMask","tumorMask"),
  njobs=4,
  overwrite=T
){
  #' Downsamples images in data structure output from `setup_data`
  #'
  #' @param sr data structure output from `setup_data`
  #' @param fact factor to downsample in each dimension, in pixels
  #' @param markers list of markers to downsample
  #' @param masks list of masks to downsample
  #' @param njobs number of CPU cores to use in `mcmapply`
  #' @param overwrite force overwrite of downsampled files if they already exist
  #'
  #' @return sr
  origmarkers <- sr  # copy object for pulling full-res file names
  ss <- 1:nrow(sr)  # get counter to use in loops
  # downsample markers using average value across pixels
  message("Downsampling marker images:")
  for(marker in markers){
    sr[,marker] <- gsub('\\.png$|\\.tif$', '_downsampled.tif', origmarkers[,marker])
    if(!overwrite){
      if(!all(file.exists(sr[ss,marker]))){
        message(marker)
        mcmapply(
          downsample,
          img=origmarkers[ss,marker],
          outimg=sr[ss,marker],
          MoreArgs=list(
            fact=fact,
            mask=NULL,
            fun="mean"
          ),
          mc.cores=njobs
        )
      }
    }else{
      message(marker)
      mcmapply(
        downsample,
        img=origmarkers[ss,marker],
        outimg=sr[ss,marker],
        MoreArgs=list(
          fact=fact,
          mask=NULL,
          fun="mean"
        ),
        mc.cores=njobs
      )
    }
  }
  # downsample masks using mode
  message("Downsampling image masks:")
  for(marker in masks){
    sr[,marker] = gsub('\\.png$|\\.tif$', '_downsampled.tif', origmarkers[,marker])
    if(!overwrite){
      if(!all(file.exists(sr[ss,marker]))){
        message(marker)
        mcmapply(
          downsample,
          img=origmarkers[ss,marker],
          outimg=sr[ss,marker],
          MoreArgs=list(
            fact=fact,
            mask=NULL,
            fun="modal"
          ),
          mc.cores=njobs
        )
      }
    }else{
      message(marker)
      mcmapply(
        downsample,
        img=origmarkers[ss,marker],
        outimg=sr[ss,marker],
        MoreArgs=list(
          fact=fact,
          mask=NULL,
          fun="modal"
        ),
        mc.cores=njobs
      )
    }
  }
  message("\nDone!")
  return(sr)
}


smooth_data <- function(
  sr,
  sigma=50,
  markers=c("COLLAGEN","DAPI","MUC2","NAKATPASE","OLFM4","PANCK","SOX9","VIMENTIN"),
  njobs=4,
  overwrite=T
){
  #' Smooths images in data structure output from `setup_data` or `downsample_data`
  #'
  #' @param sr data structure output from `setup_data` or `downsample_data`
  #' @param sigma radius of circular area to blur each pixel with, in pixels
  #' @param markers list of markers to smooth
  #' @param njobs number of CPU cores to use in `mcmapply`
  #' @param overwrite force overwrite of smoothed files if they already exist
  #'
  #' @return sr
  sMarkers <- paste(markers, 'smooth', sep="_")  # create file names
  sr[,sMarkers] <- apply(
    sr[,markers],
    2,
    function(x) gsub('\\.tif$', '_smoothed.tif', x)
  )
  ss = 1:nrow(sr)  # get counter to use in loops
  # smooth the raw data; divides by the mean to normalize, no transformation
  message("Smoothing marker images:")
  for(marker in markers){
    if(!overwrite){
      if(!all(file.exists(sr[ss,marker]))){
        message(marker)
        mcmapply(
          gaussSmooth,
          img=sr[ss,marker],
          outimg=sr[ss,paste(marker,'smooth',sep="_")],
          mask=sr[ss,'mask'],
          MoreArgs=list(
            maskThr=5,
            sigma=sigma,
            transform=function(x) x,
            invtransform=function(x) x,
            offset=0
          ),
          mc.cores=njobs
        )
      }
    }else{
      message(marker)
      mcmapply(
        gaussSmooth,
        img=sr[ss,marker],
        outimg=sr[ss,paste(marker,'smooth',sep="_")],
        mask=sr[ss,'mask'],
        MoreArgs=list(
          maskThr=5,
          sigma=sigma,
          transform=function(x) x,
          invtransform=function(x) x,
          offset=0
        ),
        mc.cores=njobs
      )
    }
  }
  message("\nDone!")
  return(sr)
}


cluster_data <- function(
  sr,
  k=4,
  subsamp=8,
  transform=function(x) x,
  offset=0,
  markers=c(
    "COLLAGEN_smooth",
    "DAPI_smooth",
    "MUC2_smooth",
    "NAKATPASE_smooth",
    "OLFM4_smooth",
    "PANCK_smooth",
    "SOX9_smooth",
    "VIMENTIN_smooth"
  ),
  seed=NA,
  out.prefix='locationLabel'
){
  #' Identifies and labels k-means clusters using data structure output from
  #' `setup_data`, `downsample_data` or `smooth_data`
  #'
  #' @param sr data structure output from `setup_data`, `downsample_data` or
  #' `smooth_data`
  #' @param k number of clusters for k-means model
  #' @param subsamp ratio of total clusters from training images to subsample on grid
  #' for k-means training
  #' @param transform function to transform pixel values with prior to k-means
  #' @param offset pseudovalue to add prior to transform to avoid `transform`(0)
  #' @param markers list of markers to use for k-means model building
  #' @param seed random seed for k-means model
  #' @param out.prefix how to name the output cluster label files
  #'
  #' @return sr
  sr[,paste(out.prefix)] = gsub('DAPI', out.prefix, sr[,'DAPI'])
  ss = 1:nrow(sr)  # get counter to use in loops
  # subsample on grid, within 'mask' to get low number of pixels from each image
  message("Downsampling pixels for k-means training")
  kmeansData <- subsampImgs(
    imgs=sr[ss, markers],
    masks=sr[ss, 'mask'],
    maskThr=5,
    subsamp=subsamp,
    transform=transform,
    offset=offset
  )
  message("Fitting k-means model")
  # fit the kmeans model
  if(!is.na(seed)){
    set.seed(seed)
  }
  X = kmeansData[,markers]
  #X = apply(X, 2, function(x) ifelse(is.na(x), min(x, na.rm=TRUE), x))
  sds = apply(X, 2, sd)
  X = sweep(X, 2, sds, '/')
  kmeansMod = kmeans(X, k, iter.max=200, algorithm='MacQueen')
  # assign k-means labels
  message("Assigning k-means labels")
  outimgs <- sapply(
    ss,
    function(rind){
      message(sr$slideRegion[rind])
      labelImage(
        unlist(sr[rind,markers]),
        centers=kmeansMod$centers[,markers],
        weights=1/sds,
        mask=sr[rind, 'mask'],
        maskThr=5,
        transform=transform,
        offset=offset,
        outname=sr[rind,out.prefix]
      )
    }
  )
  message("\nDone!")
  return(sr)
}


tissue_clustering_pipe <- function(
  wd,
  slides,
  markers=c(
    "ACTININ",
    "BCATENIN",
    "COLLAGEN",
    "DAPI",
    "ERBB2",
    "MUC2",
    "NAKATPASE",
    "OLFM4",
    "PANCK",
    "PCNA",
    "SMA",
    "SOX9",
    "VIMENTIN"
  ),
  fact=8,
  sigma=50,
  k=4,
  njobs=4,
  overwrite=F,
  outprefix="clusters"
){
  #' Perform entire clustering pipeline
  #'
  #' @param wd path to base directory containing image data
  #' @param slides list of slide names (subdirectories of `wd`) to cluster with
  #' @param markers list of markers to use for tissue clustering
  #' @param fact factor to downsample in each dimension, in pixels
  #' @param sigma radius of circular area to blur each pixel with, in pixels
  #' @param k number of clusters for k-means model
  #' @param njobs number of CPU cores to use in `mcmapply`
  #' @param overwrite force overwrite of downsampled/smoothed files if they already
  #' exist
  #' @param out.prefix how to name the output cluster label files
  #'
  #' @return writes clustered images to files
  # data setup function to create sr
  setup_data(
    wd = wd,
    slides = slides,
    markers = markers
  ) %>%
  # downsampling function to loop through sr
  downsample_data(
    .,
    fact = fact,
    markers = markers,
    masks = c("mask"),
    njobs = 4,
    overwrite = overwrite
  ) %>%
  smooth_data(
    .,
    sigma = sigma,
    markers = markers,
    njobs = njobs,
    overwrite = overwrite
  ) -> sr
  # get names of smoothed markers for k-means
  kMeansMarkers <- paste(markers, "smooth", sep="_")
  # perform k-means clustering
  cluster_data(
    sr,
    k = k,
    subsamp = fact,
    transform = function(x) x,
    offset = 0,
    markers = kMeansMarkers,
    seed = NA,
    out.prefix = outprefix
  ) -> sr
}


if(!interactive()){
  parser <- ArgumentParser(
    description="Perform pixel-level clustering of tissue regions from MxIF data"
  )
  parser$add_argument(
    "wd",
    type="character",
    help="Base working directory containing image files for each slide"
  )
  parser$add_argument(
    "slides",
    nargs="+",
    type="character",
    help="Names of slides to use for clustering"
  )
  parser$add_argument(
    "-m",
    "--markers",
    dest="markers",
    nargs="+",
    type="character",
    default=c(
      "ACTININ",
      "BCATENIN",
      "COLLAGEN",
      "DAPI",
      "ERBB2",
      "MUC2",
      "NAKATPASE",
      "OLFM4",
      "PANCK",
      "PCNA",
      "SMA",
      "SOX9",
      "VIMENTIN"
    ),
    help="Names of protein markers to use for clustering"
  )
  parser$add_argument(
    "-d",
    "--downsample-factor",
    dest="downsamplefactor",
    type="integer",
    default=8,
    help="Factor in pixels by which to downsample images"
  )
  parser$add_argument(
    "-s",
    "--sigma",
    dest="sigma",
    type="integer",
    default=50,
    help="Radius of disk in pixels used to smooth downsampled images"
  )
  parser$add_argument(
    "-k",
    "--k-clusters",
    dest="k",
    type="integer",
    default=4,
    help="Number of k-means clusters"
  )
  parser$add_argument(
    "-j",
    "--n-jobs",
    dest="njobs",
    type="integer",
    default=4,
    help="Number of CPU cores to use for parallelization"
  )
  parser$add_argument(
    "-p",
    "--prefix",
    dest="prefix",
    type="character",
    default="clusters",
    help="Prefix to cluster-labeled image files written out"
  )
  parser$add_argument(
    "-o",
    "--overwrite",
    action="store_true",
    default=FALSE,
    help="Force overwrite of all downsampled and smoothed images before clustering"
  )
  args <- parser$parse_args()

  sr <- tissue_clustering_pipe(
    wd=args$wd,
    slides=args$slides,
    markers=args$markers,
    fact=args$downsamplefactor,
    sigma=args$sigma,
    k=args$k,
    njobs=args$njobs,
    overwrite=args$overwrite,
    outprefix=args$prefix
  )
}
