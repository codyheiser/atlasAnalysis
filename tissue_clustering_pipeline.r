#' Tissue Clustering Pipeline
#' @author: C Heiser
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(source("~/git/atlasAnalysis/tissue_clustering_utils.r"))


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
  kMeansMarkers <- paste(markers, 'smooth', sep="_")
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
    description='Perform pixel-level clustering of tissue regions from MxIF data'
  )
  parser$add_argument(
    "wd",
    type="character",
    help="Base working directory containing image files for each slide"
  )
  parser$add_argument(
    "slides",
    nargs="+",
    help="Names of slides to use for clustering"
  )
  parser$add_argument(
    "-m",
    "--markers",
    dest="markers",
    nargs="+",
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
    default=8,
    help="Factor in pixels by which to downsample images"
  )
  parser$add_argument(
    "-s",
    "--sigma",
    dest="sigma",
    default=50,
    help="Radius of disk in pixels used to smooth downsampled images"
  )
  parser$add_argument(
    "-k",
    "--k-clusters",
    dest="k",
    default=4,
    help="Number of k-means clusters"
  )
  parser$add_argument(
    "-j",
    "--n-jobs",
    dest="njobs",
    default=4,
    help="Number of CPU cores to use for parallelization"
  )
  parser$add_argument(
    "-p",
    "--prefix",
    dest="prefix",
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
