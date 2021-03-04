#' @title Bootstrap dissimilarity
#'
#' @description Bootstrap dissimilarity between two vectors of interactions (or
#'   species names).
#'
#' @param lst_vect A list with two character vectors of interactions.
#'
#' @param index The name of the desired index, e.g.: "whittaker", "jaccard",
#'   "sorensen". Character type.
#'
#' @param by Integer. Number of interactions used to increase gradually the
#'   sampled vectors of interactions until all observations are sampled. If `by`
#'   is too small (e.g. 1) then the computation time is very long depending on
#'   your total number of interactions from which samples are taken. As a rule
#'   of thum set `by` to maybe 5-10\\% of your total interactions.
#'
#' @param replace Should sampling be with replacement? Passed to
#'   \code{\link[base]{sample}}.
#'
#' @param size_short Should the total sample size be set by the length of the
#'   shortest vector from `lst_vect`? Logic type (`TRUE` or `FALSE`).
#'
#' @param n_boot Number of desired bootstraps (100 can be enough).
#'
#' @param n_cpu Number of CPU-s to use for parallel processing.
#'
#' @return A matrix of dissimilarities. The number of columns corresponds to
#'   `n_boot` (number of bootstraps). The names of the rows store the sample
#'   size of each sampled chunk. The final value corresponds to the length of
#'   the shortest vector (when `size_short = TRUE`), or to the length of the
#'   longest vector (when `size_short = FALSE`).
#'
#' @examples
#'
#' library(bootdissim)
#' library(bipartite)
#'
#' net1 <- reshape_net(vazarr, seed = 1)
#' net2 <- reshape_net(vazcer, seed = 1)
#'
#' vect <- get_interactions(net1, net2)
#'
#' tst <- boot_dissim(lst_vect = vect,
#'                    index = "whittaker",
#'                    by = 50,
#'                    replace = TRUE,
#'                    size_short = FALSE,
#'                    n_boot = 10,
#'                    n_cpu = 2)
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom parallel splitIndices makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#'
#' @export
#'
#' @md
boot_dissim <- function(lst_vect,
                        index,
                        by,
                        replace,
                        size_short,
                        n_boot,
                        n_cpu){

  # Detect which of the two vectors is the short one and which is the long one.
  which_short <- which.min(sapply(X = lst_vect, FUN = length))
  which_long <- which(1:2 != which_short)
  vect_short <- lst_vect[[which_short]]
  vect_long <- lst_vect[[which_long]]

  # If the user chooses to use as total sample size the length of the shortest
  # vector, then set is so.
  size <- ifelse(isTRUE(size_short),
                 yes = length(vect_short),
                 no = length(vect_long))

  # Sample the vectors here to get the sample sizes of each chunk so that can be
  # passed later on easily to the results and will end up on the OX axis of the
  # graphs. Avoid to call this in the loop, because it stays constant anyways
  # and otherwise consumes unnecessary time.
  chunks_lst <- sample_vectors(v1 = vect_short,
                               v2 = vect_long,
                               size = size,
                               by = by,
                               seed = 42,
                               replace = replace)
  n <- length(chunks_lst[[1]])
  # same as length(chunks_lst[[2]])
  spl_size <- sapply(chunks_lst[[1]], FUN = length)
  # same as sapply(chunks_lst[[2]], FUN = length)

  # Start parallel processing
  chunks <- parallel::splitIndices(n_boot, n_cpu)
  cl <- parallel::makeCluster(n_cpu)
  doParallel::registerDoParallel(cl)
  i <- NULL # to avoid 'Undefined global functions or variables: i' in R CMD check
  # Compute dissimilarities in parallel
  boot_lst <-
    foreach::foreach(i = iterators::iter(chunks),
                     .errorhandling = 'pass',
                     .export = c("boot_dissim_once",
                                 "sample_vectors",
                                 "get_abc",
                                 "get_beta")) %dopar% {
                                   lapply(i, FUN = function(x) # note the lapply!
                                     boot_dissim_once(vect_short = vect_short,
                                                      vect_long = vect_long,
                                                      index = index,
                                                      size = size,
                                                      by = by,
                                                      replace = replace,
                                                      seed = x, # iterator is passed as seed (values are 1:n_boot)
                                                      n = n)
                                   )
                                 }
  parallel::stopCluster(cl)
  remove(cl)
  foreach::registerDoSEQ()
  # End of parallel processing

  # Prepare results as a matrix. The row names give the sample size of each
  # chunk. There are as many column as the number of bootstraps (n_boot).
  boot_lst <- unlist(boot_lst, recursive = FALSE)
  results <- do.call(what = "cbind", args = boot_lst)
  rownames(results) <- spl_size

  return(results)
}


#' @title Compute dissimilarity for a sample
#'
#' @description Compute dissimilarity between two vectors of interactions (or
#'   species names) for a single run (sample). You will rarely use this function
#'   alone. It was designed to be executed in parallel by
#'   \code{\link[bootdissim]{boot_dissim}}. See more details there.
#'
#' @param vect_short The shorter vectors of interactions.
#'
#' @param vect_long The longer vectors of interactions.
#'
#' @param index,size,by,replace See \code{\link[bootdissim]{boot_dissim}}.
#'
#' @param seed Passed to \code{\link[base]{set.seed}}. Set seed to get
#'   reproducible random results. Is give by the iterator of the loop during the
#'   parallel processing in \code{\link[bootdissim]{boot_dissim}}.
#'
#' @param n Passed from \code{\link[bootdissim]{boot_dissim}}. Is constant
#'   across all bootstrap iterations.
#'
#' @return A vector of dissimilarities. The length of this vector corresponds to
#'   `n`.
#'
#' @export
#'
#' @md
boot_dissim_once <- function(vect_short,
                             vect_long,
                             index,
                             size,
                             by,
                             replace,
                             seed,
                             n){

  chunks_lst <- sample_vectors(v1 = vect_short,
                               v2 = vect_long,
                               size = size,
                               by = by,
                               seed = seed,
                               replace = replace)

  dissim <- rep(NA, n)

  for (i in 1:n){
    v1 <- chunks_lst[["v1"]][[i]]
    v2 <- chunks_lst[["v2"]][[i]]
    abc <- get_abc(list(v1, v2))
    dissim[[i]] <- get_beta(abc, index)
  }

  return(dissim)
}
