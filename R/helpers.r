#' @title Transform network (web) matrix to data table
#'
#' @description Helper function to transform a network matrix (web) to a data
#'   frame, where each row represents an interaction between a low level and
#'   high level species (e.g. plant - pollinator).
#'
#' @param net A network matrix (a web), e.g. `Safariland` from
#'   \code{\link{bipartite}}.
#'
#' @param abund Logic type (`TRUE` or `FALSE`). Default is `TRUE`, that is
#'   abundance information is considered and an interaction is repeated as many
#'   times it was observed. If `FALSE` then it returns only the unique
#'   interactions.
#'
#' @param seed A single integer value. Passed to \code{\link[base]{set.seed}}.
#'   Set seed to get reproducible random results.
#'
#' @return Returns a data frame.
#'
#' @examples
#'
#' library(bootdissim)
#' library(bipartite)
#'
#' # Considering abundance information
#' Safariland_dt <- reshape_net(Safariland, abund = TRUE, seed = 1)
#'
#' # Without abundance information
#' Safariland_dt <- reshape_net(Safariland, abund = FALSE, seed = 1)
#'
#' @importFrom data.table as.data.table melt setDF
#' @importFrom methods is
#'
#' @export
#'
#' @md
reshape_net <- function(net, abund = TRUE, seed){

  # Reduce dependencies to a minimum!
  # Maybe soething like:

  # net <- Safariland
  # net_df <- as.data.frame(net)
  # net_df$row_names <- rownames(net_df)
  # reshape(data = net_df, idvar = "row_names", direction = "long", varying = list(1:(ncol(net_df)-1)))

  # Test data type and typos in input
  if (! is(net, "matrix")) stop("`net` must be a matrix, e.g. `bipartite::Safariland`")

  # Convert matrix to data table
  net_dt <- as.data.table(net, keep.rownames = "row_names")

  # Melt from wide to long format
  net_dt <- melt(data = net_dt,
                 id.vars = "row_names",
                 variable.name = "col_names",
                 value.name = "counts")

  # Remove "fake" interactions
  net_dt <- net_dt[counts != 0]

  # Expand/explode by number of observed interactions if abund is TRUE
  if ( isTRUE(abund) ) {
    set.seed(seed = seed)
    net_dt <- net_dt[rep(1:.N, times = counts)]
  }

  # Shuffle row-wise
  set.seed(seed = seed)
  net_dt <- net_dt[sample(.N)]
  setDF(net_dt)

  return(net_dt)
}


#' @title Get interactions
#'
#' @description Helper function which constructs two vectors of species
#'   interactions (e.g. plant - pollinator) for two network interactions
#'   provided as data frames.
#'
#' @param net1,net2 Two data frames of interactions, for example constructed
#'   with \code{\link[bootdissim]{reshape_net}}. Each data frame must have only
#'   two columns/variables, e.g. first column for plant names, second column for
#'   insect names.
#'
#' @return Returns a list of two character vectors.
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
#' @export
#'
#' @md
get_interactions <- function(net1, net2){
  # Get the vector of interactions from each data frame (network)
  v1 <- paste(net1[[1]], net1[[2]])
  v2 <- paste(net2[[1]], net2[[2]])
  return(list(v1 = v1,
              v2 = v2))
}


#' @title Sample from two vectors
#'
#' @description Prepares a list of sample values from two given vectors `v1` &
#'   `v2` for the bootstrapping workflow. This function is important for
#'   internal use.
#'
#' @param v1,v2 Vectors of interactions from which to sample.
#'
#' @param size Integer giving the total sample size which, for example, can be
#'   the length of the shortest vector of interactions or of the longest one.
#'   Passed to \code{\link[base]{sample}}.
#'
#' @param by Integer. Number of interactions used to increase gradually the
#'   sampled vectors of interactions until all observations are sampled (given
#'   in the `size` argument).
#'
#' @param seed A single integer value. Passed to \code{\link[base]{set.seed}}.
#'   Set seed to get reproducible random results.
#'
#' @param replace Should sampling be with replacement? Passed to
#'   \code{\link[base]{sample}}.
#'
#' @return A list of two lists of character vectors containing interactions.
#'
#' @export
#'
#' @md
sample_vectors <- function(v1, v2, size, by, seed, replace){

  # Sample from character vector v1
  set.seed(seed)
  v1_spl <- sample(v1, size = size, replace = replace)
  # Split into a list of vectors of the approximate size/length given in "by".
  # Inspired from https://stackoverflow.com/a/16275428/5193830
  cuts_v1 <- cut(x = seq_along(v1_spl),
                 breaks = floor(length(v1_spl)/by),
                 labels = FALSE)
  v1_lst <- split(v1_spl, cuts_v1)
  # Combine the vectors incrementally in a separate list.
  v1_lst_incr <- vector(mode = "list", length = length(v1_lst))
  v1_lst_incr[[1]] <- v1_lst[[1]]
  for (i in 2:length(v1_lst_incr)){
    v1_lst_incr[[i]] <- c(v1_lst_incr[[i-1]], v1_lst[[i]])
  }

  # Sample from character vector v2
  set.seed(seed + 1)
  v2_spl <- sample(v2, size = size, replace = replace)
  # Split with the cuts defined for v1 above, so that we get the same sampling
  # approach. This consistency is important during the bootstrap and for
  # displaying results later on.
  v2_lst <- split(v2_spl, cuts_v1)
  # Combine the vectors incrementally in a separate list.
  v2_lst_incr <- vector(mode = "list", length = length(v2_lst))
  v2_lst_incr[[1]] <- v2_lst[[1]]
  for (i in 2:length(v2_lst_incr)){
    v2_lst_incr[[i]] <- c(v2_lst_incr[[i-1]], v2_lst[[i]])
  }

  return(list(v1 = v1_lst_incr,
              v2 = v2_lst_incr))
}
