#' @title Compute turnover components.
#'
#' @description Compute the a, b, c turnover components. Given two vectors, v1
#'   and v2, of species names or interactions, component "a" is the number of
#'   elements (species or interactions) present in both vectors v1 and v2, "b"
#'   is the number of elements present only in v1 and "c" only in v2. See the
#'   section "Beta diversity estimates" in  Novotny, V. (2009). Beta diversity
#'   of plant–insect food webs in tropical forests: a conceptual framework.
#'   Insect Conservation and Diversity, 2(1), 5-9.
#'
#' @param vect A list of two vectors with species names or interactions.
#'
#' @return A list with the 3 components a, b, c. The elements of the list are
#'   numeric vectors of length 1.
#'
#' @examples
#'
#' library(bootdissim)
#' library(bipartite)
#'
#' net1 <- reshape_net(vazarr, seed = 1)
#' net2 <- reshape_net(vazcer, seed = 1)
#' vect <- get_interactions(net1, net2)
#' abc <- get_abc(vect)
#'
#' @export
#'
#' @md
get_abc <- function(vect){
  v1 <- vect[[1]]
  v2 <- vect[[2]]

  # Fast testing if elements of a vector are present in a second vector. The
  # index position doesn't matter. What matters is if for a given element from
  # the first vector, the match() function returns a value of zero or above
  # zero. A zero indicates that the respective element from the first vector
  # doesn't exist in the second vector. A positive value indicates that the
  # respective element from the first vector exists in the second vector. The
  # value itself indicates the first position where that element is encountered
  # in the second vector, but this information is not important for our
  # purposes. We used match() function because is fast.
  match_v1_in_v2 <- match(v1, v2, nomatch = 0)
  match_v2_in_v1 <- match(v2, v1, nomatch = 0)

  # Get the a, b and c components.

  # "a" is the number of elements present in both vectors v1 and v2
  a <- sum(match_v2_in_v1 > 0) # same as sum(match_v1_in_v2 > 0)
  # "b" is the number of elements present only in v1
  b <- sum(!(match_v1_in_v2 > 0))
  # "c" is the number of elements present only in v2
  c <- sum(!(match_v2_in_v1 > 0))

  return(list(a = a, b = b, c = c))
}


#' @title Compute dissimilarity index.
#'
#' @description Compute various dissimilarity indices.
#'
#' @param abc A list with the turnover components a, b and c.
#'
#' @param index The name of the desired index, e.g.: "whittaker", "jaccard",
#'   "sorensen". Character type.
#'
#' @return A dissimilarity index. Numeric vector of length one.
#'
#' @examples
#'
#' library(bootdissim)
#' library(bipartite)
#'
#' net1 <- reshape_net(vazarr, seed = 1)
#' net2 <- reshape_net(vazcer, seed = 1)
#' vect <- get_interactions(net1, net2)
#' abc <- get_abc(vect)
#' get_beta(abc, index = "whittaker")
#' get_beta(abc, index = "jaccard")
#' get_beta(abc, index = "sorensen")
#'
#' @export
#'
#' @md

get_beta <- function(abc, index){
  a <- abc[["a"]]
  b <- abc[["b"]]
  c <- abc[["c"]]

  switch(index,
         whittaker = 2*(a+b+c)/(2*a+b+c) - 1, # sometimes written as ( (a + b + c) / ((2*a + b + c)/2) ) - 1
         jaccard = a/(a+b+c),
         sorensen = (2*a)/(2*a+b+c),
         stop("Your index is not recognised or available. Typo? Check help for options.",
              call. = FALSE))
  # Switch is generally faster than if statements.
  # https://stackoverflow.com/a/7826352/5193830
}

