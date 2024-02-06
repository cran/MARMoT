#' Average rank with Deloof's approximation (parallel computing)
#'
#' @description
#' Compute the average with using Deloof's approximation, using parallel
#' computing; it can be used only on Linux and Mac systems!
#'
#' @param comparable_data
#' A dataframe or equivalent, which columns are ordered factors or numeric
#' variables.
#' @param n.cores
#' Number of cores to be used.
#'
#' @return
#' A vector containing the average rank of all the observations (it is
#' recommended to normalized it before use).
#' @export
#'
#' @references
#' Caperna, G., 2019. Approximation of AverageRank by means of a formula.
#' https://doi.org/10.5281/zenodo.2565699
#'
#' @examples
#' AR = mcdeloof(deloof_data, n.cores = 1)
#'
mcdeloof = function(comparable_data, n.cores){

  if(deloof.check.os()){

    profiles = parsec::pop2prof(comparable_data, sep = "")

    incidence_matrix_zeta = mcgetzeta(profiles, n.cores)

    profile_list = deloof.profile.list(profiles)
    n_profiles = deloof.n.profiles(profile_list)

    incomparables = parsec::incomp(incidence_matrix_zeta)

    prob_matrix = deloof.mut.prob.matrix(n_profiles, profile_list)

    row_positions = deloof.lo.row.positions(n_profiles)
    col_positions = deloof.lo.col.positions(n_profiles)

    downsets = deloof.sets(incidence_matrix_zeta, 2)
    upsets = deloof.sets(incidence_matrix_zeta, 1)

    lo_prob = deloof.lo.prob(upsets, downsets, row_positions, col_positions)

    prob_matrix = deloof.prob.matrix(prob_matrix, lo_prob, incomparables)

    downsets2 = deloof.sets.update(downsets, 1, prob_matrix)
    upsets2 = deloof.sets.update(upsets, 2, prob_matrix)

    prob_matrix2 = deloof.mut.prob.matrix(n_profiles, profile_list)

    lo_prob2 = deloof.lo.prob(upsets2, downsets2, row_positions, col_positions)

    prob_matrix2 = deloof.prob.matrix(prob_matrix2, lo_prob2, incomparables)

    result = deloof.result(downsets, prob_matrix2)

    rank_list = deloof.rank.list(comparable_data, result, profile_list, n_profiles)

    return(rank_list)

    } else{warning("You can't use multiple cores on Windows (set n.cores = 1)")}
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

deloof.check.os = function(){
  os = Sys.info()["sysname"]
  if(os != "Windows"){
    no_win = TRUE
  } else {no_win = FALSE}
  return(no_win)
}

# -------------------------------------------------------------------------

mcsapply = function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, n.cores)
{
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ..., mc.cores = n.cores)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!isFALSE(simplify))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

# -------------------------------------------------------------------------

deloof.one.vs.all = function(b, a, y){
  out = all(y[a, ] <= y[b, ])
  return(out)
}

# -------------------------------------------------------------------------

deloof.all.vs.all = function(a, y, m){
  sapply(1:m , deloof.one.vs.all, a, y)
}

# -------------------------------------------------------------------------

mcgetzeta = function (y, n.cores){
  y = y$profiles
  m = nrow(y)
  nam = rownames(y)
  res = mcsapply(1:m, deloof.all.vs.all, y, m, n.cores = n.cores)
  res = t(res)
  dimnames(res) = list(nam, nam)
  class(res) <- "incidence"
  return(res)
}

# -------------------------------------------------------------------------

deloof.strings = function(comparable_data){
  strings = apply(X = comparable_data, MARGIN = 1, FUN = toString)
  return(strings)
}

# -------------------------------------------------------------------------

deloof.label = function(strings){
  labels = gsub(", ", "", strings)
  return(labels)
}

# -------------------------------------------------------------------------

deloof.empty.rank.list = function(labels){
  empty_rank_list = rep(NA, length(labels))
  return(empty_rank_list)
}

# -------------------------------------------------------------------------

deloof.profile.list = function(profiles){
  profile_list = names(profiles$freq)
  return(profile_list)
}

# -------------------------------------------------------------------------

deloof.n.profiles = function(profile_list){
  n_profiles = length(profile_list)
  return(n_profiles)
}

# -------------------------------------------------------------------------

deloof.mut.prob.matrix = function(n_profiles, profile_list){
  prob_matrix = matrix(0, n_profiles, n_profiles)
  colnames(prob_matrix) = profile_list
  rownames(prob_matrix) = profile_list
  return(prob_matrix)
}

# -------------------------------------------------------------------------

deloof.sets = function(incidence_matrix_zeta, rowcol){
  sets = apply(incidence_matrix_zeta, rowcol, function(x) sum(x) - 1)
  return(sets)
}

# -------------------------------------------------------------------------

deloof.lo.row.positions = function(n_profiles){
  row_matrix = matrix(1:n_profiles, nrow = n_profiles, ncol = n_profiles)
  lo_row_matrix = row_matrix[lower.tri(row_matrix)]
  row_positions = c(lo_row_matrix)
  return(row_positions)
}

# -------------------------------------------------------------------------

deloof.lo.col.positions = function(n_profiles){
  col_matrix = matrix(1:n_profiles, nrow = n_profiles, ncol = n_profiles, byrow = T)
  lo_col_matrix = col_matrix[lower.tri(col_matrix)]
  col_positions = c(lo_col_matrix)
  return(col_positions)
}

# -------------------------------------------------------------------------

deloof.lo.prob = function(upsets, downsets, row_positions, col_positions){
  a = downsets[row_positions] + 1
  b = upsets[col_positions] + 1
  c = downsets[col_positions] + 1
  d = upsets[row_positions] + 1
  lo_prob = (a * b) / ((a * b) + (c * d))
  return(lo_prob)
}


# -------------------------------------------------------------------------

deloof.lo.prob.matrix = function(prob_matrix, lo_prob){
  prob_matrix[lower.tri(prob_matrix)] = lo_prob
  return(prob_matrix)
}

# -------------------------------------------------------------------------

deloof.uplo.prob.matrix = function(prob_matrix, lo_prob){
  prob_matrix[upper.tri(prob_matrix)] = t(1 - prob_matrix)[upper.tri(prob_matrix)]
  return(prob_matrix)
}

# -------------------------------------------------------------------------

deloof.nocomp.prob.matrix = function(prob_matrix, incomparables){
  prob_matrix[!incomparables] = 0
  return(prob_matrix)
}

# -------------------------------------------------------------------------

deloof.prob.matrix = function(prob_matrix, lo_prob, incomparables){
  prob_matrix = deloof.lo.prob.matrix(prob_matrix, lo_prob)
  prob_matrix = deloof.uplo.prob.matrix(prob_matrix, lo_prob)
  prob_matrix = deloof.nocomp.prob.matrix(prob_matrix, incomparables)
  return(prob_matrix)
}

# -------------------------------------------------------------------------

deloof.sets.update = function(set, rowcol, prob_matrix){
  set2 = set + apply(prob_matrix, rowcol, sum)
  return(set2)
}

# -------------------------------------------------------------------------

deloof.result = function(downsets, prob_matrix2){
  result = downsets + 1 + apply(prob_matrix2, 1, sum)
  return(result)
}

# -------------------------------------------------------------------------

deloof.rank.list.update = function(result, profile_list, n_profiles,
                                   empty_rank_list, labels){
  for (i in 1: n_profiles){
    empty_rank_list[which(labels == profile_list[i])] = result[i]
  }
  return(empty_rank_list)
}

# -------------------------------------------------------------------------

deloof.rank.list = function(comparable_data, result, profile_list, n_profiles){
  strings = deloof.strings(comparable_data)
  labels = deloof.label(strings)
  empty_rank_list = deloof.empty.rank.list(labels)
  rank_list = deloof.rank.list.update(result, profile_list, n_profiles,
                                      empty_rank_list, labels)
  return(rank_list)
}

