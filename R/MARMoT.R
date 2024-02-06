#' MARMoT balancing method
#'
#' @description
#' Matching on poset-based average rank for multiple treatments (MARMoT).
#'
#' @details
#' There are many scenarios where classic propensity score techniques are not
#' applicable (e.g. there are many treatments). In a multiple-treatment
#' framework, MARMoT is a method to balance the distribution of covariates among
#' several treatment groups. MARMoT introduces a method for achieving balance
#' among treatment groups by utilizing partially ordered set (poset) theory.
#' This approach focuses on equalizing individual characteristics without
#' relying on propensity score techniques or a dependent variable. Unlike
#' propensity score methods, poset theory doesn't require assumptions about
#' model specifications for treatment allocation. Each subject is represented by
#' a profile of their characteristics, and an average rank approximation is
#' associated with each profile. This value represents the significance of
#' individual characteristics for treatment allocation and can be normalized for
#' better interpretability.
#'
#' @param data
#' A dataframe or equivalent.
#' @param confounders
#' A vector containing the column names of the confounders to balance by.
#' @param treatment
#' A string indicating the column name of the treatment variable.
#' @param reference
#' The statistic used to determine the reference frequencies in the balancing
#' process. Default is median.
#' @param n.cores
#' Number of cores to be used (Linux and Mac systems only!); if a
#' number grater than 1 is specified the function will use a parallelized
#' version of the deloof approximation. Default set to 1.
#' @param caliper
#' Fraction of the standard deviation used to determine the closest neighbour.
#' Default is 0.25.
#' @param verbose
#' Set to FALSE to suppress any console output. Default is TRUE
#'
#' @return
#' A list of objects, also containing the balanced dataset with the same
#' structure of the input dataset.
#' @export
#'
#' @references
#' Silan, M., Boccuzzo, G. and Arpino, B., 2021. 'Matching on poset‐based average
#' rank for multiple treatments to compare many unbalanced groups'. Statistics in
#' Medicine, 40(28), pp.6443-6458. https://doi.org/10.1002/sim.9192
#'
#' Silan, M., Belloni, P. and Boccuzzo, G., 2023. 'Identification of
#' neighborhood clusters on data balanced by a poset-based approach'.
#' Statistical Methods & Applications, pp.1-22.
#' https://doi.org/10.1007/s10260-023-00695-0
#'
#' @examples
#' out = MARMoT(data = MARMoT_data, confounders = c("race", "age"),
#'              treatment = "hospital", n.cores = 1)
#' out
MARMoT = function(data, confounders, treatment, reference = "median",
                  n.cores = 1, caliper = 0.25, verbose = TRUE){



  # Start time --------------------------------------------------------------


  start_time = Sys.time()


  # Convert factor to ordinal -----------------------------------------------


  comparable_confounders = comparables(data = data, confounders = confounders)


  # ASB pre balancing -------------------------------------------------------


  if(verbose){print("Absolute standardized bias - before balancing")}
  ASB_pre = ASB(data = data, confounders = confounders,
                treatment = treatment, verbose = verbose)


  # Average rank with deloof approx -----------------------------------------


  if(verbose){print("... Computing average rank ...")}
  if(n.cores == 1){
    AR = deloof(comparable_confounders)}
  if(n.cores > 1){
    AR = mcdeloof(comparable_confounders, n.cores)}


  # Normalize AR ------------------------------------------------------------


  AR = normalize.AR(AR)


  # Balancing ---------------------------------------------------------------


  if(verbose){print("... Balancing ...")}
  balanced_data = balancing(data = data, treatment = treatment, AR = AR,
                            reference, caliper)


  # ASB post ----------------------------------------------------------------


  if(verbose){print("Absolute standardized bias - after balancing")}
  ASB_post = ASB(data = balanced_data, confounders = confounders,
                 treatment = treatment, verbose = verbose)



  # End time ----------------------------------------------------------------


  time = difftime(Sys.time(), start_time, units = "mins")


  # Output ------------------------------------------------------------------


  output = output.maker(balanced_data, ASB_pre, ASB_post, time)


  return(output)
}
