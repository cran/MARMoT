#' Absolute standardized bias
#'
#' @description
#' Compute the absolute standardized bias of given confounders and return some
#' useful statistics.
#'
#' @param data
#' A dataframe or equivalent.
#' @param confounders
#' A vector with the column names of the confounders to balance by.
#' @param treatment
#' A string with the column name of the treatment variable.
#' @param verbose
#' Set to FALSE to suppress any console output. Default is TRUE
#'
#' @return
#' A list of objects, containing the ASB matrix and some summary statistics.
#' @export
#'
#' @examples
#' ASB(data = MARMoT_data, confounders = c("race", "age"), treatment = "hospital")
#'
ASB = function(data, confounders, treatment, verbose = TRUE){

  tab_list = ASB.tab.list(data, confounders, treatment)
  tab = ASB.chain.matrix.list(tab_list)

  sel_vector = sel.vector(confounders, data)

  new_rownames_list = ASB.new.rownames.list(confounders, data)
  new_rownames = ASB.chain.vector.list(new_rownames_list)
  rownames(tab) = paste0(new_rownames, sep = "_", rownames(tab))

  tab = ASB.remove.one.level(tab, sel_vector)

  n_x_trt = ASB.n.x.trt(data, treatment)
  tab_perc = ASB.tab.perc(tab, n_x_trt)
  var_trt = ASB.var.trt(tab_perc)
  tab_tot = ASB.tab.tot(tab)
  tab_tot_perc = ASB.tab.tot.perc(tab_tot, n_x_trt)
  var_tot = ASB.var.tot(tab_tot_perc)

  numerator = ASB.numerator(tab_perc, tab_tot_perc)
  denominator = ASB.denominator(var_trt, var_tot)

  ASB_result = ASB.results.matrix(numerator, denominator)

  ASB_statistics = ASB.statistics(ASB_result)
  if(verbose){print(round(ASB_statistics, digits = 3))}

  output = ASB.out(ASB_result, ASB_statistics)

  return(output)
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

ASB.tab.single = function(x, data, treatment){
  tab_single = table(data[, x], data[, treatment])
  return(tab_single)
}

# -------------------------------------------------------------------------

ASB.tab.list = function(data, confounders, treatment){
  tab_list = lapply(confounders, ASB.tab.single, data, treatment)
  return(tab_list)
}

# -------------------------------------------------------------------------

ASB.chain.matrix.list = function(lst){
  tab = lst[[1]]
  for (j in 2:length(lst)) {
    tab = rbind(tab, lst[[j]])
  }
  return(tab)
}

# -------------------------------------------------------------------------

ASB.rm.one.level.single = function(c, data){
  lv = unique(data[, c])
  n_lv = length(lv)
  if(n_lv == 2){
    sel_vector = c(1, 0)
  } else{
    sel_vector = rep(1, length(lv))
  }
  return(sel_vector)
}

# -------------------------------------------------------------------------

ASB.rm.one.level.list = function(confounders, data){
  sel_vector_list = lapply(confounders, ASB.rm.one.level.single, data)
  return(sel_vector_list)
}

# -------------------------------------------------------------------------

ASB.chain.vector.list = function(lst){
  tab = lst[[1]]
  for (j in 2:length(lst)) {
    tab = c(tab, lst[[j]])
  }
  return(tab)
}

# -------------------------------------------------------------------------

sel.vector = function(confounders, data){
  sel_vector_list = ASB.rm.one.level.list(confounders, data)
  sel_vector = ASB.chain.vector.list(sel_vector_list)
  sel_vector = sel_vector*c(1:length(sel_vector))
  return(sel_vector)
}

# -------------------------------------------------------------------------

ASB.new.rownames.single = function(c, data){
  lv = levels(data[, c])
  n_lv = length(lv)
  name_vec = rep(c, n_lv)
  return(name_vec)
}

# -------------------------------------------------------------------------

ASB.new.rownames.list = function(confounders, data){
  names_list = lapply(confounders, ASB.new.rownames.single, data)
  return(names_list)
}

# -------------------------------------------------------------------------

ASB.remove.one.level = function(tab, sel_vector){
  tab = tab[sel_vector, ]
  return(tab)
}

# -------------------------------------------------------------------------

ASB.n.x.trt = function(data, treatment){
  n_x_trt = c(table(data[, treatment]))
  return(n_x_trt)
}

# -------------------------------------------------------------------------

ASB.tab.perc = function(tab, n_x_trt){
  tab_perc = t(t(tab) / n_x_trt)
  return(tab_perc)
}

# -------------------------------------------------------------------------

ASB.var.trt = function(tab_perc){
  var_trt = tab_perc*(1-tab_perc)
  return(var_trt)
}

# -------------------------------------------------------------------------

ASB.tab.tot = function(tab){
  tab_tot = rowSums(tab)
  return(tab_tot)
}

# -------------------------------------------------------------------------

ASB.tab.tot.perc = function(tab_tot, n_x_trt){
  tab_tot_perc = tab_tot / sum(n_x_trt)
  return(tab_tot_perc)
}

# -------------------------------------------------------------------------

ASB.var.tot = function(tab_tot_perc){
  var_tot = tab_tot_perc*(1-tab_tot_perc)
  return(var_tot)
}

# -------------------------------------------------------------------------

ASB.numerator = function(tab_perc, tab_tot_perc){
  num = abs(tab_perc - tab_tot_perc)
  return(num)
}

# -------------------------------------------------------------------------

ASB.denominator = function(var_trt, var_tot){
  den = (sqrt((var_trt + var_tot) / 2))
  return(den)
}

# -------------------------------------------------------------------------

ASB.results.matrix = function(num, den){
  ASB = num / den * 100
  return(ASB)
}

# -------------------------------------------------------------------------

ASB.quantiles = function(ASB_result){
  quant = stats::quantile(ASB_result, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(quant) = c("Min", "1st quartile", "Median", "3rd quartile", "Max")
  return(quant)
}

# -------------------------------------------------------------------------

ASB.mean = function(ASB_result){
  mn = mean(ASB_result)
  names(mn) = "Mean"
  return(mn)
}

# -------------------------------------------------------------------------

ASB.over5 = function(ASB_result){
  over5 = length(which(ASB_result>5))
  names(over5) = "Over 5%"
  return(over5)
}

# -------------------------------------------------------------------------

ASB.over10 = function(ASB_result){
  over10 = length(which(ASB_result>10))
  names(over10) = "Over 10%"
  return(over10)
}

# -------------------------------------------------------------------------

ASB.statistics = function(ASB_result){
  quant = ASB.quantiles(ASB_result)
  mn = ASB.mean(ASB_result)
  over5 = ASB.over5(ASB_result)
  over10 = ASB.over10(ASB_result)
  ASB_statistics = c(quant, mn, over5, over10)
  return(ASB_statistics)
}

# -------------------------------------------------------------------------

ASB.out = function(ASB_result, ASB_statistics){
  output = list("Matrix" =  ASB_result, "Stat" = ASB_statistics)
  return(output)
}


