# From factor to ordinal to make it comparables


# -------------------------------------------------------------------------

select.confounders.from.data = function(data, confounders){
  data_confounders = data[, confounders]
  return(data_confounders)
}

# -------------------------------------------------------------------------

check.ord.tf = function(c){
  tf = is.factor(c) & !is.ordered(c)
  return(tf)
}

# -------------------------------------------------------------------------

check.ord = function(data_confounders){
  check_ord = lapply(data_confounders, check.ord.tf)
  return(check_ord)
}

# -------------------------------------------------------------------------


check.ord.col = function(data_confounders, check_ord){
  col_to_order = colnames(data_confounders)[unlist(check_ord)]
  return(col_to_order)
}

# -------------------------------------------------------------------------

comparable.all = function(data_confounders, col_to_order){
  comparable_confounders = data_confounders
  comparable_confounders[col_to_order] = lapply(comparable_confounders[col_to_order], factor, ordered = T)
  return(comparable_confounders)
}

# -------------------------------------------------------------------------

comparables = function(data, confounders){

  data_confounders = select.confounders.from.data(data, confounders)

  check_ord = check.ord(data_confounders)

  col_to_order = check.ord.col(data_confounders, check_ord)

  comparable_confounders = comparable.all(data_confounders, col_to_order)

  return(comparable_confounders)
}




