# MARMoT balancing --------------------------------------------------------

balancing = function(data, treatment, AR, reference, caliper){

  data$AR = AR

  AR_sorted = bal.sort.AR(AR)

  tabella = bal.better.tab(AR_sorted, data, treatment)

  caliper = bal.caliper(data, caliper)

  ref = bal.ref.selector(tabella, reference)

  # Non match removing -----------------------------------------------------

  AR = AR_sorted
  no_match = c()

  for (trt in unique(data[, treatment])) {

    data_trt = bal.treatment.data(data, treatment, trt)

    zeroes_match = bal.zero.match(tabella, trt, AR)

    for(single_zero_match in zeroes_match){

      elegible = bal.elegible.data(single_zero_match, caliper, data_trt)
      no_match = bal.no.elegible(single_zero_match, elegible, no_match)

    }
  }

  no_match = unique(no_match)

  # Balancing -----------------------------------------------------------

  balanced_data = data[0,]

  for (trt in unique(data[, treatment])) {

    data_trt = bal.treatment.data(data, treatment, trt)

    # Zeroes ----------------------------------------------------------------

    zeroes = bal.zeroes(data_trt, AR, tabella, trt, ref, caliper, no_match, AR_sorted)

    # Exact match -----------------------------------------------------------

    exact = bal.exact(data_trt, AR, tabella, trt, ref, no_match)

    # Inexact match (but non zero) ------------------------------------------

    inexact = bal.inexact(data_trt, AR, tabella, trt, ref, no_match, AR_sorted)

    # Merge -----------------------------------------------------------------

    balanced_data = rbind(balanced_data, exact, inexact, zeroes)
  }

  balanced_data[, "AR"] = NULL
  return(balanced_data)
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

bal.better.tab.row = function(single_AR, data, treatment){
  tab_row = table(data[data[, "AR"] == single_AR, treatment])
  return(tab_row)
}

# -------------------------------------------------------------------------

bal.better.tab = function(AR_sorted, data, treatment){
  data[, treatment] = factor(data[, treatment], levels = unique(data[, treatment]))
  tab_list = lapply(AR_sorted, bal.better.tab.row, data, treatment)
  tab = tab_list[[1]]
  for (i in 2:length(tab_list)) {
    tab = rbind(tab, tab_list[[i]])
  }
  rownames(tab) = NULL
  return(tab)
}

# -------------------------------------------------------------------------

bal.caliper = function(data, caliper){
  caliper = stats::sd(data[, "AR"])*caliper
  return(caliper)
}

# -------------------------------------------------------------------------

bal.median = function(tabella){
  median = ifelse(ceiling(apply(tabella, 1, median)) == 0, 1,
                  ceiling(apply(tabella, 1, median)))
  return(median)
}

# -------------------------------------------------------------------------

bal.mean = function(tabella){
  mean = ceiling(apply(tabella, 1, mean))
  return(mean)
}

# -------------------------------------------------------------------------

bal.mean0 = function(tabella){
  mean = floor(apply(tabella, 1, mean))
  return(mean)
}

# -------------------------------------------------------------------------

bal.median0 = function(tabella){
  median = ceiling(apply(tabella, 1, median))
  return(median)
}

# -------------------------------------------------------------------------

bal.min = function(tabella){
  median = ifelse(ceiling(apply(tabella, 1, min)) == 0, 1,
                  ceiling(apply(tabella, 1, min)))
  return(median)
}

# -------------------------------------------------------------------------

bal.max = function(tabella){
  median = ceiling(apply(tabella, 1, max))
  return(median)
}

# -------------------------------------------------------------------------

#alternative statistics

# -------------------------------------------------------------------------

bal.ref.selector = function(tabella, reference){
  if(reference == "median"){ref = bal.median(tabella)}
  if(reference == "mean"){ref = bal.mean(tabella)}
  if(reference == "mean0"){ref = bal.mean0(tabella)}
  if(reference == "median0"){ref = bal.median0(tabella)}
  if(reference == "min"){ref = bal.min(tabella)}
  if(reference == "max"){ref = bal.max(tabella)}
  return(ref)
}

# -------------------------------------------------------------------------

bal.sort.AR = function(AR){
  AR = sort(unique(AR))
  return(AR)
}

# -------------------------------------------------------------------------

bal.treatment.data = function(data, treatment, trt){
  data_trt = data[data[, treatment] == trt, ]
  return(data_trt)
}

# -------------------------------------------------------------------------

bal.zero.match = function(tabella, trt, AR){
  zeroes_match = AR[tabella[, trt] == 0]
  return(zeroes_match)
}

# -------------------------------------------------------------------------

bal.single.ref = function(ref, AR_sorted, single_x_match){
  single_ref = as.numeric(ref[AR_sorted == single_x_match])
  return(single_ref)
}

# -------------------------------------------------------------------------

bal.AR.ceiling = function(single_zero_match, caliper){
  top_AR = single_zero_match + caliper
  return(top_AR)
}

# -------------------------------------------------------------------------

bal.AR.floor = function(single_zero_match, caliper){
  bottom_AR = single_zero_match - caliper
  return(bottom_AR)
}

# -------------------------------------------------------------------------

bal.closest.elegible = function(elegible, single_zero_match){
  AR_elegible = unique(elegible[, "AR"])
  closest_elegible_AR = AR_elegible[which.min(abs(AR_elegible -
                                                    single_zero_match))]
  closest_elegible = elegible[elegible[, "AR"] == closest_elegible_AR, ]
  return(closest_elegible)
}

# -------------------------------------------------------------------------

bal.elegible.substitutes = function(data_trt, top_AR, bottom_AR){
  elegible = data_trt[data_trt[, "AR"] <= top_AR &
                        data_trt[, "AR"] >= bottom_AR,]
  return(elegible)
}

# -------------------------------------------------------------------------

bal.elegible.data = function(single_zero_match, caliper, data_trt){
  top_AR = bal.AR.ceiling(single_zero_match, caliper)
  bottom_AR = bal.AR.floor(single_zero_match, caliper)
  elegible = bal.elegible.substitutes(data_trt, top_AR, bottom_AR)
  return(elegible)
}

# -------------------------------------------------------------------------

bal.exact.match = function(AR, tabella, trt, ref){
  exact_match = AR[tabella[, trt] == ref]
  return(exact_match)
}

# -------------------------------------------------------------------------

bal.exact.data = function(data_trt, exact_match){
  exact = data_trt[data_trt[, "AR"] %in% exact_match,]
  return(exact)
}

# -------------------------------------------------------------------------

bal.exact = function(data_trt, AR, tabella, trt, ref, no_match){
  exact_match = bal.exact.match(AR, tabella, trt, ref)
  exact_match = bal.rm.no.match(exact_match, no_match)
  exact = bal.exact.data(data_trt, exact_match)
  return(exact)
}

# -------------------------------------------------------------------------

bal.inexact.match = function(AR, tabella, trt, ref){
  inexact_match = AR[tabella[, trt] != ref & tabella[, trt] != 0]
  return(inexact_match)
}

# -------------------------------------------------------------------------

bal.single.inexact.data = function(data_trt, single_inexact_match){
  single_inexact_data = data_trt[data_trt[, "AR"] == single_inexact_match,]
  return(single_inexact_data)
}

# -------------------------------------------------------------------------

bal.single.inexact.sampling = function(single_inexact_data, single_ref){
  single_inexact = single_inexact_data[sample(nrow(single_inexact_data),
                                              single_ref, replace = T), ]
  return(single_inexact)
}

# -------------------------------------------------------------------------

bal.inexact = function(data_trt, AR, tabella, trt, ref, no_match, AR_sorted){
  inexact_match = bal.inexact.match(AR, tabella, trt, ref)
  inexact = data_trt[0,]

  for (single_inexact_match in inexact_match) {

    if(single_inexact_match %in% no_match){next}
    single_ref = bal.single.ref(ref, AR_sorted, single_inexact_match)
    single_inexact_data = bal.single.inexact.data(data_trt,
                                                  single_inexact_match)
    single_inexact = bal.single.inexact.sampling(single_inexact_data,
                                                 single_ref)
    inexact = rbind(inexact, single_inexact)
  }
  return(inexact)
}

# -------------------------------------------------------------------------

bal.no.elegible = function(single_zero_match, elegible, no_match){
  if(nrow(elegible) == 0){
    no_match = c(no_match, single_zero_match)
  }
  return(no_match)
}

# -------------------------------------------------------------------------

bal.rm.no.match = function(exact_match, no_match){
  exact_match = exact_match[!exact_match %in% no_match]
  return(exact_match)
}

# -------------------------------------------------------------------------

bal.zeroes = function(data_trt, AR, tabella, trt, ref, caliper, no_match, AR_sorted){
  zeroes_match = bal.zero.match(tabella, trt, AR)
  zeroes = data_trt[0, ]

  for(single_zero_match in zeroes_match){
    if(single_zero_match %in% no_match){next}
    single_ref = bal.single.ref(ref, AR_sorted, single_zero_match)
    elegible = bal.elegible.data(single_zero_match, caliper, data_trt)
    elegible = bal.closest.elegible(elegible, single_zero_match)
    single_zeroes = elegible[sample(nrow(elegible),
                                    single_ref, replace = T), ]
    zeroes = rbind(zeroes, single_zeroes)
  }
  return(zeroes)
}


