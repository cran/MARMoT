normalize.AR = function(AR){
  AR = (AR - min(AR)) / (max(AR) - min(AR))
  return(AR)
}