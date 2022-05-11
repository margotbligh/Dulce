ppm_to_mz = function(mz, noise){
  ppm = mz / 1000000 * noise
  return(ppm)
}

calc_distances = function(x){
  df = reshape2::melt(as.matrix(dist(x@mz)), varnames=c("f1","f2"), value.name="mzDifference")
  df = df[as.numeric(df$f1)>as.numeric(df$f2),]
  df$f1 = x@mz[df$f1]
  df$f2 = x@mz[df$f2]
  return(df)
}
