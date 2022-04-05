ppm_to_mz = function(mz, noise){
  ppm = mz / 1000000 * noise
  return(ppm)
}
