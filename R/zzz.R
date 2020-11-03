

.onLoad  <-  function(libname, pkgname)  {
  library.dynam("spcopula", pkgname, libname)
}

.onUnload  <- function(libpath)  {
  library.dynam.unload("spcopula", libpath)
}