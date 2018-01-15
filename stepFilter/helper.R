# Helper function for string wrapping. 
# Default 20 character target width.
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)