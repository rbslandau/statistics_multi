pkgs <- readLines(file("https://www.uni-koblenz-landau.de/en/campus-landau/faculty7/environmental-sciences/landscape-ecology/Teaching/inst_pkgs/at_download/file", "r"))
str(pkgs)
install.packages(pkgs)