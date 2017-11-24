pkgs <- readLines(file("https://raw.githubusercontent.com/rbslandau/statistics_multi/master/Data/installed_pkgs.txt", "r"))
str(pkgs)
install.packages(pkgs)