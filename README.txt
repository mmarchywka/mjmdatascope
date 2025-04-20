
https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf

needed to bring from Beaver to Focal, ontop of prior
glib etc, 

 2014  . cdcpp datascope
 2015  ./run_datascope -compile
 2016  ls
 2017  apt-cache search libavcodec
 2018  sudo apt-get install libavcodec-extra5B
 2019  sudo apt-get install libavcodec-extra58
 2020  ./run_datascope -compile
 2021  apt-cache search libavcodec
 2022  sudo apt-get install libavcodec-dev
 2023  ./run_datascope -compile
 2024  apt-cache search ibswscale
 2025  sudo apt-get install ibswscale-dev ibswscale5
 2026  sudo apt-get install libswscale-dev libswscale5
 2027  ./run_datascope -compile
 2028  find -mtime -1
 2029  ./datascope.out

possible use in a headless mode for automation, 

https://stackoverflow.com/questions/6281998/can-i-run-glu-opengl-on-a-headless-server

Xvfb :5 -screen 0 800x600x24 &
export DISPLAY=:5
glxgears 
