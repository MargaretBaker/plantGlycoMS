#!/usr/bin/env bash
Rscript -e "roxygen2::roxygenise('.');"

R CMD build .
R CMD check plantGlycoMS_0.1.tar.gz --as-cran
# R CMD Rd2pdf . -o "manual.pdf" --no-preview --force

# Rscript -e "rmarkdown::render('vignettes/vignette.Rmd', 'all');"

git add .
git commit -am 'init'
git push origin master
