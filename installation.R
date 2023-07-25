

# install required packages
requirements <- c('remotes', 'fpp2', 'tseries', 'TSA', 'forecast', 'seastests', 'forecast', 
                  'libcurl', 'httr', 'plotly', 'prettydoc', 'stringi', 'stringr', 
                  'polynom', 'parallel')
installed <- .packages(all.available = TRUE)

for (package in requirements) {
    if (!(package %in% installed)) {
        install.packages(package)
    }
}
remotes::install_github("rstudio/reticulate")

# create a miniconda environment to save static plotly figures
reticulate::install_miniconda(force=T)
reticulate::conda_install('r-reticulate', 'python-kaleido==0.1.0')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')
