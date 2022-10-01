

# install required packages
requirements <- c('fpp2', 'tseries', 'TSA', 'forecast', 'seastests', 'forecast', 
                  'plotly', 'prettydoc', 'stringi', 'stringr', 'polynom', 
                  'parallel')
installed <- .packages(all.available = TRUE)

for (package in requirements) {
    if (!(package %in% installed)) {
        install.packages(package)
    }
}

# create a miniconda environment to save static plotly figures
install.packages('reticulate')
reticulate::install_miniconda()
reticulate::conda_install('r-reticulate', 'python-kaleido==0.1.0')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')