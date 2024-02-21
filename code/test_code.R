library(reticulate)

use_python("C:/Py27/python.exe")
# py_install(envname = "C:/Users/peterman.73/AppData/Local/Continuum/miniconda2/envs/",c("numpy", 'scipy'))
py_config()


import('numpy')
import('scipy')

test.py <- py_run_string("def my_function(): print('Hello from a function')")

test.py()
setwd("C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/")


py_run_string("cd 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/")

py_run_string('C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/
              inputvars.csv output_test')



# Make Python function ----------------------------------------------------
system(paste("python C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/ inputvars.csv output_test"))

CDPOP_run <- function(CDPOP.py,
                      data.dir,
                      input,
                      output){
  
  system(paste("python", CDPOP.py, data.dir, input, output, sep = " "))
  
}


CDPOP_run(CDPOP.py = 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py',
          data.dir = 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/data/',
          input = 'inputvars.csv',
          output = 'output')

cdpop <- function(CDPOP.py,
                  data.dir,
                  input,
                  output) {
  
}

def CDPOP():
  script = sys.argv[1]
filename = sys.argv[1]
data = np.loadtxt(filename, delimiter=',')
for m in data.mean(axis=1):
  print m


# Example -----------------------------------------------------------------
# create a new environment 
conda_create("r-reticulate")

# install SciPy
conda_install("r-reticulate", "scipy", "numpy")

# import SciPy (it will be automatically discovered in "r-reticulate")
scipy <- import("scipy")
numpy <- import("numpy")

# https://rstudio-pubs-static.s3.amazonaws.com/407460_396f867ce3494d479fd700960879e22c.html
py_config()

py_run_file()