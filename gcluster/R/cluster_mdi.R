library(matlabr)
# need to change the hard written filename
#mdi_func('myeloma_data/polya_100.csv','myeloma_data/ribosome_100.csv','myeloma_data/te_100.csv', 'Gaussian'); 

cluster_MDI <- function(datasets, dataType,
                        matlab_path = "~/Document/MATLAB/MDI") {
  
  add_path <- paste("addpath(genpath('", matlab_path , "'))", sep = "")
  script <- paste(datasets,collapse = ",")
  mdi_function <- paste("mdi_func(",script,",",dataType,")",sep = "")
  code <- c(add_path, mdi_function)
  run_matlab_code(code)
}

#datasets <- c("myeloma_data/polya_100.csv","myeloma_data/ribosome_100.csv","myeloma_data/te_100.csv")

#dataset_1 <- 'myeloma_data/polya_100.csv'
#dataset_2 <- 'myeloma_data/ribosome_100.csv'
#dataset_3 <- 'myeloma_data/te_100.csv'
#dataType <- 'Gaussian'
#mdi_func(dataset1, dataset_2, dataset_3,dataType)
#cluster_MDI(datasets, dataType)
