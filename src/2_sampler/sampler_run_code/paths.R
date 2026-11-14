
args <- commandArgs(trailingOnly = TRUE) #true means doesn't include script name


# If this is not being run on the cluster
# use paths purely internal to the repo

if (length(args) == 0){
  src_path = "../../"
  data_path = "../../../data/"
  data_object_path = "../../../data/DataObjects/"
  obit_data_path = "../../../data/ObitData/"
  save_sampler_output_path = "obitdata_outputs/"
}

#If this is being run on the cluster, use paths fed from there
if (length(args) > 0){
  if (length(args) != 5){
    stop("Must specify 5 paths as arguments: order is src_path, data_path,
         data_object_path, obit_data_path, and save_sampler_output_path")
  }
  src_path = args[1]
  data_path = args[2]
  data_object_path = args[3]
  obit_data_path = args[4]
  save_sampler_output_path = args[5]

}
