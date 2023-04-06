#!/usr/bin/env nextflow_multiomics
nextflow.enable.dsl=2

//nextflow script to run OmicsFold multi-omics integration pipeline

//tuning #1: tune optimal number of components
process findncomp {
  
  //use OmicsFold conda environment
  conda "~/.conda/envs/OmicsFold"
  
  //set output directory
  publishDir "omicsfold_output"
  
  input:
  path data
  path data_labels
  
  output:
  path data
  path data_labels
  path "1_find_ncomp_output.txt"
  path "1_find_ncomp_plots.pdf"
  path "1_find_ncomp_perf_result.rds"
  
  """
  Rscript 1_find_ncomp.R $data $data_labels
  """
}

//tuning #2: tune optimal number of features per component
process findkeepX {
  
  //use OmicsFold conda environment
  conda "~/.conda/envs/OmicsFold"
  
  //set output directory
  publishDir "omicsfold_output"
  
  input:
  path data
  path data_labels
  path perf_untrained
 
  
  output:
  path data
  path data_labels
  path "2_find_keepX_output.txt"
  path "2_find_keepX_plots.pdf"
  path "2_find_keepX_tune_result.rds"
  
  """
  Rscript 2_find_keepx.R $data $data_labels $perf_untrained
  """
}

//final performance check
 process perf_check {
   
  //use OmicsFold conda environment
  conda "~/.conda/envs/OmicsFold"
   
  //set output directory
  publishDir "omicsfold_output"
  
  input:
  path data
  path data_labels
  path tune_diablo
 
  
  output:
  path "3_loadings_stability_all.csv"
  path "3_perf_check_output.txt"
  path "3_perf_check_plots.pdf"
  path "3_perf_result_object.rds"
  path "3_trained_model_object.rds"
  path "3_blockrank.pdf"

  """
  Rscript 3_perf_check.R $data $data_labels $tune_diablo
  """
}

//define workflow 
workflow {
findncomp(params.data, params.data_labels)
findkeepX(findncomp.out[0], findncomp.out[1], findncomp.out[4])
perf_check(findkeepX.out[0], findkeepX.out[1], findkeepX.out[4])
}
