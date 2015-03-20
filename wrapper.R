
###################### SIMULATION WRAPPER ######################

# this file controls the entire simulation process: 
# 1. generates dataset with Xs
# 2. generates Ys
# 3. collapse to "observed" dataset
# 4. splits on follow-up times
# 5. fits right and wrong models
# 6. compare the 2 models

# args in sbatch file for genSurv
-- arg 1: source location
-- arg 2: beta location
-- arg 3: covariate data path
-- arg 4: number cores
-- arg 5: number segments
-- arg 6: segment index
-- arg 7: results write directory
-- arg 8: sim results name
-- arg 9: log write directory
-- arg 10: data write directory

# n = number subjects

simulation_wrapper = function() {
  
  # 
}
