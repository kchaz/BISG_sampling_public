# This code is meant to be run only once to generate a list
# of fake surnames that will then be used across all simulated data

library(glue)
set.seed(25)

k = 10000

surnames <- replicate(k,
                      paste0(sample(letters,
                                    6,
                                    replace = T),
                                    collapse = ""))
surnames = unique(surnames)
stopifnot(length(surnames) == k)

save_path <- paste0("../../data/simulated/")
if (!dir.exists(save_path)){
  dir.create(save_path)
}

saveRDS(surnames[1:10000], glue("{save_path}surnames{k}.RDS"))
