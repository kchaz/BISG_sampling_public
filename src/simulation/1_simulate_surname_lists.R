# NOTE: this code is meant to be run only once
# to generate lists of fake surnames that will then
# be used across all simulated data
set.seed(25)
k = 10000

surnames <- replicate(k, paste0(sample(letters, 6, replace = T),
                                      collapse = ""))
stopifnot(length(unique(surnames)) == k)


#saving in different overlapping chunks
saveRDS(surnames[1:100], "simulated_data/surnames100.RDS")
saveRDS(surnames[1:500], "simulated_data/surnames500.RDS")
saveRDS(surnames[1:1000], "simulated_data/surnames1000.RDS")
saveRDS(surnames[1:2000], "simulated_data/surnames2000.RDS")
saveRDS(surnames[1:5000], "simulated_data/surnames5000.RDS")
saveRDS(surnames[1:10000], "simulated_data/surnames10000.RDS")


