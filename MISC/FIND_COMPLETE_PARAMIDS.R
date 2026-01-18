# # File naming: control_sterile_macspec_tregs_ros_level_pat_level_trnd
# # Dynamically create file lists for all ros/pat combinations
# all_files = list()
# all_indices = list()
# for (ros_in in ros_vals) {
#   for (pat_in in pat_vals) {
#     var_name = paste0("files_", ros_in, "_", pat_in)
#     pattern_str = paste0("^longitudinal_df_param_set_id_\\d+\\_sterile_",sterile_in,
#                          "_macspec_",macspec_in,
#                          "_tregs_",tregs_on_in,
#                          "_ros_level_",ros_in,
#                          "_pat_level_",pat_in,
#                          "_trnd_",trnd_in,".rds$")
#     files_temp = list.files(path_data, pattern = pattern_str, full.names = TRUE)
#     all_files[[var_name]] = files_temp
# 
#     # Extract indices
#     indices_name = paste0("indices_", ros_in, "_", pat_in)
#     all_indices[[indices_name]] = str_extract(basename(files_temp), "\\d+") |> as.numeric()
#   }
# }
# 
# # Find common indices across all combinations
# indices = Reduce(intersect, all_indices)