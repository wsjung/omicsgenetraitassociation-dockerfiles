library(corrmeta)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

input_dir <- args[1]
#output_dir <- args[2]
#trait_ <- args[3]
#category_ <- args[4]

#input_dir <- file.path(input_parent_dir, trait_, category_)
#output_dir <- file.path(output_parent_dir, trait_, category_)

# create output directory if it doesn't exist
#dir.create(output_dir, showWarnings=FALSE)

input_files <- list.files(input_dir)
input_vars <- lapply(input_files, function(x) strsplit(x, "\\.")[[1]][1]) %>% unlist
input_files <- file.path(input_dir, input_files) # add back input_dir

print(paste0("Found ", length(input_vars), " input files."))

# read only relevant columns from files "markname,pval"
df_input_list <- lapply(input_files, function(x) select(read.csv(x), markname, pval) )
# merge by markname
df_input <- Reduce(function(x,y) merge(x,y, by='markname', all=TRUE), df_input_list)
# remove duplicates
df_input <- df_input[!duplicated(df_input),]
colnames(df_input)[-1] = input_vars
print(head(df_input))
print(paste0("# num rows: ", nrow(df_input)))


# ------------- corrmeta --------------- #
# calculate tetrachoric correlations
result_tetracorr <- corrmeta::tetracorr(df_input, input_vars)


# calcaulte fisher's p-value
result_fisher <- corrmeta::fishp(df_input, input_vars, result_tetracorr$sigma, result_tetracorr$sum_sigma)
print(paste0("# fisher result nrows: ", nrow(result_fisher)))
print(head(result_fisher))

# merge final output file
result_fisher$corr <- (result_fisher$sum_sigma_var - 2) / 2

# calculate nlog10p for input columns
result_fisher <- result_fisher %>%
    mutate(across(all_of(input_vars), ~ -log10(.), .names = "{.col}_nlogp"))

# subset final output file columns
corrmeta_output <- result_fisher %>%
    select(c("markname","meta_nlog10p","meta_p"), contains(input_vars, ignore.case=TRUE), "corr") %>%
    arrange(meta_p)
print(paste0("# corrmeta_output result nrows: ", nrow(corrmeta_output)))
print(head(corrmeta_output))

# write tetrachoric correlations to disk
write.table(result_tetracorr$sigma, "tetrachor_sigma.txt", quote=FALSE, sep="\t", row.names=FALSE, na="NA")
# write corrmeta result to disk
write.csv(corrmeta_output, "CMA_meta.csv", quote=FALSE, row.names=FALSE, na="NA")

print("DONE")
