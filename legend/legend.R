#######################################################
# Example key table, and legend
#######################################################

#######################################################
# Make a table that is the abbreviations key
#######################################################
# Here, the file with the "key" is a tab-delimited text file -- 
# change to your file location if you don't want Hawaii
keyfn = "./abbreviations_table.txt"

keydf = read.table(keyfn, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
keydf

keydf2 = keydf[,c("ab1", "desc")]
names(keydf2) = c("Abbr", "Description")
keydf2

library(gridExtra)    # for grid.table() function
grid.table(keydf2, rows=NULL)

# Defaults, if you want to use the Hawaii names/abbreviations from above
max_range_size = 8
include_null_range = TRUE
areanames = keydf$ab1
areanames

# You could paste the below directly into your script, if it is
# set up like in the PhyloWiki example script

#######################################################
# Plot legend for ALL states/ranges (there may be a 
# ton, and getting them all to display is hard)
#######################################################
pdffn = "colors_legend_all_v1.pdf"
pdf(pdffn, width=40, height=8)

areanames = names(tipranges@df)
areanames

include_null_range = TRUE

library(cladoRcpp)
states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range)

statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range, split_ABC=FALSE)
statenames

relprobs_matrix = res$ML_marginal_prob_each_state_at_branch_top_AT_node
MLprobs = get_ML_probs(relprobs_matrix)
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")

colors_matrix = get_colors_for_numareas(length(areanames))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=include_null_range)
colors_list_for_states

possible_ranges_list_txt = areas_list_to_states_list_new(areas=areanames,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

legend_ncol=NULL
legend_cex=1.5
colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)

dev.off()

#######################################################
# Plot legend for SOME states/ranges 
# (especially since the widespread ranges are just
#  combinations of the primary colors used for 
#  single-area ranges)
#######################################################

pdffn = "colors_legend_some_v1.pdf"
pdf(pdffn, width=6, height=6)

# Subset to just some ranges (since there are sooo many combinations)
states_to_put_in_legend = c(1,2,3,4,5,6,7,8)
colors_list_for_states_subset = colors_list_for_states[states_to_put_in_legend]
possible_ranges_list_txt_subset = possible_ranges_list_txt[states_to_put_in_legend]

legend_ncol=NULL
legend_cex=1.5
colors_legend(possible_ranges_list_txt_subset, colors_list_for_states_subset, legend_ncol=legend_ncol, legend_cex=legend_cex)

dev.off()

