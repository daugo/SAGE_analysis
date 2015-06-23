#understanding the ocurrences of restriction sites in the potential 3' UTR part

freq <- read.csv("./phytophthora_data/count_of_cutting_sequence_in_extended_part.csv")
names(freq) <- c("id", "added_nucleotides", "number_restr_sites")
hist(freq$number_restr_sites)
hist(freq$added_nucleotides)
plot(freq$added_nucleotides, freq$number_restr_sites)
plot(freq$number_restr_sites, freq$added_nucleotides)


hist(freq$number_restr_sites[freq$added_nucleotides>980])
