library(ape)

# script to convert Nexus files to newick, for us in AU testing with iqtree
for (geno in c("B3", "D8")) {
  # Whole-genome timetree as nexus, convert to newick format
  tr_wgs <- read.nexus(file.path(geno, paste0(geno, "_timetree_withprivate"), "timetree.nexus"))
  write.tree(
    tr_wgs,
    file = file.path(geno, paste0(geno, "_timetree_withprivate"), paste0(geno, "_timetree.nwk"))
  )

  # Same for N450 timetree (already includes private, so no need to specify that in the path)
  tr_n450 <- read.nexus(file.path("..", "n450_tree", geno, paste0(geno, "_timetree"), "timetree.nexus"))
  write.tree(
    tr_n450,
    file = file.path("..", "n450_tree", geno, paste0(geno, "_timetree"),
                     paste0("n450_", geno, "_timetree.nwk"))
  )
}