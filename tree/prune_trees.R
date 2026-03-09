library(ape)
# For comparisons between N450 and WG trees (PS/AU tests),
# we can only use sequences that are present in both regions
# So, we trim to only consider those tips

for (geno in c("B3", "D8")) {
  # Read trees
  wgs <- read.tree(file.path(geno, paste0(geno, "_timetree_withprivate"), paste0(geno, "_timetree.nwk")))
  reg <- read.tree(file.path("..", "n450_tree", geno, paste0(geno, "_timetree"),
                             paste0("n450_", geno, "_timetree.nwk")))

  # Common taxa
  common <- intersect(wgs$tip.label, reg$tip.label)

  # Prune
  wgs_p <- drop.tip(wgs, setdiff(wgs$tip.label, common))
  reg_p <- drop.tip(reg, setdiff(reg$tip.label, common))

  # Write pruned trees
  write.tree(
    wgs_p,
    file.path(geno, paste0(geno, "_timetree_withprivate"), paste0(geno, "_timetree_pruned.nwk"))
  )
  write.tree(
    reg_p,
    file.path("..", "n450_tree", geno, paste0(geno, "_timetree"),
              paste0("n450_", geno, "_timetree_pruned.nwk"))
  )

  # Write common taxa list
  writeLines(common, paste0(geno, "_common_taxa.txt"))
}