library(ape)
library(phangorn)

## script for calculating parsimony scores, comparing WG and N450 sequencing approaches. 

paired_perm_test_ps <- function(tree_wgs, tree_n450, trait_wgs, nperm = 5000, seed = 1) {
  set.seed(seed)
  lev <- sort(unique(trait_wgs))  # keep levels fixed across both trees

  ps_from_trait <- function(tree, trait_named) {
    parsimony(tree, phyDat(trait_named, type = "USER", levels = lev))
  }

  # observed
  ps_w_obs <- ps_from_trait(tree_wgs, trait_wgs)
  trait_n450 <- trait_wgs[tree_n450$tip.label]
  names(trait_n450) <- tree_n450$tip.label
  ps_n_obs <- ps_from_trait(tree_n450, trait_n450)
  delta_obs <- ps_n_obs - ps_w_obs

  # null
  delta_null <- replicate(nperm, {
    perm <- sample(trait_wgs)
    names(perm) <- names(trait_wgs)

    ps_w <- ps_from_trait(tree_wgs, perm)

    perm_n <- perm[tree_n450$tip.label]
    names(perm_n) <- tree_n450$tip.label
    ps_n <- ps_from_trait(tree_n450, perm_n)

    ps_n - ps_w
  })

  # empirical p value testing, using resampling
  p <- (sum(delta_null >= delta_obs) + 1) / (nperm + 1)

  list(delta_obs = delta_obs, delta_null = delta_null, p = p,
       ps_wgs = ps_w_obs, ps_n450 = ps_n_obs)
}


results <- list()

for (geno in c("D8", "B3")) {
  message("Running ", geno)

  meta <- read.delim(file.path(geno, paste0(geno, "_all_meta.tsv")), stringsAsFactors = FALSE)
  meta$country <- sapply(strsplit(meta$strain, "\\|"), `[`, 2)

  tree_wgs  <- read.tree(file.path(geno, paste0(geno, "_timetree_withprivate"), paste0(geno, "_timetree_pruned.nwk")))
  tree_n450 <- read.tree(file.path("..", "n450_tree", geno, paste0(geno, "_timetree"),
                                  paste0("n450_", geno, "_timetree_pruned.nwk")))

  meta$strain  <- as.character(meta$strain)
  meta$country <- as.character(meta$country)

  trait_all <- setNames(meta$country, meta$strain)

  trait_wgs <- trait_all[tree_wgs$tip.label]
  names(trait_wgs) <- tree_wgs$tip.label

  trait_n450 <- trait_all[tree_n450$tip.label]
  names(trait_n450) <- tree_n450$tip.label

  message("running paired PS test")
  ps_paired <- paired_perm_test_ps(tree_wgs, tree_n450, trait_wgs, nperm = 10000, seed = 1)
  message(sprintf("paired PS delta (n450 - wgs): %s ; p=%s", ps_paired$delta_obs, ps_paired$p))

}

