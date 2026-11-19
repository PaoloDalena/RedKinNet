# Synthetic example cohort shipped with RedKinNet as `rkn_example`: a table of pairwise
# kinship estimates, filtered to pairs with a kinship of at least 0.02, and a table of sample
# covariates. The cohort is simulated from pedigrees; no real individual is involved.
# Run from the package root.

set.seed(20261120)

# Kinship within a pedigree whose members are ordered so that parents precede their children
# (0 for an unknown parent).
kinship_pedigree <- function(father, mother) {
  n <- length(father)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    f <- father[i]
    m <- mother[i]
    K[i, i] <- if (f > 0 && m > 0) (1 + K[f, m]) / 2 else 0.5
    if (i > 1) {
      j <- seq_len(i - 1)
      v <- (if (f > 0) K[f, j] else 0) / 2 + (if (m > 0) K[m, j] else 0) / 2
      K[i, j] <- v
      K[j, i] <- v
    }
  }
  K
}

# A founder couple with their children, grandchildren and great-grandchildren; unrelated
# spouses enter as founders. Generation 1 is the founders.
one_clan <- function(mean_children = 2.3, p_spouse = 0.85) {
  father <- c(0L, 0L)
  mother <- c(0L, 0L)
  sex <- c("M", "F")
  generation <- c(1L, 1L)
  couples <- matrix(c(1L, 2L), 1)
  for (gen in 2:4) {
    next_couples <- matrix(integer(0), 0, 2)
    for (k in seq_len(nrow(couples))) {
      for (kid in seq_len(rpois(1, mean_children))) {
        s <- sample(c("M", "F"), 1)
        father <- c(father, couples[k, 1])
        mother <- c(mother, couples[k, 2])
        sex <- c(sex, s)
        generation <- c(generation, gen)
        child <- length(father)
        if (gen < 4 && runif(1) < p_spouse) {
          father <- c(father, 0L)
          mother <- c(mother, 0L)
          sex <- c(sex, if (s == "M") "F" else "M")
          generation <- c(generation, gen)
          spouse <- length(father)
          next_couples <- rbind(next_couples, if (s == "M") c(child, spouse) else c(spouse, child))
        }
      }
    }
    couples <- next_couples
    if (nrow(couples) == 0) break
  }
  list(father = father, mother = mother, sex = sex, generation = generation)
}

# A small isolated village: n_per_gen founders and n_gen generations of couples formed at
# random within the village, each with a Poisson number of children so that the size stays
# about constant. Returns the pedigree of all generations.
village <- function(n_per_gen, n_gen) {
  father <- integer(n_per_gen)
  mother <- integer(n_per_gen)
  sex <- sample(c("M", "F"), n_per_gen, replace = TRUE)
  generation <- rep(1L, n_per_gen)
  current <- seq_len(n_per_gen)
  for (gen in 2:n_gen) {
    males <- current[sex[current] == "M"]
    females <- current[sex[current] == "F"]
    n_couples <- min(length(males), length(females))
    fa <- males[sample.int(length(males), n_couples)]
    mo <- females[sample.int(length(females), n_couples)]
    # Siblings and half siblings do not marry: their partners are reshuffled.
    for (attempt in 1:50) {
      sibs <- which((father[fa] > 0 & father[fa] == father[mo]) | (mother[fa] > 0 & mother[fa] == mother[mo]))
      if (length(sibs) == 0) break
      swap <- sample(n_couples, length(sibs))
      tmp <- mo[sibs]
      mo[sibs] <- mo[swap]
      mo[swap] <- tmp
    }
    kids <- rpois(n_couples, n_per_gen / n_couples)
    father <- c(father, rep(fa, kids))
    mother <- c(mother, rep(mo, kids))
    sex <- c(sex, sample(c("M", "F"), sum(kids), replace = TRUE))
    generation <- c(generation, rep(gen, sum(kids)))
    current <- which(generation == gen)
  }
  list(father = father, mother = mother, sex = sex, generation = generation)
}

# Estimated kinship: the pedigree kinship plus the variation of the realised relatedness, with
# variance proportional to the kinship (standard deviation 0.02 for full siblings, none for parent
# and child), and an estimation error.
estimate <- function(phi, parent_child) {
  sd <- sqrt(0.0016 * phi * (!parent_child) + 0.003^2)
  round(phi + rnorm(length(phi), 0, sd), 4)
}

# Ages by generation below the youngest (0 = youngest sampled generation).
ages <- function(gen_from_youngest) {
  lo <- c(18, 44, 68)[gen_from_youngest + 1]
  as.integer(lo + floor(runif(length(gen_from_youngest), 0, 24)))
}

samples <- list()
pairs <- list()
add_group <- function(ped, picked, youngest, label) {
  ids <- sprintf("%s_%03d", label, picked)
  samples[[length(samples) + 1]] <<- data.frame(key = ids, sex = ped$sex[picked],
                                                age = ages(youngest - ped$generation[picked]),
                                                stringsAsFactors = FALSE)
  if (length(picked) > 1) {
    K <- kinship_pedigree(ped$father, ped$mother)[picked, picked, drop = FALSE]
    ij <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
    if (nrow(ij) > 0) {
      a <- picked[ij[, 1]]
      b <- picked[ij[, 2]]
      parent_child <- ped$father[a] == b | ped$mother[a] == b | ped$father[b] == a | ped$mother[b] == a
      pairs[[length(pairs) + 1]] <<- data.frame(key1 = ids[ij[, 1]], key2 = ids[ij[, 2]],
                                                kinship = estimate(K[ij], parent_child), stringsAsFactors = FALSE)
    }
  }
}

# Families: clans sampled sparsely, so that most samples have no relatives in the cohort.
n_clan_samples <- 0
clan <- 0
while (n_clan_samples < 720) {
  clan <- clan + 1
  ped <- one_clan()
  eligible <- which(ped$generation >= 2)
  picked <- eligible[runif(length(eligible)) < 0.12]
  if (length(picked) == 0) next
  add_group(ped, picked, 4L, sprintf("c%04d", clan))
  n_clan_samples <- n_clan_samples + length(picked)
}

# The village: most members of its last two generations are sampled.
vil <- village(n_per_gen = 40, n_gen = 5)
eligible <- which(vil$generation >= 4)
add_group(vil, eligible[runif(length(eligible)) < 0.4], 5L, "v")

samples <- do.call(rbind, samples)
pairs <- do.call(rbind, pairs)

# One sample genotyped twice (a member of a family) and one contaminated sample, whose kinship
# with every other sample is raised.
related <- unique(c(pairs$key1, pairs$key2))
twice <- sample(setdiff(related, grep("^v_", related, value = TRUE)), 1)
contaminated <- sample(setdiff(samples$key, related), 1)
copy <- paste0(twice, "_b")
hit <- pairs$key1 == twice | pairs$key2 == twice
partner <- ifelse(pairs$key1[hit] == twice, pairs$key2[hit], pairs$key1[hit])
pairs <- rbind(pairs,
               data.frame(key1 = copy, key2 = partner, kinship = round(pairs$kinship[hit] + rnorm(sum(hit), 0, 0.003), 4)),
               data.frame(key1 = twice, key2 = copy, kinship = round(0.5 - abs(rnorm(1, 0, 0.005)), 4)))
samples <- rbind(samples, transform(samples[samples$key == twice, ], key = copy))
others <- setdiff(samples$key, contaminated)
pairs <- rbind(pairs, data.frame(key1 = contaminated, key2 = others,
                                 kinship = round(0.05 + rnorm(length(others), 0, 0.006), 4)))

# Identifiers in random order, covariates with some missing values, and the filter of the
# kinship table at 0.02.
samples <- samples[sample.int(nrow(samples)), ]
samples$id <- sprintf("S%03d", seq_len(nrow(samples)))
samples$call_rate <- round(1 - rbeta(nrow(samples), 1, 150), 4)
samples$sex[runif(nrow(samples)) < 0.02] <- NA
samples$age[runif(nrow(samples)) < 0.03] <- NA
pairs$id1 <- samples$id[match(pairs$key1, samples$key)]
pairs$id2 <- samples$id[match(pairs$key2, samples$key)]
pairs <- pairs[pairs$kinship >= 0.02, ]
swap <- pairs$id1 > pairs$id2
pairs[swap, c("id1", "id2")] <- pairs[swap, c("id2", "id1")]
pairs <- pairs[order(pairs$id1, pairs$id2), ]

rkn_example <- list(
  kinship = data.frame(id1 = pairs$id1, id2 = pairs$id2, kinship = pairs$kinship, stringsAsFactors = FALSE),
  samples = data.frame(id = samples$id, sex = samples$sex, age = samples$age, call_rate = samples$call_rate,
                       stringsAsFactors = FALSE)
)
rownames(rkn_example$samples) <- NULL
rownames(rkn_example$kinship) <- NULL

usethis::use_data(rkn_example, overwrite = TRUE, compress = "xz")
