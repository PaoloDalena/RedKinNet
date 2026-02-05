# rules_db.R
# Definition of structural rules for connected components of size <= 5.
# Each rule block encodes:
#   - a matcher (match_fn)
#   - nodes always kept (kept_rules)
#   - nodes always discarded (discarded_rules)
#   - mutually exclusive candidate groups (conflict_rules)

rules_db <- list(

  ## size 1 ---------------------------------------------------------------
  list(
    name = "single_node",
    size_allowed = 1,
    match_fn = function(fam_df) {
      nrow(fam_df) == 1
    },
    kept_rules = list(
      c("TRUE")
    ),
    discarded_rules = list(),
    conflict_rules = list()
  ),

  ## size 2 ---------------------------------------------------------------
  ## single edge: maximum independent set size = 1
  list(
    name = "pair_edge",
    size_allowed = 2,
    match_fn = function(fam_df) {
      nrow(fam_df) == 2 &&
        all(fam_df$degree >= 1)
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    )
  ),

  ## size 3: path (degrees 1,1,2) ----------------------------------------
  list(
    name = "path3",
    size_allowed = 3,
    match_fn = function(fam_df) {
      nrow(fam_df) == 3 &&
        all(sort(fam_df$degree) == c(1,1,2))
    },
    kept_rules = list(
      c("degree == 1")
    ),
    discarded_rules = list(
      c("degree == 2")
    ),
    conflict_rules = list()
  ),

  ## size 3: triangle / K3 (degrees 2,2,2) --------------------------------
  list(
    name = "triangle3",
    size_allowed = 3,
    match_fn = function(fam_df) {
      nrow(fam_df) == 3 &&
        length(unique(fam_df$degree)) == 1 &&
        unique(fam_df$degree) == 2
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    )
  ),

  ## size 4: star (degrees 1,1,1,3) --------------------------------------
  list(
    name = "star4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        all(sort(fam_df$degree) == c(1,1,1,3))
    },
    kept_rules = list(
      c("degree == 1")
    ),
    discarded_rules = list(
      c("degree == 3")
    ),
    conflict_rules = list()
  ),

  ## size 4: path (degrees 1,1,2,2) --------------------------------------
  ## multiple maximum independent sets, no fixed keep/discard
  list(
    name = "path4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        all(sort(fam_df$degree) == c(1,1,2,2)) &&
        sum(fam_df$degree == 1) == 2 &&
        sum(fam_df$degree == 2) == 2
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    ),
    alpha = 2
  ),

  ## size 4: 4-cycle / square (degrees 2,2,2,2) --------------------------
  ## checkerboard-like solutions, no fixed keep/discard
  list(
    name = "cycle4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        length(unique(fam_df$degree)) == 1 &&
        unique(fam_df$degree) == 2
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    ),
    alpha = 2
  ),

  ## size 4: clique K4 (degrees 3,3,3,3) ----------------------------------
  ## maximum independent set size = 1
  list(
    name = "clique4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        length(unique(fam_df$degree)) == 1 &&
        unique(fam_df$degree) == 3
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    )
  ),

  ## size 4: triangle + tail (degrees 1,2,2,3) ---------------------------
  list(
    name = "triangle_tail4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        all(sort(fam_df$degree) == c(1,2,2,3)) &&
        sum(fam_df$degree == 1) == 1 &&
        sum(fam_df$degree == 2) == 2 &&
        sum(fam_df$degree == 3) == 1
    },
    kept_rules = list(
      c("degree == 1")
    ),
    discarded_rules = list(
      c("degree == 3")
    ),
    conflict_rules = list(
      list("degree == 2")
    )
  ),

  ## size 4: "diamond"/"kite" (degrees 2,2,3,3) --------------------------
  ## structurally dense core, no universal keep
  list(
    name = "diamond4",
    size_allowed = 4,
    match_fn = function(fam_df) {
      nrow(fam_df) == 4 &&
        all(sort(fam_df$degree) == c(2,2,3,3)) &&
        sum(fam_df$degree == 2) == 2 &&
        sum(fam_df$degree == 3) == 2
    },
    kept_rules = list(list("degree == 2") ),
    discarded_rules = list(list("degree == 3")),
    conflict_rules = list()
  ),

  ## size 5: star (degrees 1,1,1,1,4) ------------------------------------
  list(
    name = "star5",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        all(sort(fam_df$degree) == c(1,1,1,1,4)) &&
        sum(fam_df$degree == 1) == 4 &&
        sum(fam_df$degree == 4) == 1
    },
    kept_rules = list(
      c("degree == 1")
    ),
    discarded_rules = list(
      c("degree == 4")
    ),
    conflict_rules = list()
  ),

  ## size 5: path (degrees 1,1,2,2,2) ------------------------------------
  ## odd path: unique maximum independent set (deterministic)
  list(
    name = "path5",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        all(sort(fam_df$degree) == c(1,1,2,2,2)) &&
        sum(fam_df$degree == 1) == 2 &&
        sum(fam_df$degree == 2) == 3
    },
    kept_rules = list(
      c("degree == 1"),
      c("degree == 2 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("max_path_length == 3")
    ),
    conflict_rules = list()
  ),

  ## size 5: 5-cycle / pentagon (degrees 2,2,2,2,2) ----------------------
  ## rotational symmetry, multiple equivalent max sets
  list(
    name = "cycle5",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        all(fam_df$degree == 2L) &&                # all degree 2
        all(fam_df$max_path_length == 2)           # distinguish from path5
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list( list("TRUE") ),         # all candidates
    alpha = 2                                      # uses largest_ivs()
  ),

  ## size 5: clique K5 (degrees 4,4,4,4,4) --------------------------------
  ## maximum independent set size = 1
  list(
    name = "clique5",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        length(unique(fam_df$degree)) == 1 &&
        unique(fam_df$degree) == 4
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list("TRUE")
    )
  ),

  ## size 5: triangle with two tails (degrees 1,1,2,2,3) -----------------
  list(
    name = "triangle_two_tails5",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        all(sort(fam_df$degree) == c(1,1,2,2,3)) &&
        sum(fam_df$degree == 1) == 2 &&
        sum(fam_df$degree == 2) == 2 &&
        sum(fam_df$degree == 3) == 1
    },
    kept_rules = list(
      c("degree == 1")
    ),
    discarded_rules = list(
      c("degree == 3")
    ),
    conflict_rules = list(
      list("degree == 2")
    )
  ),

  list(
    name = "hub_tail5_iso17_determ",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,2,2,2,3" &&
        paste(
          sort(paste(fam_df$degree, fam_df$max_path_length, sep=":")),
          collapse = ","
        ) == "1:3,2:2,2:2,2:3,3:2"
    },

    kept_rules = list(
      c("(degree == 2 & max_path_length == 2) | degree == 1")
    ),

    discarded_rules = list(
      c("(degree == 3) | (degree == 2 & max_path_length == 3)")
    ),

    conflict_rules = list()
  ),



  ## size 5: asymmetric hub + tails (degrees 1,2,2,2,3) ------------------
  list(
    name = "hub_tail5_like_iso15",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,2,2,2,3"
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list( list("TRUE") ),   # all candidates
    alpha = 2                                # uses largest_ivs
  ),

  ## size 5: degree signature 1,1,1,2,3 ----------------------------------
  ## structurally ambiguous leaves; fallback to global resolution
  list(
    name = "size5_sig_1,1,1,2,3",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,1,1,2,3"
    },
    # a priori tabellari
    kept_rules = list(),
    discarded_rules = list(
      c("degree == 3 & max_path_length == 2")  # discard the center (node 1)
    ),
    # a priori graph-dependent: leaves attached to grade 3 node
    kept_fn = function(df, g) {
      deg <- igraph::degree(g)
      leaves <- df$.name[df$degree == 1]
      keep <- leaves[vapply(leaves, function(v) {
        any(deg[igraph::neighbors(g, v)] == 3)
      }, logical(1))]
      keep
    },
    conflict_rules = list( list("TRUE") )
  ),

  ## size 5: degree signature 1,1,2,2,4 ----------------------------------
  list(
    name = "size5_sig_1,1,2,2,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,1,2,2,4"
    },
    kept_rules = list(
      c("degree == 1 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list(
      list("degree == 2 & max_path_length == 2")
    )
  ),

  ## size 5: degree signature 1,1,2,3,3 ----------------------------------
  list(
    name = "size5_sig_1,1,2,3,3",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,1,2,3,3"
    },
    kept_rules = list(
      c("degree == 2 & max_path_length == 2"),
      c("degree == 1 & max_path_length == 3")
    ),
    discarded_rules = list(
      c("degree == 3 & max_path_length == 2")
    ),
    conflict_rules = list()
  ),

  ## size 5: degree signature 1,2,2,3,4 ----------------------------------
  list(
    name = "size5_sig_1,2,2,3,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,2,2,3,4"
    },
    kept_rules = list(
      c("degree == 2 & max_path_length == 2"),
      c("degree == 1 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1"),
      c("degree == 3 & max_path_length == 2")
    ),
    conflict_rules = list()
  ),

  ## size 5: degree signature 1,2,3,3,3 ----------------------------------
  ## no fixed keep/discard; single global conflict group
  list(
    name = "size5_sig_1,2,3,3,3",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,2,3,3,3"
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list( list("TRUE") ),   # all candidates
    alpha = 2                                # largest_ivs
  ),

  ## size 5: degree signature 1,3,3,3,4 ----------------------------------
  list(
    name = "size5_sig_1,3,3,3,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "1,3,3,3,4"
    },
    kept_rules = list(
      c("degree == 1 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list(
      list("degree == 3 & max_path_length == 2")
    )
  ),

  ## size 5: degree signature 2,2,2,3,3 ----------------------------------
  list(
    name = "size5_sig_2,2,2,3,3",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,2,2,3,3"
    },
    kept_rules = list(
      c("degree == 2")
    ),
    discarded_rules = list(
      c("degree == 3")
    ),
    conflict_rules = list()
  ),


  ## size 5: degree signature 2,2,2,4,4 ----------------------------------
  list(
    name = "size5_sig_2,2,2,4,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,2,2,4,4"
    },
    kept_rules = list(
      c("degree == 2 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list()
  ),

  ## size 5: degree signature 2,2,2,2,4 ----------------------------------
  ## two disjoint conflict groups, same local pattern
  list(
    name = "size5_sig_2,2,2,2,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,2,2,2,4"
    },
    kept_rules = list(),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list(
      list("degree == 2 & max_path_length == 2"),
      list("degree == 2 & max_path_length == 2")
    )
  ),

  ## size 5: degree signature 2,2,2,3,3 (alternate form) -----------------
  ## no fixed keep/discard; global conflict
  list(
    name = "size5_sig_2,2,2,3,3_alt",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,2,2,3,3"
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list(
      list(
        "degree == 3 & max_path_length == 2",
        "degree == 2 & max_path_length == 2"
      )
    )
  ),

  ## size 5: degree signature 2,2,3,3,4 ----------------------------------
  list(
    name = "size5_sig_2,2,3,3,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,2,3,3,4"
    },
    kept_rules = list(),
    discarded_rules = list( list("degree == 4") ),
    conflict_rules = list( list("TRUE") ),
    alpha = 2
  ),

  ## size 5: degree signature 2,3,3,3,3 ----------------------------------
  ## no fixed keep/discard; global conflict
  list(
    name = "size5_sig_2,3,3,3,3",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,3,3,3,3"
    },
    kept_rules = list(),
    discarded_rules = list(),
    conflict_rules = list( list("TRUE") ),   # all candidates
    alpha = 2                                # largest_ivs()
  ),

  ## size 5: degree signature 2,3,3,4,4 ----------------------------------
  list(
    name = "size5_sig_2,3,3,4,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "2,3,3,4,4"
    },
    kept_rules = list(
      c("degree == 2 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list(
      list("degree == 3 & max_path_length == 2")
    )
  ),

  ## size 5: degree signature 3,3,3,3,4 ----------------------------------
  list(
    name = "size5_sig_3,3,3,3,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "3,3,3,3,4"
    },
    kept_rules = list(),
    discarded_rules = list(
      list("degree == 4")                      # discard the hub (1)
    ),
    conflict_rules = list(
      list("degree == 3")                      # candidates = {2,3,4,5}
    ),
    alpha = 2                                  # largest_ivs()
  ),

  ## size 5: degree signature 3,3,4,4,4 ----------------------------------
  list(
    name = "size5_sig_3,3,4,4,4",
    size_allowed = 5,
    match_fn = function(fam_df) {
      nrow(fam_df) == 5 &&
        paste(sort(fam_df$degree), collapse = ",") == "3,3,4,4,4"
    },
    kept_rules = list(
      c("degree == 3 & max_path_length == 2")
    ),
    discarded_rules = list(
      c("degree == 4 & max_path_length == 1")
    ),
    conflict_rules = list()
  )

)
