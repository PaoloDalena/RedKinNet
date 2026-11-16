#' Selections of unrelated individuals
#'
#' @description
#' \code{\link{rkn_select_random}}, \code{\link{rkn_select_criteria}} and
#' \code{\link{rkn_select_global}} return an object of class \code{rkn_selection}: a
#' list holding the selected individuals, the network they were selected from and,
#' when computed, diagnostics on how each individual and each family was resolved. The object has \code{\link[=print.rkn_selection]{print}},
#' \code{\link[=summary.rkn_selection]{summary}} and
#' \code{\link[=plot.rkn_selection]{plot}} methods.
#'
#' @section Components:
#' \describe{
#'   \item{\code{kept}}{character vector, the identifiers of the selected individuals;}
#'   \item{\code{discarded}}{character vector, the identifiers of the other individuals;}
#'   \item{\code{family_sizes}}{integer vector, the size of each family (connected
#'     component of the network);}
#'   \item{\code{method}}{\code{"random"}, \code{"criteria"} or \code{"global"};}
#'   \item{\code{uniform}}{logical, whether the selection was drawn uniformly among the
#'     maximum sets (\code{rkn_select_random(uniform = TRUE)});}
#'   \item{\code{optimal}}{logical, whether the criteria were evaluated over all the
#'     maximum sets of each family (\code{rkn_select_criteria(optimal = TRUE)});}
#'   \item{\code{criteria}}{for criteria-based selections, a list with the
#'     \code{primary} and, if given, the \code{secondary} criterion; \code{NULL}
#'     otherwise;}
#'   \item{\code{seed}}{the seed of the selection, or \code{NULL};}
#'   \item{\code{diagnostics}}{logical, whether diagnostics were computed;}
#'   \item{\code{node_diagnostics}}{a data frame with one row per individual, or
#'     \code{NULL} without diagnostics. Columns: \code{id}; \code{pre_class}, how the
#'     individual was decided (\code{"kept_apriori"}, \code{"discarded_apriori"},
#'     \code{"candidate_simple"} or \code{"candidate_complex"}, described in the section
#'     Types of decision below); \code{post_class}, \code{"kept"} or
#'     \code{"discarded"}; \code{is_selected}, the outcome as a logical;
#'     \code{pre_reason}, an internal label of the rule or procedure that decided the
#'     individual; \code{family_id} and \code{family_size}; in uniform selections,
#'     \code{inclusion_probability}, the exact probability that the individual is
#'     selected, i.e. the proportion of the maximum sets of its family that contain it
#'     (\code{NA} in the families too intricate for uniform sampling). From these
#'     probabilities \code{print()} reports the expected number of selected individuals
#'     that differ between two independent uniform draws, the sum of \eqn{p(1 - p)} over
#'     the individuals: the part of the selection left to chance;}
#'   \item{\code{family_summary}}{a data frame with one row per family, or \code{NULL}
#'     without diagnostics. Columns: \code{family_id} and \code{family_size};
#'     \code{selector}, \code{"upto5"} for families resolved by the structural rules in
#'     \code{\link{rules_db}}, \code{"6plus"} for larger families, \code{"uniform"}
#'     for families drawn uniformly and \code{"optimal"} for families where the criteria
#'     were evaluated over all maximum sets;
#'     \code{structure_name}, an internal label of the recognised structure;
#'     \code{mis_method}, how the individuals left undecided by the structural rules
#'     were resolved (see the section Search over maximum sets): \code{NA} when the
#'     rules decided everyone, \code{"enumeration"}, \code{"exact"},
#'     \code{"exact_reduced"} or \code{"exact_sample"}; \code{uniform}, whether the
#'     family was drawn uniformly among its maximum sets; \code{optimal}, whether the
#'     criteria were evaluated over all its maximum sets; \code{n_max_sets}, in the
#'     families drawn uniformly or evaluated over all maximum sets, the number of maximum
#'     sets of the family (\code{NA} otherwise, and after \code{"exact_reduced"}).
#'     Since the choices in different families are independent, the maximum sets of the
#'     network are all the combinations of one maximum set per family; their number,
#'     reported by \code{print()}, is the product of \code{n_max_sets} over the families
#'     and grows exponentially with the number of families that have a choice (each
#'     related pair doubles it), so \code{print()} also reports how many families have
#'     more than one maximum set; the numbers of kept and discarded individuals, overall
#'     (\code{n_kept}, \code{n_discarded}) and by type of decision
#'     (\code{n_kept_apriori}, \code{n_discarded_apriori}, \code{n_candidate_simple},
#'     \code{n_candidate_complex}). In global selections the columns are
#'     \code{family_id}, \code{family_size}, \code{selector} (\code{"global"}),
#'     \code{exact} (\code{FALSE} for the families beyond the computational budget,
#'     see \code{\link{rkn_select_global}}), \code{n_kept}, \code{n_discarded},
#'     \code{n_category} (selected members in the category) and \code{n_missing}
#'     (selected members with a missing value), and the types of decision of
#'     \code{node_diagnostics} are not available;}
#'   \item{\code{graph}}{the igraph network the selection was made on;}
#'   \item{\code{call}}{the call that produced the selection;}
#'   \item{\code{global}}{for global selections, a list describing the composition
#'     reached: \code{covariate}, \code{category} and \code{target}; \code{cohort},
#'     whether the target is the proportion of the network; \code{proportion}, the
#'     proportion reached among the selected individuals with a known value,
#'     \code{n_category} of \code{n_known}; \code{n_missing}; \code{size},
#'     \code{max_size} and \code{max_loss}; \code{log_count}, the natural logarithm of
#'     the number of selections reaching the same proportion and size (\code{NA} when
#'     \code{exact} is \code{FALSE}); \code{exact}, whether every family was computed
#'     exactly; \code{n_inexact}, the number of families that were not. \code{NULL}
#'     otherwise.}
#' }
#'
#' @section Family structures:
#' A family is a connected component of the relatedness network: a group of
#' individuals linked, directly or through other members, by relatedness above the
#' threshold used to build the network. \code{\link{summary.rkn_selection}} groups the
#' families into seven structures:
#' \describe{
#'   \item{single individual}{an individual related to no one else in the cohort;}
#'   \item{related pair}{two related individuals;}
#'   \item{chain of three}{three individuals of whom two are not related to each other
#'     but both are related to the third, typically two parents and their child;}
#'   \item{three mutually related}{three individuals all related to each other, for
#'     instance three siblings or a parent and two children;}
#'   \item{family of 4, of 5, and of 6 or more members}{larger families. Families of up
#'     to 5 members are resolved by the structural rules in \code{\link{rules_db}},
#'     larger families by reducing the family step by step with structural rules and
#'     resolving what remains.}
#' }
#'
#' @section Search over maximum sets:
#' The structural rules decide many individuals directly. The individuals they leave
#' undecided are resolved by a search over the maximum sets of unrelated individuals,
#' reported in the column \code{mis_method} of \code{family_summary}:
#' \describe{
#'   \item{\code{"enumeration"}}{in families of up to 5 members, all the maximum sets
#'     are listed and one of them is chosen;}
#'   \item{\code{"exact"}}{in families of 6 or more members, an exact algorithm
#'     (dynamic programming along an elimination ordering of the family) finds the
#'     maximum size and draws one of the maximum sets uniformly at random. Its cost
#'     grows with how intricate the relationships are, not with the size of the family,
#'     so it also handles the large families of genetic isolates;}
#'   \item{\code{"exact_reduced"}}{the same after removing the individuals who cannot
#'     enlarge the selection (an individual is removed when one of their relatives has
#'     no relative outside their own). This is used only for families too intricate for
#'     the exact algorithm: the set chosen is still of maximum size, but not drawn
#'     uniformly. In criteria-based selections an individual is removed only if that
#'     relative meets the criteria \code{level} or \code{cutoff} at least as well, so
#'     these criteria are still optimised exactly; \code{target_mean} and
#'     \code{category_balance} are evaluated on the maximum sets that remain.}
#'   \item{\code{"exact_sample"}}{criteria-based selections only: see below.}
#' }
#' In criteria-based selections (\code{\link{rkn_select_criteria}}), the criteria
#' \code{level} and \code{cutoff} are optimised exactly by the same algorithm when no
#' undecided member of the family lacks the covariate, because among sets of the same size
#' they rank sets by the number of members meeting the condition. The criteria
#' \code{target_mean} and \code{category_balance}, and \code{level} or \code{cutoff}
#' with missing covariate values, are evaluated on all the maximum sets when there are
#' at most 10,000; above that number, as in some large families of genetic isolates,
#' they are evaluated on a uniform random sample of 10,000 maximum sets, and the best of
#' the sample is kept (\code{"exact_sample"}). The number of such families is reported
#' by \code{print()}.
#' If a family remains too intricate even after this reduction, the selection stops
#' with an error that suggests excluding the samples flagged by
#' \code{\link{rkn_check_samples}} or using a stricter threshold.
#'
#' @section Types of decision:
#' Within each family, every individual is decided in one of four ways, reported by
#' \code{\link{summary.rkn_selection}} and shown by the \code{"diagnostics"} plot (see
#' \code{\link{plot.rkn_selection}}):
#' \describe{
#'   \item{kept without a choice}{the individual is kept with no alternative to choose
#'     from: a rule based on the structure of the family keeps them (for instance the
#'     leaves of a star-shaped family), or they belong to every maximum set of the
#'     individuals left undecided by the rules;}
#'   \item{discarded without a choice}{likewise, the individual is discarded by a
#'     structural rule (for instance the centre of a star) or belongs to no maximum set
#'     of the individuals left undecided;}
#'   \item{decided among single candidates}{the remaining choice is among
#'     interchangeable candidates of a single group, such as the two members of a
#'     related pair, and one of them is selected;}
#'   \item{decided among alternative sets}{several maximum sets remain that differ in
#'     more than one individual, and one of them is selected.}
#' }
#' The choices among candidates are made at random by
#' \code{\link{rkn_select_random}} and by the criteria in
#' \code{\link{rkn_select_criteria}}.
#'
#' "Without a choice" describes the procedure, not the family: to simplify the search,
#' the structural rules fix some individuals that other maximum sets of the family leave
#' out (for instance, in a chain of relatives they keep the individual at the end), so an
#' individual kept without a choice is not necessarily in every maximum set of the
#' family. The exact classes over the whole family are given by uniform selections.
#'
#' In uniform selections (\code{rkn_select_random(uniform = TRUE)}) no structural rule
#' is applied, and the same codes have an exact meaning, computed from the inclusion
#' probabilities: \code{"kept_apriori"}, in every maximum set (probability 1);
#' \code{"discarded_apriori"}, in no maximum set (probability 0);
#' \code{"candidate_complex"}, chosen at random among the maximum sets (probability
#' strictly between 0 and 1). \code{summary()} and \code{plot()} label them
#' accordingly.
#'
#' In criteria-based selections over all maximum sets
#' (\code{rkn_select_criteria(optimal = TRUE)}) no structural rule is applied either, so
#' "kept without a choice" and "discarded without a choice" mean in every and in no
#' maximum set of the family, and "decided among alternative sets" means chosen by the
#' criteria among maximum sets that differ.
#'
#' @example man-roxygen/example-network.R
#' @examples
#' sel <- rkn_select_random(g, seed = 1)
#' names(sel)
#' sel$kept
#' head(sel$node_diagnostics)
#' sel$family_summary[, c("family_id", "family_size", "structure_name")]
#'
#' @name rkn_selection
NULL

.rkn_structure_classes <- c("single individual", "related pair",
                            "chain of three (e.g. parents and child)", "three mutually related",
                            "family of 4 members", "family of 5 members",
                            "family of 6 or more members")

.rkn_decision_labels <- c(kept_apriori = "kept without a choice",
                          discarded_apriori = "discarded without a choice",
                          candidate_simple = "decided among single candidates",
                          candidate_complex = "decided among alternative sets")

# In uniform selections the same codes have an exact meaning (see ?rkn_selection).
.rkn_decision_labels_uniform <- c(kept_apriori = "in every maximum set",
                                  discarded_apriori = "in no maximum set",
                                  candidate_simple = "decided among single candidates",
                                  candidate_complex = "chosen at random among maximum sets")

#' Internal helper: decision labels of a selection
#'
#' @keywords internal
#' @noRd
.rkn_labels_of <- function(x) {
  if (isTRUE(x$uniform)) .rkn_decision_labels_uniform else .rkn_decision_labels
}

#' Internal constructor for rkn_selection objects
#'
#' @keywords internal
#' @noRd
.new_rkn_selection <- function(kept, discarded, family_sizes, method, criteria, seed,
                               node_diagnostics, family_summary, graph, call, uniform = FALSE,
                               optimal = FALSE, global = NULL) {
  structure(
    list(
      kept = kept,
      discarded = discarded,
      family_sizes = family_sizes,
      method = method,
      uniform = uniform,
      optimal = optimal,
      criteria = criteria,
      seed = seed,
      diagnostics = !is.null(node_diagnostics),
      node_diagnostics = node_diagnostics,
      family_summary = family_summary,
      graph = graph,
      call = call,
      global = global
    ),
    class = "rkn_selection"
  )
}

#' Internal helper: one-line description of how a selection was made
#'
#' @keywords internal
#' @noRd
.rkn_selection_header <- function(x) {
  header <- if (isTRUE(x$uniform)) "Uniform random selection" else "Random selection"
  if (identical(x$method, "criteria")) {
    crit <- x$criteria$primary$type
    if (!is.null(x$criteria$secondary)) crit <- paste0(crit, ", then ", x$criteria$secondary$type)
    header <- sprintf(if (isTRUE(x$optimal)) "Criteria-based selection over all maximum sets (%s)" else
      "Criteria-based selection (%s)", crit)
  }
  if (identical(x$method, "global")) {
    header <- sprintf("Global selection (proportion of %s = %s)", x$global$covariate, x$global$category)
  }
  if (!is.null(x$seed)) header <- paste0(header, ", seed ", x$seed)
  header
}

#' Print an rkn_selection object
#'
#' Prints a brief summary of a selection of unrelated individuals: how the selection
#' was made, how many individuals were retained, and how the network decomposes into
#' families.
#'
#' @param x An \code{rkn_selection} object, as returned by \code{\link{rkn_select_random}}
#'   or \code{\link{rkn_select_criteria}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rkn_selection}}
#'
#' @example man-roxygen/example-network.R
#' @examples
#' sel <- rkn_select_random(g, seed = 1)
#' sel # the same as print(sel)
#'
#' print(rkn_select_random(g, seed = 1, diagnostics = FALSE))
#'
#' @export
print.rkn_selection <- function(x, ...) {
  n_ind <- length(x$kept) + length(x$discarded)
  sizes <- x$family_sizes

  cat(sprintf("<rkn_selection> %s\n", .rkn_selection_header(x)))
  cat(sprintf("Call: %s\n", paste(trimws(deparse(x$call)), collapse = " ")))
  cat(sprintf("Individuals: %d | selected: %d (%.1f%%)\n",
              n_ind, length(x$kept), 100 * length(x$kept) / n_ind))
  cat(sprintf("Families: %d | size min/mean/max: %d / %.2f / %d\n",
              length(sizes), min(sizes), mean(sizes), max(sizes)))
  cat(sprintf("Families with up to 5 members: %d | with 6 or more: %d\n",
              sum(sizes <= 5), sum(sizes >= 6)))

  if (identical(x$method, "global")) {
    gl <- x$global
    cat(sprintf("Proportion of %s = %s among the selected: %.4f (%d of %d%s; target %.4g%s)\n",
                gl$covariate, gl$category, gl$proportion, gl$n_category, gl$n_known,
                if (gl$n_missing > 0) " with a known value" else "", gl$target,
                if (gl$cohort) ", the proportion of the network" else ""))
    if (gl$n_missing > 0) cat(sprintf("Selected with a missing value of %s: %d\n", gl$covariate, gl$n_missing))
    if (gl$max_loss > 0) {
      cat(sprintf("Selected: %d of a largest possible %d (max_loss = %d)\n", gl$size, gl$max_size, gl$max_loss))
    }
    if (gl$exact) {
      cat(sprintf("Selections reaching this proportion%s: %s, one drawn uniformly\n",
                  if (gl$max_loss > 0) " and size" else "", .rkn_format_log10(gl$log_count / log(10))))
    } else {
      cat(sprintf(paste0("Families beyond the computational budget: %d; the proportion may not be the best ",
                         "possible, and the choice within them is not uniform\n"), gl$n_inexact))
    }
  }

  if (x$diagnostics && !identical(x$method, "global")) {
    fs <- x$family_summary
    if (isTRUE(x$uniform) || isTRUE(x$optimal)) {
      whole <- if (isTRUE(x$uniform)) fs$uniform else fs$optimal
      cat(sprintf(if (isTRUE(x$uniform)) "Families sampled uniformly among their maximum sets: %d of %d\n" else
        "Families where the criteria were evaluated over all maximum sets: %d of %d\n", sum(whole), length(sizes)))
      if (!all(whole)) {
        cat(sprintf("Families too intricate for %s, selected by the structural rules: %d\n",
                    if (isTRUE(x$uniform)) "uniform sampling" else "this evaluation", sum(!whole)))
      }
      counted <- !is.na(fs$n_max_sets)
      if (any(counted)) {
        n_sets <- fs$n_max_sets[counted]
        cat(sprintf("Maximum sets: %s in total%s; families with more than one: %d of %d\n",
                    .rkn_format_log10(sum(log10(n_sets))),
                    if (all(counted)) "" else " over the families where they were counted",
                    sum(n_sets > 1), sum(counted)))
      }
      if (isTRUE(x$uniform) && any(whole)) {
        # An individual is in exactly one of two independent draws with probability
        # 2p(1 - p); the draws have the same size, so half of these differ on each side.
        nd <- x$node_diagnostics
        p <- nd$inclusion_probability
        sampled <- !is.na(p)
        n_selected <- sum(nd$post_class[sampled] == "kept")
        expected <- sum(p[sampled] * (1 - p[sampled]))
        cat(sprintf("Two uniform draws differ on average in %.1f of the %d selected individuals%s (%.1f%%)\n",
                    expected, n_selected, if (all(sampled)) "" else " of the families sampled uniformly",
                    100 * expected / n_selected))
      }
    } else {
      n_mis <- sum(!is.na(fs$mis_method))
      cat(sprintf("Families resolved by a search over maximum independent sets: %d\n", n_mis))
    }
    n_sample <- sum(fs$mis_method %in% "exact_sample")
    if (n_sample > 0) {
      cat(sprintf("Families where the criteria were evaluated on a sample of maximum sets: %d\n", n_sample))
    }
  } else if (!x$diagnostics) {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); ",
            "rerun with diagnostics = TRUE for a richer output.")
  }
  invisible(x)
}

#' Internal helper: a number given by its base-10 logarithm, as text
#'
#' Formats a count known only through its logarithm, since counts of sets of
#' individuals easily exceed the largest double: in full below a million, in
#' scientific notation above.
#'
#' @keywords internal
#' @noRd
.rkn_format_log10 <- function(total) {
  if (total < 6) return(formatC(round(10^total), format = "f", digits = 0, big.mark = ","))
  exponent <- floor(total)
  mantissa <- round(10^(total - exponent), 1)
  if (mantissa >= 10) {
    mantissa <- 1
    exponent <- exponent + 1
  }
  sprintf("%.1fe+%02d", mantissa, exponent)
}

#' Summarise an rkn_selection object
#'
#' Computes summary tables of a selection of unrelated individuals: the selection by
#' family structure and, when diagnostics were computed, how each individual was
#' decided. The family structures and the types of decision are described in the
#' sections Family structures and Types of decision of \code{\link{rkn_selection}}.
#'
#' @param object An \code{rkn_selection} object, as returned by
#'   \code{\link{rkn_select_random}} or \code{\link{rkn_select_criteria}}.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.rkn_selection}: a list with elements
#' \describe{
#'   \item{\code{header}}{a one-line description of how the selection was made;}
#'   \item{\code{overall}}{a one-row data frame with the numbers of individuals, of
#'     selected individuals and of families;}
#'   \item{\code{by_structure}}{a data frame with one row per family structure present,
#'     giving the family size (a range for families of 6 or more members), the
#'     structure, and the numbers of families, of individuals and of selected
#'     individuals;}
#'   \item{\code{decisions}}{a contingency table of individuals by type of decision and
#'     outcome (selected or discarded), or \code{NULL} if diagnostics were not
#'     computed.}
#' }
#'
#' @seealso \code{\link{rkn_selection}}
#'
#' @example man-roxygen/example-network.R
#' @examples
#' sel <- rkn_select_random(g, seed = 1)
#' summary(sel)
#'
#' # The tables as data
#' s <- summary(sel)
#' s$by_structure
#' s$decisions
#'
#' @export
summary.rkn_selection <- function(object, ...) {
  g <- object$graph
  comp <- igraph::components(g)
  size <- comp$csize
  n_edges <- tabulate(comp$membership[igraph::ends(g, igraph::E(g), names = FALSE)[, 1]],
                      nbins = comp$no)
  kept_per_family <- tabulate(comp$membership[match(object$kept, igraph::V(g)$name)],
                              nbins = comp$no)

  # Index into .rkn_structure_classes: sizes 1 and 2, then three individuals as a chain
  # (two relationships) or all related (three), then 4, 5 and 6 or more members.
  class_id <- ifelse(size <= 2, size, pmin(size, 6) + 1)
  class_id[size == 3 & n_edges == 2] <- 3

  present <- sort(unique(class_id))
  family_size <- character(length(present))
  n_families <- integer(length(present))
  n_individuals <- integer(length(present))
  n_selected <- integer(length(present))
  for (k in seq_along(present)) {
    in_class <- class_id == present[k]
    r <- range(size[in_class])
    family_size[k] <- if (r[1] == r[2]) as.character(r[1]) else paste0(r[1], "-", r[2])
    n_families[k] <- sum(in_class)
    n_individuals[k] <- sum(size[in_class])
    n_selected[k] <- sum(kept_per_family[in_class])
  }
  by_structure <- .rkn_df(
    family_size = family_size,
    structure = .rkn_structure_classes[present],
    n_families = n_families,
    n_individuals = n_individuals,
    n_selected = n_selected
  )

  overall <- .rkn_df(
    n_individuals = length(object$kept) + length(object$discarded),
    n_selected = length(object$kept),
    n_families = comp$no
  )

  decisions <- NULL
  if (object$diagnostics && !identical(object$method, "global")) {
    nd <- object$node_diagnostics
    labels <- .rkn_labels_of(object)
    if ((isTRUE(object$uniform) || isTRUE(object$optimal)) && !any(nd$pre_class == "candidate_simple")) {
      labels <- labels[names(labels) != "candidate_simple"]
    }
    decisions <- table(
      decision = factor(labels[nd$pre_class], levels = labels),
      outcome = factor(ifelse(nd$post_class == "kept", "selected", "discarded"),
                       levels = c("selected", "discarded"))
    )
  }

  structure(
    list(header = .rkn_selection_header(object), overall = overall,
         by_structure = by_structure, decisions = decisions, method = object$method),
    class = "summary.rkn_selection"
  )
}

#' Print a summary of an rkn_selection object
#'
#' @param x A \code{summary.rkn_selection} object, as returned by
#'   \code{\link{summary.rkn_selection}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @method print summary.rkn_selection
#' @export
print.summary.rkn_selection <- function(x, ...) {
  o <- x$overall
  cat(sprintf("<rkn_selection summary> %s\n", x$header))
  cat(sprintf("Individuals: %d | selected: %d (%.1f%%) | families: %d\n\n",
              o$n_individuals, o$n_selected, 100 * o$n_selected / o$n_individuals, o$n_families))

  # Short headers keep the table within 80 columns.
  by_structure <- x$by_structure
  names(by_structure) <- c("size", "structure", "families", "individuals", "selected")
  cat("Selection by family structure:\n")
  print(by_structure, row.names = FALSE)

  if (identical(x$method, "global")) return(invisible(x))
  if (is.null(x$decisions)) {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); ",
            "rerun with diagnostics = TRUE for a richer output.")
    return(invisible(x))
  }

  cat("\nIndividuals by type of decision and outcome:\n")
  print(x$decisions)
  invisible(x)
}

#' @param ... Not used.
#'
#' @rdname rkn_plot_pre
#' @order 1
#' @export
plot.rkn_selection <- function(x, type = c("post", "diagnostics", "pre"),
                               layout = c("auto", "components", "grid"), seed = 1,
                               palette = c("colorblind", "classic"), show_outcome = TRUE,
                               colour_by = c("size", "family"), min_family_size = 1,
                               vertex_size = NULL, edge_width = 1, vertex_label = FALSE,
                               show_legend = TRUE, main = NULL, ...) {
  type <- match.arg(type)
  layout <- match.arg(layout)
  palette <- match.arg(palette)
  colour_by <- match.arg(colour_by)
  if (identical(type, "diagnostics") && identical(x$method, "global")) {
    message("The diagnostics plot is not available for global selections; plotting selected vs ",
            "discarded instead.")
    type <- "post"
  }
  if (identical(type, "diagnostics") && !x$diagnostics) {
    message("Node-level diagnostics were not computed (diagnostics = FALSE); plotting ",
            "selected vs discarded instead. Rerun with diagnostics = TRUE for the diagnostics plot.")
    type <- "post"
  }

  common <- list(layout = layout, seed = seed, min_family_size = min_family_size,
                 vertex_size = vertex_size, edge_width = edge_width,
                 vertex_label = vertex_label, show_legend = show_legend, main = main)
  if (identical(type, "post")) {
    do.call(.rkn_plot_post, c(list(x$graph, kept = x$kept, palette = palette), common))
  } else if (identical(type, "diagnostics")) {
    do.call(.rkn_plot_diagnostics, c(list(x$graph, res = x, palette = palette,
                                          show_outcome = show_outcome), common))
  } else {
    do.call(rkn_plot_pre, c(list(x$graph, colour_by = colour_by), common))
  }
  invisible(x)
}
