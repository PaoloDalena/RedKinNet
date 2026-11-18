.rkn_king_cutoffs <- c(duplicate = 0.354, first = 0.177, second = 0.0884, third = 0.0442)

#' Internal helper: numeric kinship threshold from a number or a KING label
#'
#' @keywords internal
#' @noRd
.rkn_resolve_threshold <- function(threshold) {
  labels <- names(.rkn_king_cutoffs)
  if (is.character(threshold) && length(threshold) == 1 && !is.na(threshold) && threshold %in% labels) {
    threshold <- .rkn_king_cutoffs[[threshold]]
  }
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold)) {
    stop("`threshold` must be a number or one of \"duplicate\", \"first\", \"second\" and \"third\".",
         call. = FALSE)
  }
  threshold
}

#' Internal helper: validate kinship estimates given as a matrix or a pair table
#'
#' Returns the sample identifiers and either the matrix (\code{matrix}) or the listed
#' pairs (\code{id1}, \code{id2}, \code{kinship}), with \code{complete} telling whether
#' every pair of samples is listed. Pairs with a missing estimate raise a warning.
#'
#' @keywords internal
#' @noRd
.rkn_read_kinship <- function(x, samples) {
  if (is.matrix(x)) {
    ids <- rownames(x)
    if (!is.numeric(x) || is.null(ids) || !identical(ids, colnames(x))) {
      stop("A kinship matrix must be numeric, with the sample identifiers as identical row and column names.",
           call. = FALSE)
    }
    if (!isSymmetric(unname(x))) stop("A kinship matrix must be symmetric.", call. = FALSE)
    out <- list(ids = ids, matrix = x, complete = TRUE)
    n_missing <- (sum(is.na(x)) - sum(is.na(diag(x)))) / 2
  } else if (is.data.frame(x)) {
    nm <- tolower(names(x))
    cols <- if (all(c("iid1", "iid2") %in% nm)) c("iid1", "iid2") else c("id1", "id2")
    if (!all(c(cols, "kinship") %in% nm)) {
      stop("A pair table needs two identifier columns (ID1 and ID2, IID1 and IID2, or id1 and id2) ",
           "and a kinship column.", call. = FALSE)
    }
    a <- as.character(x[[match(cols[1], nm)]])
    b <- as.character(x[[match(cols[2], nm)]])
    k <- x[[match("kinship", nm)]]
    listed <- unique(c(a, b))
    ids <- if (is.null(samples)) listed else unique(as.character(samples))
    absent <- setdiff(listed, ids)
    if (length(absent) > 0) {
      stop("The pair table refers to samples missing from `samples`, for instance ",
           paste(absent[seq_len(min(3, length(absent)))], collapse = ", "), ".", call. = FALSE)
    }
    n_pairs <- nrow(unique(data.frame(lo = pmin(a, b), hi = pmax(a, b))))
    complete <- n_pairs >= length(ids) * (length(ids) - 1) / 2
    if (is.null(samples) && !complete) {
      message("The pair table does not list every pair of samples, so samples without relatives ",
              "may be missing from the network; supply all sample identifiers through `samples`.")
    }
    out <- list(ids = ids, id1 = a, id2 = b, kinship = k, complete = complete)
    n_missing <- sum(is.na(k))
  } else {
    stop("`x` must be a kinship matrix or a data frame of pairs.", call. = FALSE)
  }
  if (n_missing > 0) {
    warning(sprintf(if (n_missing == 1) "%d pair has a missing kinship estimate and is treated as unrelated." else
                      "%d pairs have a missing kinship estimate and are treated as unrelated.", n_missing),
            call. = FALSE)
  }
  out
}

#' Internal helper: remove samples from kinship input read by .rkn_read_kinship()
#'
#' @keywords internal
#' @noRd
.rkn_exclude_samples <- function(input, exclude) {
  keep <- !input$ids %in% exclude
  if (!is.null(input$matrix)) {
    input$matrix <- input$matrix[keep, keep, drop = FALSE]
  } else {
    drop <- input$id1 %in% exclude | input$id2 %in% exclude
    input$id1 <- input$id1[!drop]
    input$id2 <- input$id2[!drop]
    input$kinship <- input$kinship[!drop]
  }
  input$ids <- input$ids[keep]
  input
}

#' Internal helper: one warning summarising the anomalies found by the sample check
#'
#' Returns \code{NULL} when the check found nothing.
#'
#' @keywords internal
#' @noRd
.rkn_check_warning <- function(chk) {
  count <- function(n, what) sprintf("%d %s%s", n, what, if (n > 1) "s" else "")
  parts <- character(0)
  ids <- chk$contamination$id
  if (length(ids) > 0) {
    shown <- paste(ids[seq_len(min(5, length(ids)))], collapse = ", ")
    if (length(ids) > 5) shown <- paste0(shown, ", ...")
    parts <- c(parts, sprintf("%s flagged as possible contamination (%s)", count(length(ids), "sample"), shown))
  }
  if (nrow(chk$duplicates) > 0) {
    parts <- c(parts, paste(count(nrow(chk$duplicates), "pair"), "of possible duplicates or monozygotic twins"))
  }
  if (length(parts) == 0) return(NULL)
  paste0("Sample check: ", paste(parts, collapse = ", "), ". See rkn_check_samples() for details ",
         "and the argument `exclude` to leave samples out; check = FALSE skips this check.")
}

#' Internal helper: pairs of samples whose kinship is at least the threshold
#'
#' @keywords internal
#' @noRd
.rkn_kinship_edges <- function(input, threshold) {
  if (!is.null(input$matrix)) {
    hit <- which(input$matrix >= threshold, arr.ind = TRUE)
    hit <- hit[hit[, 1] < hit[, 2], , drop = FALSE]
    return(list(id1 = input$ids[hit[, 1]], id2 = input$ids[hit[, 2]], kinship = input$matrix[hit]))
  }
  keep <- !is.na(input$kinship) & input$kinship >= threshold
  list(id1 = input$id1[keep], id2 = input$id2[keep], kinship = input$kinship[keep])
}

#' Build a relatedness network from kinship estimates
#'
#' @description
#' Builds the relatedness network used by \code{\link{rkn_select_random}} and
#' \code{\link{rkn_select_criteria}}: one vertex per sample, and an edge between two
#' samples whose kinship coefficient is at least \code{threshold}. The families of the
#' selectors are the connected components of this network.
#'
#' @param x Kinship estimates, either a square symmetric matrix with the sample
#'   identifiers as row and column names (the diagonal is ignored), or a data frame with
#'   one row per pair of samples, such as the \code{.kin0} output of KING or the output of
#'   \code{plink2 --make-king-table}. In a data frame, the two identifier columns are
#'   recognised by the names \code{ID1} and \code{ID2} (KING), \code{IID1} and
#'   \code{IID2} (PLINK 2) or \code{id1} and \code{id2}, and the kinship column by the
#'   name \code{kinship}, regardless of case.
#' @param threshold The kinship coefficient from which two samples are considered
#'   related: a number, or one of \code{"duplicate"}, \code{"first"}, \code{"second"} and
#'   \code{"third"}, which select the cutoff of KING for relatives of that degree or
#'   closer (see Details). There is no default: the choice depends on the study.
#' @param samples Optional character vector with the identifiers of all genotyped
#'   samples, for instance the IID column of the PLINK \code{.fam} or \code{.psam} file.
#'   Pair tables usually list only related pairs, so samples without relatives appear in
#'   no row; \code{samples} adds them to the network as unrelated individuals. Ignored
#'   for a matrix.
#' @param exclude Optional character vector of sample identifiers to leave out of the
#'   network, for instance samples flagged by a quality check.
#' @param covariates Optional data frame of sample covariates, with the sample
#'   identifiers in the column named by \code{covariates_id}. Its other columns become
#'   vertex attributes of the network, which \code{\link{rkn_select_criteria}} reads;
#'   samples missing from \code{covariates} get \code{NA}. Covariates cannot be named
#'   \code{name}, which holds the sample identifiers, nor \code{degree},
#'   \code{max_path_length} or \code{.name}, which are reserved for the structural
#'   measures computed during the selection.
#' @param covariates_id Name of the identifier column of \code{covariates}.
#' @param check Logical. If \code{TRUE} (default), the samples left after
#'   \code{exclude} are checked with \code{\link{rkn_check_samples}}, with the same
#'   threshold and default settings, and a single warning summarises the possible
#'   duplicates and contamination it finds. Nothing is removed. Samples of possibly
#'   different ancestry, which do not affect the selection, are reported only by
#'   \code{rkn_check_samples()}.
#'   \code{check = FALSE} skips the check, for instance when rebuilding a network whose
#'   samples were already checked.
#'
#' @details
#' The threshold labels select the boundaries between the degrees of relationship
#' inferred by KING (Manichaikul et al., 2010):
#' \tabular{lll}{
#'   \strong{label} \tab \strong{cutoff} \tab \strong{relatives joined, for instance} \cr
#'   \code{"duplicate"} \tab 0.354 \tab the same individual genotyped twice, monozygotic twins \cr
#'   \code{"first"} \tab 0.177 \tab parent and child, full siblings, and the above \cr
#'   \code{"second"} \tab 0.0884 \tab grandparent and grandchild, aunt or uncle and nephew
#'     or niece, half-siblings, and the above \cr
#'   \code{"third"} \tab 0.0442 \tab first cousins, and the above
#' }
#' Pairs with a missing kinship estimate are treated as unrelated, with a warning.
#'
#' The sample check runs here, rather than in the selectors, because it needs the kinship
#' of every pair of samples, including the unrelated ones, which the network no longer
#' holds. A network built in another way can be checked by calling
#' \code{\link{rkn_check_samples}} on the kinship estimates it was built from.
#'
#' @return An undirected igraph object with one vertex per sample: the vertex attribute
#'   \code{name} holds the sample identifiers and the other vertex attributes the
#'   covariates, the edge attribute \code{kinship} holds the kinship coefficient, and
#'   the graph attribute \code{kinship_threshold} holds the cutoff used.
#'
#' @references
#' Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010). Robust
#' relationship inference in genome-wide association studies. \emph{Bioinformatics},
#' 26(22), 2867--2873. \doi{10.1093/bioinformatics/btq559}
#'
#' @seealso \code{\link{rkn_check_samples}}, \code{\link{rkn_select_random}},
#'   \code{\link{rkn_select_criteria}}, \code{\link{rkn_plot_pre}}.
#'
#' @examples
#' # Kinship estimates for six samples: a parent and child (A, B), first cousins (C, D)
#' # and two unrelated samples (E, F).
#' K <- matrix(-0.01, 6, 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
#' K["A", "B"] <- K["B", "A"] <- 0.25
#' K["C", "D"] <- K["D", "C"] <- 0.06
#' diag(K) <- 0.5
#'
#' # Relatives of the third degree or closer: two families of two, two singletons
#' g <- rkn_graph_from_kinship(K, threshold = "third")
#' igraph::components(g)$csize
#'
#' # First-degree relatives only: the cousins count as unrelated
#' igraph::components(rkn_graph_from_kinship(K, threshold = "first"))$csize
#'
#' # The same data as a pair table that lists only the related pairs, as KING does;
#' # samples adds the unrelated E and F
#' pairs <- data.frame(ID1 = c("A", "C"), ID2 = c("B", "D"), Kinship = c(0.25, 0.06))
#' g <- rkn_graph_from_kinship(pairs, threshold = "third", samples = LETTERS[1:6])
#'
#' # Covariates for a criteria-based selection
#' covariates <- data.frame(id = LETTERS[1:6], sex = c("M", "F", "F", "M", "M", "F"))
#' g <- rkn_graph_from_kinship(K, threshold = "third", covariates = covariates)
#' rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1)$kept
#'
#' @export
rkn_graph_from_kinship <- function(x, threshold, samples = NULL, exclude = NULL,
                                   covariates = NULL, covariates_id = "id", check = TRUE) {
  threshold <- .rkn_resolve_threshold(threshold)
  .rkn_check_flag(check, "check")
  input <- .rkn_read_kinship(x, samples)

  if (!is.null(exclude)) {
    unknown <- setdiff(exclude, input$ids)
    if (length(unknown) > 0) {
      stop("`exclude` contains identifiers absent from the data, for instance ",
           paste(unknown[seq_len(min(3, length(unknown)))], collapse = ", "), ".", call. = FALSE)
    }
    input <- .rkn_exclude_samples(input, exclude)
  }

  if (check) {
    msg <- .rkn_check_warning(.rkn_check_input(input, threshold))
    if (!is.null(msg)) warning(msg, call. = FALSE)
  }

  edges <- .rkn_kinship_edges(input, threshold)
  ids <- input$ids
  id1 <- edges$id1
  id2 <- edges$id2
  kinship <- edges$kinship

  vertices <- data.frame(name = ids, stringsAsFactors = FALSE)
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates) || !covariates_id %in% names(covariates)) {
      stop("`covariates` must be a data frame with an identifier column named \"", covariates_id, "\".",
           call. = FALSE)
    }
    cid <- as.character(covariates[[covariates_id]])
    if (anyDuplicated(cid) > 0) stop("`covariates` contains duplicated identifiers.", call. = FALSE)
    vars <- setdiff(names(covariates), covariates_id)
    if ("name" %in% vars) {
      stop("A covariate cannot be called \"name\", which holds the sample identifiers.", call. = FALSE)
    }
    reserved <- intersect(vars, .rkn_reserved_attributes)
    if (length(reserved) > 0) stop(.rkn_reserved_message(reserved), call. = FALSE)
    row <- match(ids, cid)
    if (anyNA(row)) {
      message(sprintf("Covariates are missing for %d of %d samples and set to NA.", sum(is.na(row)), length(ids)))
    }
    for (v in vars) vertices[[v]] <- covariates[[v]][row]
  }

  edges <- data.frame(from = id1, to = id2, kinship = kinship, stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  g <- igraph::simplify(g, edge.attr.comb = "max")
  g$kinship_threshold <- threshold
  g
}
