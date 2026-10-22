#' Check samples for duplicates, possible contamination and different ancestry
#'
#' @description
#' Summarises, for each sample, its relatives and its mean kinship with the rest of the
#' cohort, and flags three kinds of anomalies before a selection:
#' \itemize{
#'   \item \strong{possible duplicates or monozygotic twins}: pairs with a kinship of at
#'     least 0.354. They need no action for the selection, which keeps at most one of
#'     them, but they are worth a quality check;
#'   \item \strong{possible contamination}: samples related to unusually many others
#'     with an unusually high mean kinship with everyone. A sample containing DNA from
#'     more than one person tends to appear slightly related to everyone, and can merge
#'     many families into one; excluding it (argument \code{exclude} of
#'     \code{\link{rkn_graph_from_kinship}}) restores the real families;
#'   \item \strong{possible different ancestry}: samples with an unusually low mean
#'     kinship. KING estimates are negative for pairs of different ancestry
#'     (Manichaikul et al., 2010). These samples do not affect the selection, but matter
#'     for downstream association analyses.
#' }
#' Nothing is removed: the decision is left to the user.
#'
#' @inheritParams rkn_graph_from_kinship
#' @param z_kinship,z_relatives A sample is flagged as possible contamination when the
#'   robust z-score of its mean kinship exceeds \code{z_kinship} and that of its number
#'   of relatives exceeds \code{z_relatives}.
#' @param z_relatives_only When the mean kinship cannot be computed, a sample is flagged
#'   as possible contamination when the robust z-score of its number of relatives
#'   exceeds \code{z_relatives_only}.
#' @param z_ancestry A sample is flagged as possible different ancestry when the robust
#'   z-score of its mean kinship is below \code{z_ancestry}.
#'
#' @details
#' The robust z-score of a value is its distance from the median of the cohort in units
#' of the median absolute deviation (MAD), so that the anomalies being sought do not
#' inflate the scale. When the MAD is zero, as for the number of relatives in a cohort
#' where most samples have none, the mean absolute deviation from the median is used
#' instead, rescaled to be comparable with the MAD. The number of relatives enters on the
#' logarithmic scale, \code{log(1 + relatives)}.
#'
#' The mean kinship of a sample is computed over all the other samples, so it requires a
#' kinship matrix or a pair table that lists every pair. With a table of related pairs
#' only, it is not available: contamination is then judged on the number of relatives
#' alone (see \code{z_relatives_only}) and ancestry is not assessed.
#'
#' Relatives are counted by degree with the KING cutoffs listed in
#' \code{\link{rkn_graph_from_kinship}}.
#'
#' @return An object of class \code{rkn_sample_check}, with a \code{print} method: a
#'   list with elements
#' \describe{
#'   \item{\code{samples}}{a data frame with one row per sample: \code{id};
#'     \code{n_relatives}, the number of samples with a kinship of at least
#'     \code{threshold}; \code{n_duplicate}, \code{n_first}, \code{n_second} and
#'     \code{n_third}, the relatives by degree; \code{mean_kinship}; the robust z-scores
#'     \code{z_relatives} and \code{z_kinship}; \code{flag}, \code{"possible
#'     contamination"}, \code{"possible different ancestry"} or an empty string;}
#'   \item{\code{duplicates}}{a data frame of the pairs with a kinship of at least 0.354;}
#'   \item{\code{contamination}}{a data frame of the samples flagged as possible
#'     contamination, with \code{largest_component_without}, the size of the largest
#'     family of the network after removing the sample;}
#'   \item{\code{threshold}, \code{largest_component}}{the threshold used and the size
#'     of the largest family of the network;}
#'   \item{\code{mean_kinship_available}}{whether the mean kinship could be computed.}
#' }
#'
#' @references
#' Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010). Robust
#' relationship inference in genome-wide association studies. \emph{Bioinformatics},
#' 26(22), 2867--2873. \doi{10.1093/bioinformatics/btq559}
#'
#' @seealso \code{\link{rkn_graph_from_kinship}}
#'
#' @examples
#' # A synthetic cohort of 200 samples: 30 parent-child pairs, one sample genotyped
#' # twice (S199 and S200), one contaminated sample (S001) and 15 samples of a
#' # different ancestry (S181 to S195).
#' set.seed(1)
#' ids <- sprintf("S%03d", 1:200)
#' K <- matrix(rnorm(200^2, -0.015, 0.004), 200, 200, dimnames = list(ids, ids))
#' K[upper.tri(K)] <- t(K)[upper.tri(K)]
#' ancestry <- 181:195
#' K[ancestry, -ancestry] <- rnorm(15 * 185, -0.15, 0.03)
#' K[-ancestry, ancestry] <- t(K[ancestry, -ancestry])
#' for (i in seq(2, 60, by = 2)) K[i, i + 1] <- K[i + 1, i] <- 0.25
#' K["S199", "S200"] <- K["S200", "S199"] <- 0.49
#' K["S001", 62:160] <- K[62:160, "S001"] <- runif(99, 0.045, 0.09)
#' diag(K) <- 0.5
#'
#' chk <- rkn_check_samples(K, threshold = "third")
#' chk
#' head(chk$samples)
#'
#' # Build the network without the possibly contaminated sample
#' g <- rkn_graph_from_kinship(K, threshold = "third", exclude = chk$contamination$id)
#'
#' @export
rkn_check_samples <- function(x, threshold, samples = NULL, z_kinship = 5, z_relatives = 3,
                              z_relatives_only = 5, z_ancestry = -5) {
  threshold <- .rkn_resolve_threshold(threshold)
  input <- .rkn_read_kinship(x, samples)
  ids <- input$ids
  n <- length(ids)
  lower <- c(.rkn_king_cutoffs[c("duplicate", "first", "second", "third")])
  upper <- c(Inf, .rkn_king_cutoffs[c("duplicate", "first", "second")])

  n_class <- matrix(0L, n, 4)
  if (!is.null(input$matrix)) {
    K <- input$matrix
    diag(K) <- NA
    n_relatives <- as.integer(rowSums(K >= threshold, na.rm = TRUE))
    for (j in 1:4) n_class[, j] <- as.integer(rowSums(K >= lower[j] & K < upper[j], na.rm = TRUE))
    mean_kinship <- unname(rowMeans(K, na.rm = TRUE))
  } else {
    a <- match(input$id1, ids)
    b <- match(input$id2, ids)
    k <- input$kinship
    ok <- !is.na(k) & a != b
    count_ends <- function(sel) tabulate(c(a[sel], b[sel]), nbins = n)
    n_relatives <- count_ends(ok & k >= threshold)
    for (j in 1:4) n_class[, j] <- count_ends(ok & k >= lower[j] & k < upper[j])
    mean_kinship <- rep(NA_real_, n)
    if (input$complete) {
      sums <- rowsum(c(k[ok], k[ok]), c(a[ok], b[ok]), reorder = TRUE)
      mean_kinship[as.integer(rownames(sums))] <- sums[, 1] / count_ends(ok)[as.integer(rownames(sums))]
    }
  }
  available <- !all(is.na(mean_kinship))

  z_rel <- .rkn_robust_z(log1p(n_relatives))
  z_kin <- if (available) .rkn_robust_z(mean_kinship) else rep(NA_real_, n)
  contaminated <- if (available) z_kin > z_kinship & z_rel > z_relatives else z_rel > z_relatives_only
  ancestry <- available & !contaminated & z_kin < z_ancestry
  flag <- rep("", n)
  flag[contaminated %in% TRUE] <- "possible contamination"
  flag[ancestry %in% TRUE] <- "possible different ancestry"

  samples_df <- data.frame(
    id = ids, n_relatives = n_relatives,
    n_duplicate = n_class[, 1], n_first = n_class[, 2], n_second = n_class[, 3], n_third = n_class[, 4],
    mean_kinship = mean_kinship, z_relatives = z_rel, z_kinship = z_kin, flag = flag,
    stringsAsFactors = FALSE
  )

  dup <- .rkn_kinship_edges(input, .rkn_king_cutoffs[["duplicate"]])
  duplicates <- data.frame(id1 = dup$id1, id2 = dup$id2, kinship = dup$kinship, stringsAsFactors = FALSE)

  edges <- .rkn_kinship_edges(input, threshold)
  g <- igraph::graph_from_data_frame(data.frame(from = edges$id1, to = edges$id2),
                                     directed = FALSE, vertices = data.frame(name = ids))
  flagged <- which(flag == "possible contamination")
  without <- integer(length(flagged))
  for (k in seq_along(flagged)) {
    without[k] <- as.integer(max(igraph::components(igraph::delete_vertices(g, flagged[k]))$csize))
  }
  contamination <- samples_df[flagged, c("id", "n_relatives", "mean_kinship", "z_relatives", "z_kinship")]
  contamination$largest_component_without <- without
  rownames(contamination) <- NULL

  structure(
    list(samples = samples_df, duplicates = duplicates, contamination = contamination,
         threshold = threshold, largest_component = as.integer(max(igraph::components(g)$csize)),
         mean_kinship_available = available),
    class = "rkn_sample_check"
  )
}

#' Internal helper: robust z-scores
#'
#' Distance from the median in units of the MAD, or of the mean absolute deviation from
#' the median (rescaled by sqrt(pi / 2), which makes it consistent with the standard
#' deviation under normality, as the MAD is) when the MAD is zero; all zeros when the
#' values do not vary.
#'
#' @keywords internal
#' @noRd
.rkn_robust_z <- function(x) {
  centre <- stats::median(x, na.rm = TRUE)
  scale <- stats::mad(x, na.rm = TRUE)
  if (scale == 0) scale <- sqrt(pi / 2) * mean(abs(x - centre), na.rm = TRUE)
  if (scale == 0) return(rep(0, length(x)))
  (x - centre) / scale
}

#' Print an rkn_sample_check object
#'
#' @param x An \code{rkn_sample_check} object, as returned by
#'   \code{\link{rkn_check_samples}}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.rkn_sample_check <- function(x, ...) {
  s <- x$samples
  cat(sprintf("<rkn_sample_check> %d samples, relatedness threshold %g\n\n", nrow(s), x$threshold))

  cat(sprintf("Possible duplicates or monozygotic twins (kinship >= 0.354): %d pair(s)\n", nrow(x$duplicates)))
  if (nrow(x$duplicates) > 0) print(x$duplicates, row.names = FALSE, digits = 3)

  cat(sprintf("\nPossible contamination: %d sample(s)\n(related to unusually many samples%s)\n",
              nrow(x$contamination),
              if (x$mean_kinship_available) ", with an unusually high mean kinship" else ""))
  if (nrow(x$contamination) > 0) {
    print(x$contamination, row.names = FALSE, digits = 3)
    cat(sprintf("Largest family of the network: %d samples.\n", x$largest_component))
  }

  if (x$mean_kinship_available) {
    cat(sprintf("\nPossible different ancestry (mean kinship far below the cohort): %d sample(s)\n",
                sum(s$flag == "possible different ancestry")))
    cat("These do not affect the selection; see $samples for the values.\n")
  } else {
    cat("\nAncestry not assessed: the mean kinship needs every pair of samples.\n")
  }
  invisible(x)
}
