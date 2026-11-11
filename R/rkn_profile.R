#' Internal helper: size of the tables of a family profile
#'
#' Computes from the elimination ordering alone, before anything is allocated, the
#' number of cells of the largest table built by \code{.rkn_family_profile()} and the
#' number of cell operations of its products, a proxy of the run time. The value range
#' of the table of a vertex grows with the values of the vertices eliminated before it.
#'
#' @param elim The result of \code{.rkn_elimination_order()}.
#' @param x Non-negative integer value of each vertex.
#' @param k Largest number of individuals below the maximum size.
#' @return A named numeric vector with \code{cells} and \code{work}.
#'
#' @keywords internal
#' @noRd
.rkn_profile_budget <- function(elim, x, k) {
  n <- length(elim$order)
  pos <- integer(n)
  pos[elim$order] <- seq_len(n)
  inbox <- vector("list", n)
  r_sub <- numeric(n)
  cells <- 0
  work <- 0
  for (v in elim$order) {
    scope <- elim$later[[v]]
    n_code <- 2^(length(scope) + 1)
    r <- x[v]
    for (w in inbox[[v]]) {
      work <- work + n_code * (k + 1) * (k + 2) / 2 * (min(r, r_sub[w]) + 1) * (max(r, r_sub[w]) + 1)
      r <- r + r_sub[w]
    }
    r_sub[v] <- r
    cells <- max(cells, n_code * (k + 1) * (r + 1))
    if (length(scope) > 0) {
      target <- scope[which.min(pos[scope])]
      inbox[[target]] <- c(inbox[[target]], v)
    }
  }
  c(cells = cells, work = work)
}

#' Internal helper: product of two profile tables
#'
#' Multiplies two tables with one row per assignment, one column per number of
#' individuals below the largest size (0 to \code{k}) and one layer per value sum:
#' sizes and value sums add, so the product is a convolution over both, truncated at
#' \code{k}. The loop runs over the cells of the table with the shorter value range.
#'
#' @keywords internal
#' @noRd
.rkn_profile_product <- function(a, b, k) {
  if (dim(a)[3] > dim(b)[3]) {
    swap <- a
    a <- b
    b <- swap
  }
  n_b <- dim(b)[3]
  out <- array(0, c(dim(a)[1], k + 1, dim(a)[3] + n_b - 1))
  for (j in 0:k) {
    jj <- seq_len(k + 1 - j)
    for (r in seq_len(dim(a)[3])) {
      w <- a[, j + 1, r]
      if (!any(w > 0)) next
      cols <- r - 1 + seq_len(n_b)
      out[, j + jj, cols] <- out[, j + jj, cols, drop = FALSE] + w * b[, jj, , drop = FALSE]
    }
  }
  out
}

#' Internal helper: profile of a family by size and value
#'
#' For a connected family with adjacency matrix \code{A} and non-negative integer values
#' \code{x}, counts the independent sets of size \code{alpha - j}, for \code{j} from 0
#' to \code{k}, by the sum of the values of their members. The elimination follows
#' \code{.rkn_mis_eliminate()} with tables in place of numbers: for each assignment of a
#' vertex and its neighbours at elimination, the number of ways to complete it for each
#' size and each value sum. A partial set more than \code{k} below the largest size
#' reachable with the same assignment cannot belong to a set within \code{k} of the
#' maximum, so smaller sizes are dropped. Each table is rescaled to a largest count of
#' 1. The computation is not attempted when the tables would exceed \code{max_cells}
#' cells or their products \code{max_work} cell operations, and it stops when the
#' counts of a table span more than 300 natural orders of magnitude: beyond that, the
#' products of rare combinations could fall below double precision and be lost.
#'
#' @param A Symmetric logical adjacency matrix of a connected family.
#' @param x Non-negative integer value of each vertex.
#' @param k Largest number of individuals below the maximum size.
#' @param elim The result of \code{.rkn_elimination_order(A)}.
#' @param max_cells,max_work Budget of the computation.
#' @param keep If \code{TRUE}, the messages are returned as well, for
#'   \code{.rkn_profile_sample()}.
#' @return A list with \code{alpha}, the maximum size; \code{counts}, a matrix with one
#'   row per \code{j} and one column per value sum from 0, proportional to the numbers
#'   of sets; and \code{log_scale}, the logarithm of the factor that restores them. Or
#'   \code{NULL} when the family exceeds the budget.
#'
#' @keywords internal
#' @noRd
.rkn_family_profile <- function(A, x, k, elim = .rkn_elimination_order(A), max_cells = 4e6, max_work = 4e9,
                                keep = FALSE) {
  budget <- .rkn_profile_budget(elim, x, k)
  if (budget[["cells"]] > max_cells || budget[["work"]] > max_work) return(NULL)
  n <- nrow(A)
  pos <- integer(n)
  pos[elim$order] <- seq_len(n)
  msg <- vector("list", n)
  inbox <- vector("list", n)
  log_scale <- 0
  rescale <- function(tab) {
    top <- max(tab)
    log_scale <<- log_scale + log(top)
    tab <- tab / top
    if (-log(min(tab[tab > 0])) > 300) return(NULL)
    tab
  }

  for (v in elim$order) {
    scope <- elim$later[[v]]
    n_code <- 2^(length(scope) + 1)
    code <- 0:(n_code - 1)
    # Bit 0 of an assignment code is v, bit j is scope[j].
    bit <- function(j) (code %/% 2^j) %% 2
    xv <- bit(0)
    mx <- xv
    tab <- array(0, c(n_code, k + 1, x[v] + 1))
    tab[cbind(seq_len(n_code), 1, xv * x[v] + 1)] <- 1
    relatives <- scope[A[v, scope]]
    if (length(relatives) > 0) {
      clash <- rep(FALSE, n_code)
      for (u in relatives) clash <- clash | bit(match(u, scope)) == 1
      clash <- clash & xv == 1
      mx[clash] <- -Inf
      tab[clash, , ] <- 0
    }
    for (w in inbox[[v]]) {
      m <- msg[[w]]
      index <- 0
      for (i in seq_along(m$vars)) {
        j <- if (m$vars[i] == v) 0 else match(m$vars[i], scope)
        index <- index + bit(j) * 2^(i - 1)
      }
      tab <- rescale(.rkn_profile_product(tab, m$tab[index + 1, , , drop = FALSE], k))
      if (is.null(tab)) return(NULL)
      mx <- mx + m$mx[index + 1]
    }

    # Summing v out aligns the two choices of v on the larger maximum size.
    i0 <- which(xv == 0)
    i1 <- i0 + 1
    best <- pmax(mx[i0], mx[i1])
    out <- array(0, c(length(i0), k + 1, dim(tab)[3]))
    for (side in list(i0, i1)) {
      gap <- best - mx[side]
      for (g in 0:k) {
        rows <- which(gap == g)
        if (length(rows) == 0) next
        out[rows, (g + 1):(k + 1), ] <- out[rows, (g + 1):(k + 1), , drop = FALSE] +
          tab[side[rows], seq_len(k + 1 - g), , drop = FALSE]
      }
    }
    out <- rescale(out)
    if (is.null(out)) return(NULL)
    if (length(scope) == 0) {
      prof <- list(alpha = best, counts = matrix(out, k + 1, dim(out)[3]), log_scale = log_scale)
      if (keep) prof <- c(prof, list(A = A, x = x, k = k, elim = elim, msg = msg, inbox = inbox))
      return(prof)
    }
    msg[[v]] <- list(vars = scope, mx = best, tab = out)
    target <- scope[which.min(pos[scope])]
    inbox[[target]] <- c(inbox[[target]], v)
  }
}

#' Internal helper: draw an independent set of given size and value from a profile
#'
#' Draws uniformly at random one of the independent sets of size \code{alpha - j} and
#' value sum \code{r} counted by \code{.rkn_family_profile(keep = TRUE)}. Goes through
#' the vertices in reverse elimination order, as \code{.rkn_mis_sample()}: each vertex
#' is reached with its neighbours at elimination already decided and with the size and
#' value that it and the vertices eliminated before it must reach. Recomputing its
#' table for the two choices of the vertex gives the number of ways to reach the target
#' with each choice; the chosen one is then split among the messages the vertex
#' received, in proportion to the number of ways to reach each split, which sets the
#' targets of the vertices that sent them.
#'
#' @param prof The result of \code{.rkn_family_profile(keep = TRUE)}.
#' @param j Number of individuals below the maximum size.
#' @param r Value sum.
#' @return A logical vector, \code{TRUE} for the selected vertices.
#'
#' @keywords internal
#' @noRd
.rkn_profile_sample <- function(prof, j, r) {
  A <- prof$A
  x <- prof$x
  k <- prof$k
  later <- prof$elim$later
  n <- nrow(A)
  sel <- integer(n)
  tgt_j <- numeric(n)
  tgt_r <- numeric(n)
  root <- prof$elim$order[n]
  tgt_j[root] <- j
  tgt_r[root] <- r
  pick <- function(w) {
    w <- w / sum(w)
    sample.int(length(w), 1, prob = w)
  }

  for (v in rev(prof$elim$order)) {
    scope <- later[[v]]
    # Rows 1 and 2 of every table are the choices 0 and 1 of v.
    tab <- array(0, c(2, k + 1, x[v] + 1))
    tab[1, 1, 1] <- 1
    mx <- c(0, 1)
    relatives <- scope[A[v, scope]]
    if (length(relatives) > 0 && any(sel[relatives] == 1)) mx[2] <- -Inf else tab[2, 1, x[v] + 1] <- 1
    chain <- list(tab)
    rows <- list()
    for (w in prof$inbox[[v]]) {
      m <- prof$msg[[w]]
      idx <- integer(2)
      for (xv in 0:1) {
        bits <- sel[m$vars]
        bits[m$vars == v] <- xv
        idx[xv + 1] <- sum(bits * 2^(seq_along(bits) - 1)) + 1
      }
      tab <- .rkn_profile_product(tab, m$tab[idx, , , drop = FALSE], k)
      tab <- tab / max(tab)
      mx <- mx + m$mx[idx]
      chain[[length(chain) + 1]] <- tab
      rows[[length(rows) + 1]] <- list(w = w, idx = idx)
    }

    best <- max(mx)
    weight <- numeric(2)
    for (xv in 0:1) {
      jd <- tgt_j[v] - (best - mx[xv + 1])
      if (is.finite(jd) && jd >= 0 && jd <= k && tgt_r[v] < dim(tab)[3]) {
        weight[xv + 1] <- tab[xv + 1, jd + 1, tgt_r[v] + 1]
      }
    }
    xs <- pick(weight) - 1
    jc <- tgt_j[v] - (best - mx[xs + 1])
    rc <- tgt_r[v]

    # Split the target between the table before each message and the message.
    for (t in rev(seq_along(rows))) {
      prev <- matrix(chain[[t]][xs + 1, , ], k + 1)
      m_row <- matrix(prof$msg[[rows[[t]]$w]]$tab[rows[[t]]$idx[xs + 1], , ], k + 1)
      split <- matrix(0, jc + 1, ncol(m_row))
      for (j2 in 0:jc) {
        r2 <- 0:min(rc, ncol(m_row) - 1)
        ok <- rc - r2 < ncol(prev)
        r2 <- r2[ok]
        split[j2 + 1, r2 + 1] <- prev[jc - j2 + 1, rc - r2 + 1] * m_row[j2 + 1, r2 + 1]
      }
      cell <- pick(as.vector(split)) - 1
      j2 <- cell %% (jc + 1)
      r2 <- cell %/% (jc + 1)
      tgt_j[rows[[t]]$w] <- j2
      tgt_r[rows[[t]]$w] <- r2
      jc <- jc - j2
      rc <- rc - r2
    }
    sel[v] <- xs
  }
  sel == 1
}

#' Internal helper: partial profile of a family at its maximum size
#'
#' For a family beyond the budget of \code{.rkn_family_profile()}: the smallest and
#' largest value sums among its maximum independent sets, found exactly by the
#' weighted engine (integer weights that rank sets by size first and by value sum
#' second), and the value sums of \code{n_draws} uniform draws. A graph too wide for the
#' engine is reduced by weighted dominance first, as in
#' \code{.rkn_criteria_maximum_set()}.
#'
#' @param A Symmetric logical adjacency matrix of a connected family.
#' @param x Non-negative integer value of each vertex.
#' @param n_draws Number of uniform draws.
#' @param max_width Largest elimination width accepted, as in \code{.rkn_mis_exact()}.
#' @return A list with \code{alpha}; \code{values}, the distinct value sums found, in
#'   increasing order; and \code{sets}, for each value a logical vector marking a maximum
#'   independent set that reaches it. Or \code{NULL} when the family is too intricate
#'   even after the reduction.
#'
#' @keywords internal
#' @noRd
.rkn_partial_profile <- function(A, x, n_draws = 100, max_width = 20) {
  n <- nrow(A)
  solve <- function(weight) {
    keep <- rep(TRUE, n)
    sol <- .rkn_mis_exact(A, max_width, weight)
    if (!sol$feasible) {
      keep <- .rkn_dominance_reduce(A, weight)
      sol <- .rkn_mis_exact(A[keep, keep, drop = FALSE], max_width, weight[keep])
    }
    if (!sol$feasible) return(NULL)
    list(sol = sol, keep = keep)
  }
  full <- function(s, selected) {
    out <- rep(FALSE, n)
    out[which(s$keep)[selected]] <- TRUE
    out
  }
  base <- n * max(x) + 1
  plain <- solve(rep(1, n))
  largest <- solve(base + x)
  smallest <- solve(base + max(x) - x)
  if (is.null(plain) || is.null(largest) || is.null(smallest)) return(NULL)

  sets <- list(full(largest, .rkn_mis_sample(largest$sol)), full(smallest, .rkn_mis_sample(smallest$sol)))
  draws <- .rkn_mis_sample_many(plain$sol, n_draws)
  for (d in seq_len(n_draws)) sets[[length(sets) + 1]] <- full(plain, draws[, d])
  values <- numeric(length(sets))
  for (i in seq_along(sets)) values[i] <- sum(x[sets[[i]]])
  first <- !duplicated(values)
  o <- order(values[first])
  list(alpha = plain$sol$alpha, values = values[first][o], sets = sets[first][o])
}
