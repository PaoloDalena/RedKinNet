#' Internal helper: greedy elimination ordering of a graph
#'
#' Eliminates the vertices one at a time: the chosen vertex is the one whose neighbours
#' need the fewest new edges to become a clique (min-fill, ties broken by the fewest
#' neighbours), or simply the one with the fewest neighbours (min-degree) when the graph
#' has more than \code{fill_limit} vertices, where computing the fill at every step is
#' slow. The neighbours of the chosen vertex are joined into a clique and the vertex is
#' removed. The neighbours of each vertex at its elimination form the scope of the
#' message it sends in \code{.rkn_mis_eliminate()}; the largest scope is the width of
#' the ordering, an upper bound on the treewidth of the graph.
#'
#' @param A Symmetric logical adjacency matrix without self-loops.
#' @return A list with \code{order}, the elimination order; \code{later}, the neighbours
#'   of each vertex at its elimination; and \code{width}.
#'
#' @keywords internal
#' @noRd
.rkn_elimination_order <- function(A, fill_limit = 1000) {
  n <- nrow(A)
  alive <- rep(TRUE, n)
  order <- integer(n)
  later <- vector("list", n)
  use_fill <- n <= fill_limit
  for (step in seq_len(n)) {
    idx <- which(alive)
    sub <- A[idx, idx, drop = FALSE]
    deg <- rowSums(sub)
    if (use_fill) {
      fill <- numeric(length(idx))
      for (k in seq_along(idx)) {
        nb <- which(sub[k, ])
        m <- length(nb)
        if (m > 1) fill[k] <- m * (m - 1) / 2 - sum(sub[nb, nb]) / 2
      }
      candidates <- which(fill == min(fill))
      pick <- candidates[which.min(deg[candidates])]
    } else {
      pick <- which.min(deg)
    }
    v <- idx[pick]
    nb <- idx[sub[pick, ]]
    later[[v]] <- nb
    A[nb, nb] <- TRUE
    A[cbind(nb, nb)] <- FALSE
    A[v, ] <- FALSE
    A[, v] <- FALSE
    alive[v] <- FALSE
    order[step] <- v
  }
  list(order = order, later = later, width = max(c(0L, lengths(later))))
}

#' Internal helper: maximum independent sets by variable elimination
#'
#' Eliminates the vertices of \code{A} in the order given by
#' \code{.rkn_elimination_order()}. Eliminating a vertex v combines the messages
#' addressed to it into a table over v and its neighbours at elimination: for each
#' assignment of these vertices (selected or not), the largest number of selected
#' vertices among v and the vertices eliminated before it, and the number of ways to
#' reach it. Assignments selecting v together with one of its relatives in \code{A} are
#' infeasible. Summing v out, keeping for each assignment of the other vertices the
#' larger of the two maxima and adding the counts that reach it, gives the message sent
#' to the first of these vertices to be eliminated. The last vertex of each connected
#' component sends no message: its result gives the maximum size and the number of
#' maximum independent sets of the component. Counts are rescaled when they grow very
#' large, and their total is kept on the log scale.
#'
#' @param A Symmetric logical adjacency matrix without self-loops.
#' @param elim The result of \code{.rkn_elimination_order(A)}.
#' @return A list with \code{alpha}, the maximum size; \code{log_count}, the natural
#'   logarithm of the number of maximum independent sets; and the messages needed by
#'   \code{.rkn_mis_sample()}.
#'
#' @keywords internal
#' @noRd
.rkn_mis_eliminate <- function(A, elim) {
  n <- nrow(A)
  order <- elim$order
  later <- elim$later
  pos <- integer(n)
  pos[order] <- seq_len(n)
  msg <- vector("list", n)
  inbox <- vector("list", n)
  alpha <- 0
  log_count <- 0

  for (v in order) {
    scope <- later[[v]]
    code <- 0:(2^(length(scope) + 1) - 1)
    # Bit 0 of an assignment code is v, bit j is scope[j].
    bit <- function(j) (code %/% 2^j) %% 2
    xv <- bit(0)
    S <- xv
    C <- rep(1, length(code))
    relatives <- scope[A[v, scope]]
    if (length(relatives) > 0) {
      clash <- rep(FALSE, length(code))
      for (u in relatives) clash <- clash | bit(match(u, scope)) == 1
      clash <- clash & xv == 1
      S[clash] <- -Inf
      C[clash] <- 0
    }
    for (w in inbox[[v]]) {
      m <- msg[[w]]
      index <- 0
      for (k in seq_along(m$vars)) {
        j <- if (m$vars[k] == v) 0 else match(m$vars[k], scope)
        index <- index + bit(j) * 2^(k - 1)
      }
      S <- S + m$s[index + 1]
      C <- C * m$c[index + 1]
    }

    s0 <- S[xv == 0]
    s1 <- S[xv == 1]
    best <- pmax(s0, s1)
    count <- (s0 == best) * C[xv == 0] + (s1 == best) * C[xv == 1]
    top <- max(count)
    if (top > 1e100) {
      count <- count / top
      log_count <- log_count + log(top)
    }
    if (length(scope) == 0) {
      alpha <- alpha + best
      log_count <- log_count + log(count)
    } else {
      msg[[v]] <- list(vars = scope, s = best, c = count)
      target <- scope[which.min(pos[scope])]
      inbox[[target]] <- c(inbox[[target]], v)
    }
  }
  list(alpha = alpha, log_count = log_count, A = A, order = order, later = later, msg = msg, inbox = inbox)
}

#' Internal helper: draw a maximum independent set uniformly at random
#'
#' Goes through the vertices in reverse elimination order. When a vertex is reached, the
#' vertices of its scope are already decided; the messages addressed to it give, for
#' each choice of the vertex, the largest size reachable by the vertices eliminated
#' before it and the number of ways to reach it. The vertex is selected or not with
#' probability proportional to the number of ways among the choices that reach the
#' largest size, so that every maximum independent set has the same probability.
#'
#' @param sol The result of \code{.rkn_mis_eliminate()}.
#' @return A logical vector, \code{TRUE} for the selected vertices.
#'
#' @keywords internal
#' @noRd
.rkn_mis_sample <- function(sol) {
  x <- integer(nrow(sol$A))
  for (v in rev(sol$order)) {
    scope <- sol$later[[v]]
    score <- c(0, 1)
    count <- c(1, 1)
    relatives <- scope[sol$A[v, scope]]
    if (length(relatives) > 0 && any(x[relatives] == 1)) {
      score[2] <- -Inf
      count[2] <- 0
    }
    for (w in sol$inbox[[v]]) {
      m <- sol$msg[[w]]
      values <- x[m$vars]
      is_v <- m$vars == v
      weights <- 2^(seq_along(values) - 1)
      for (choice in 0:1) {
        values[is_v] <- choice
        i <- sum(values * weights) + 1
        score[choice + 1] <- score[choice + 1] + m$s[i]
        count[choice + 1] <- count[choice + 1] * m$c[i]
      }
    }
    p <- (score == max(score)) * count
    x[v] <- if (p[2] == 0) 0L else if (p[1] == 0) 1L else as.integer(stats::runif(1) < p[2] / sum(p))
  }
  x == 1
}

#' Internal helper: exact maximum independent sets of a graph, within a width limit
#'
#' @param A Symmetric logical adjacency matrix without self-loops.
#' @param max_width Largest width of the elimination ordering accepted. Time and memory
#'   grow as 2 to the power of the width: at width 20 the largest table holds about two
#'   million values (some 32 MB), which keeps the computation within the memory of an
#'   ordinary laptop. Above the limit the computation is not attempted.
#' @return The result of \code{.rkn_mis_eliminate()} with \code{feasible = TRUE} and
#'   \code{width}, or \code{list(feasible = FALSE, width = ...)}.
#'
#' @keywords internal
#' @noRd
.rkn_mis_exact <- function(A, max_width = 20) {
  elim <- .rkn_elimination_order(A)
  if (elim$width > max_width) return(list(feasible = FALSE, width = elim$width))
  sol <- .rkn_mis_eliminate(A, elim)
  sol$feasible <- TRUE
  sol$width <- elim$width
  sol
}
