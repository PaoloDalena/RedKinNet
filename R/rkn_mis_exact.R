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
  # Rows and columns of eliminated vertices are cleared, so A always holds the graph of
  # the vertices still to eliminate; degree and fill are updated only for the vertices
  # whose neighbourhood changes at each step.
  fill_of <- function(x) {
    nb <- which(A[x, ])
    m <- length(nb)
    if (m < 2) 0 else m * (m - 1) / 2 - sum(A[nb, nb]) / 2
  }
  deg <- rowSums(A)
  fill <- numeric(n)
  if (use_fill) for (x in seq_len(n)) fill[x] <- fill_of(x)
  for (step in seq_len(n)) {
    candidates <- which(alive)
    if (use_fill) {
      candidates <- candidates[fill[candidates] == min(fill[candidates])]
    }
    v <- candidates[which.min(deg[candidates])]
    nb <- unname(which(A[v, ]))
    later[[v]] <- nb
    A[nb, nb] <- TRUE
    A[cbind(nb, nb)] <- FALSE
    A[v, ] <- FALSE
    A[, v] <- FALSE
    alive[v] <- FALSE
    order[step] <- v
    if (length(nb) > 0) {
      deg[nb] <- rowSums(A[nb, , drop = FALSE])
      if (use_fill) {
        # The fill changes only for the neighbours of v and for the vertices adjacent to
        # two of them, which may have gained an edge among their neighbours.
        touched <- unique(c(nb, which(colSums(A[nb, , drop = FALSE]) > 0)))
        for (x in touched) fill[x] <- fill_of(x)
      }
    }
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
#' large, and their total is kept on the log scale. With \code{weight}, the size of a set
#' is replaced by the sum of the weights of its vertices: integer weights such as
#' (n + 1) plus a 0/1 indicator rank sets by size first and by the indicator second.
#'
#' @param A Symmetric logical adjacency matrix without self-loops.
#' @param elim The result of \code{.rkn_elimination_order(A)}.
#' @param weight Weight of each vertex, 1 by default.
#' @return A list with \code{alpha}, the maximum size (or total weight);
#'   \code{log_count}, the natural logarithm of the number of sets that reach it; and
#'   the messages needed by \code{.rkn_mis_sample()}.
#'
#' @keywords internal
#' @noRd
.rkn_mis_eliminate <- function(A, elim, weight = rep(1, nrow(A))) {
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
    S <- xv * weight[v]
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
  list(alpha = alpha, log_count = log_count, A = A, weight = weight, order = order, later = later, msg = msg,
       inbox = inbox)
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
    score <- c(0, sol$weight[v])
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
#' @param weight Weight of each vertex, as in \code{.rkn_mis_eliminate()}.
#' @return The result of \code{.rkn_mis_eliminate()} with \code{feasible = TRUE} and
#'   \code{width}, or \code{list(feasible = FALSE, width = ...)}.
#'
#' @keywords internal
#' @noRd
.rkn_mis_exact <- function(A, max_width = 20, weight = rep(1, nrow(A))) {
  elim <- .rkn_elimination_order(A)
  if (elim$width > max_width) return(list(feasible = FALSE, width = elim$width))
  sol <- .rkn_mis_eliminate(A, elim, weight)
  sol$feasible <- TRUE
  sol$width <- elim$width
  sol
}

#' Internal helper: remove dominated vertices
#'
#' A vertex v is dominated by a relative u when every relative of u is also a relative of
#' v, that is when the closed neighbourhood of u (u and its relatives) is contained in
#' that of v. Some maximum independent set then
#' avoids v, so removing v leaves the maximum size unchanged; it does change the set of
#' maximum independent sets, so the rule suits the search for one maximum set but not
#' counting or uniform sampling. All the dominated vertices are removed at once, keeping
#' the one with the smallest index among vertices with identical closed neighbourhoods,
#' so that every removed vertex keeps a dominator; the rule is applied again until no
#' vertex is dominated.
#'
#' @param A Symmetric logical adjacency matrix without self-loops.
#' @return A logical vector, \code{TRUE} for the vertices kept.
#'
#' @keywords internal
#' @noRd
.rkn_dominance_reduce <- function(A) {
  keep <- rep(TRUE, nrow(A))
  repeat {
    idx <- which(keep)
    closed <- A[idx, idx, drop = FALSE]
    diag(closed) <- TRUE
    size <- rowSums(closed)
    drop <- rep(FALSE, length(idx))
    for (v in seq_along(idx)) {
      nb <- setdiff(which(closed[v, ]), v)
      if (length(nb) == 0) next
      outside <- rowSums(closed[nb, , drop = FALSE] & !rep(closed[v, ], each = length(nb)))
      inside <- nb[outside == 0]
      if (any(size[inside] < size[v] | inside < v)) drop[v] <- TRUE
    }
    if (!any(drop)) return(keep)
    keep[idx[drop]] <- FALSE
  }
}

#' Internal helper: draw a maximum independent set of a graph
#'
#' Draws a maximum independent set uniformly at random with the exact engine. When the
#' elimination width exceeds \code{max_width}, the dominated vertices are removed first:
#' the set drawn is still of maximum size but no longer uniform. If the graph is still
#' too wide, the function stops before allocating the tables, rather than exhausting the
#' memory.
#'
#' @param g An igraph object with vertex names.
#' @param max_width Largest elimination width accepted, as in \code{.rkn_mis_exact()}.
#' @return A list with \code{kept}, the names of the selected vertices, and
#'   \code{method}, \code{"exact"} or \code{"exact_reduced"} (after removing dominated
#'   vertices).
#'
#' @keywords internal
#' @noRd
.rkn_draw_maximum_set <- function(g, max_width = 20) {
  ids <- igraph::V(g)$name
  A <- .rkn_adjacency(g)
  keep <- rep(TRUE, length(ids))
  method <- "exact"
  sol <- .rkn_mis_exact(A, max_width)
  if (!sol$feasible) {
    keep <- .rkn_dominance_reduce(A)
    sol <- .rkn_mis_exact(A[keep, keep, drop = FALSE], max_width)
    method <- "exact_reduced"
  }
  if (!sol$feasible) {
    stop(sprintf(paste0("A family with %d individuals left after the structural rules is too intricate for ",
                        "an exact selection (elimination width %d, limit %d). Exclude the samples flagged by ",
                        "rkn_check_samples(), or use a stricter threshold such as \"second\"."),
                 length(ids), sol$width, max_width), call. = FALSE)
  }
  list(kept = ids[keep][.rkn_mis_sample(sol)], method = method)
}

#' Internal helper: logical adjacency matrix of a graph
#'
#' @keywords internal
#' @noRd
.rkn_adjacency <- function(g) as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE)) > 0

#' Internal helper: scores of the two choices of a vertex in many partial assignments
#'
#' For each column of \code{x}, whose vertices later than \code{v} in the elimination
#' order are already decided, the best total weight reachable by the vertices eliminated
#' up to \code{v} and the number of ways to reach it, when \code{v} is left out
#' (\code{s0}, \code{c0}) or selected (\code{s1}, \code{c1}).
#'
#' @keywords internal
#' @noRd
.rkn_mis_choices <- function(sol, v, x) {
  n_col <- ncol(x)
  scope <- sol$later[[v]]
  s0 <- numeric(n_col)
  s1 <- rep(sol$weight[v], n_col)
  c0 <- rep(1, n_col)
  c1 <- rep(1, n_col)
  relatives <- scope[sol$A[v, scope]]
  if (length(relatives) > 0) {
    clash <- colSums(x[relatives, , drop = FALSE]) > 0
    s1[clash] <- -Inf
    c1[clash] <- 0
  }
  for (w in sol$inbox[[v]]) {
    m <- sol$msg[[w]]
    powers <- 2^(seq_along(m$vars) - 1)
    is_v <- m$vars == v
    others <- m$vars[!is_v]
    base <- if (length(others) > 0) colSums(x[others, , drop = FALSE] * powers[!is_v]) else numeric(n_col)
    i0 <- base + 1
    i1 <- base + powers[is_v] + 1
    s0 <- s0 + m$s[i0]
    c0 <- c0 * m$c[i0]
    s1 <- s1 + m$s[i1]
    c1 <- c1 * m$c[i1]
  }
  list(s0 = s0, s1 = s1, c0 = c0, c1 = c1)
}

#' Internal helper: draw many maximum independent sets at once
#'
#' Vectorised counterpart of \code{.rkn_mis_sample()}: the draws are independent and
#' each is uniform over the sets of maximum size (or total weight).
#'
#' @param sol The result of \code{.rkn_mis_eliminate()}.
#' @param n_draws Number of draws.
#' @return A logical matrix with one row per vertex and one column per draw.
#'
#' @keywords internal
#' @noRd
.rkn_mis_sample_many <- function(sol, n_draws) {
  x <- matrix(0L, nrow(sol$A), n_draws)
  for (v in rev(sol$order)) {
    ch <- .rkn_mis_choices(sol, v, x)
    best <- pmax(ch$s0, ch$s1)
    p0 <- (ch$s0 == best) * ch$c0
    p1 <- (ch$s1 == best) * ch$c1
    x[v, ] <- as.integer(p1 > 0 & (p0 == 0 | stats::runif(n_draws) < p1 / (p0 + p1)))
  }
  x == 1
}

#' Internal helper: list all the maximum independent sets
#'
#' Goes through the vertices in reverse elimination order, keeping every choice that
#' reaches the maximum size (or total weight), so that each set is listed once.
#'
#' @param sol The result of \code{.rkn_mis_eliminate()}.
#' @param limit Largest number of sets listed.
#' @return A logical matrix with one row per vertex and one column per set, or
#'   \code{NULL} when there are more than \code{limit} sets.
#'
#' @keywords internal
#' @noRd
.rkn_mis_enumerate <- function(sol, limit) {
  if (sol$log_count > log(limit) + 1e-9) return(NULL)
  x <- matrix(0L, nrow(sol$A), 1)
  for (v in rev(sol$order)) {
    ch <- .rkn_mis_choices(sol, v, x)
    best <- pmax(ch$s0, ch$s1)
    keep0 <- which(ch$s0 == best & ch$c0 > 0)
    keep1 <- which(ch$s1 == best & ch$c1 > 0)
    x1 <- x[, keep1, drop = FALSE]
    x1[v, ] <- 1L
    x <- cbind(x[, keep0, drop = FALSE], x1)
  }
  x == 1
}

#' Internal helper: exact inclusion probabilities
#'
#' For each vertex, the proportion of the maximum independent sets (or the sets of
#' maximum total weight) that contain it, i.e. its probability of being selected by a
#' uniform draw. The elimination of \code{.rkn_mis_eliminate()} passes information from
#' the first vertices eliminated to the last; this function passes it back. Going
#' through the vertices in reverse order, the table of a vertex combines its own
#' constraint, the messages it received and the message from the vertices eliminated
#' after it, so it accounts for the whole graph: its largest value is the maximum size,
#' and the counts that reach it, split by the choice of the vertex, give the
#' probability. The message passed back to each vertex that sent a message removes that
#' message from the table (messages are always finite, since leaving every vertex out is
#' always possible) and keeps the maximum over the vertices outside the scope of the
#' receiver.
#'
#' @param sol The result of \code{.rkn_mis_eliminate()}.
#' @return A numeric vector with one probability per vertex.
#'
#' @keywords internal
#' @noRd
.rkn_mis_inclusion <- function(sol) {
  n <- nrow(sol$A)
  prob <- numeric(n)
  down <- vector("list", n)
  for (p in rev(sol$order)) {
    vars <- c(p, sol$later[[p]])
    code <- 0:(2^length(vars) - 1)
    bit <- function(j) (code %/% 2^j) %% 2
    index_of <- function(m) {
      index <- 0
      for (k in seq_along(m$vars)) index <- index + bit(match(m$vars[k], vars) - 1) * 2^(k - 1)
      index + 1
    }
    xp <- bit(0)
    S <- xp * sol$weight[p]
    C <- rep(1, length(code))
    relatives <- sol$later[[p]][sol$A[p, sol$later[[p]]]]
    if (length(relatives) > 0) {
      clash <- rep(FALSE, length(code))
      for (u in relatives) clash <- clash | bit(match(u, vars) - 1) == 1
      clash <- clash & xp == 1
      S[clash] <- -Inf
      C[clash] <- 0
    }
    received <- vector("list", length(sol$inbox[[p]]))
    for (k in seq_along(sol$inbox[[p]])) {
      m <- sol$msg[[sol$inbox[[p]][k]]]
      received[[k]] <- index_of(m)
      S <- S + m$s[received[[k]]]
      C <- C * m$c[received[[k]]]
    }
    if (!is.null(down[[p]])) {
      i <- index_of(down[[p]])
      S <- S + down[[p]]$s[i]
      C <- C * down[[p]]$c[i]
    }
    top <- S == max(S)
    prob[p] <- sum(C[top & xp == 1]) / sum(C[top])

    for (k in seq_along(sol$inbox[[p]])) {
      w <- sol$inbox[[p]][k]
      s <- S - sol$msg[[w]]$s[received[[k]]]
      cnt <- C / sol$msg[[w]]$c[received[[k]]]
      kept_vars <- vars
      # Drop the variables outside the scope of w, from the highest bit down, so that the
      # positions of the remaining lower bits do not change.
      for (j in rev(seq_along(vars))) {
        if (vars[j] %in% sol$later[[w]]) next
        step <- 2^(j - 1)
        i0 <- which(((seq_along(s) - 1) %/% step) %% 2 == 0)
        i1 <- i0 + step
        best <- pmax(s[i0], s[i1])
        cnt <- (s[i0] == best) * cnt[i0] + (s[i1] == best) * cnt[i1]
        s <- best
        kept_vars <- kept_vars[-j]
      }
      if (max(cnt) > 1e100) cnt <- cnt / max(cnt)
      down[[w]] <- list(vars = kept_vars, s = s, c = cnt)
    }
  }
  prob
}
