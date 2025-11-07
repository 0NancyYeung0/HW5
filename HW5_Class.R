## HW5 Class/Methods

setClass(
  Class = "sparse_numeric",
  slots = c(
    value = "numeric",
    pos = "integer",
    length = "integer"
  )
)


## Validity method
setValidity("sparse_numeric", function(object) {
  
  # Single non-negative int
  if (length(object@length) != 1 || object@length < 0) {
    return("'length' must be a single non-negative integer.")
  }
  
  # Same length
  if (length(object@value) != length(object@pos)) {
    return("'value' and 'pos' must have the same length")
  }
  
  # pos = integer vector
  if (!is.integer(object@pos)) {
    return("'pos' must be an integer vector")
  } else {
    # Check bounds and uniqueness/ordering
    if (length(object@pos) > 0L) {
      if (any(is.na(object@pos))) {
        return("'pos' contains NA")
      }
      if (any(object@pos < 1L)) {
        return("all 'pos' entries must be >= 1")
      }
      if (any(object@pos > as.integer(object@length))) {
        return("'pos' entries cannot exceed 'length'")
      }
      if (any(duplicated(object@pos))) {
        return("'pos' must not contain duplicates")
      }
    }
  }
  
  # Not containing zeros
  if (length(object@value) > 0L) {
    if (any(is.na(object@value))) {
      return("'value' contains NA")
    }
    if (any(object@value == 0)) {
      return("'value' must not contain zeros; zero entries should not be stored")
    }
  }
})


## Coercion methods
setAs("numeric", "sparse_numeric", function(from) {
  nonzero_idx <- which(from != 0)
  new("sparse_numeric",
      value = if (length(nonzero_idx) > 0) as.numeric(from[nonzero_idx]) else numeric(0),
      pos = if (length(nonzero_idx) > 0) as.integer(nonzero_idx) else integer(0),
      length = as.integer(length(from))
  )
})

setAs("sparse_numeric", "numeric", function(from) {
  x <- as.integer(from@length)
  out <- numeric(x)
  if (length(from@pos) > 0L) {
    out[from@pos] <- from@value
  }
  out
})


## Generics
setGeneric("sparse_add", function(x, y, ...) standardGeneric("sparse_add"))
setGeneric("sparse_mult", function(x, y, ...) standardGeneric("sparse_mult"))
setGeneric("sparse_sub", function(x, y, ...) standardGeneric("sparse_sub"))
setGeneric("sparse_crossprod", function(x, y, ...) standardGeneric("sparse_crossprod"))
setGeneric("sparse_norm", function(x, ...) standardGeneric("sparse_norm"))


## Helper: sanitize result
.sanitize_sparse <- function(vals, poss, len) {
  if (length(vals) == 0L) {
    return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = as.integer(len)))
  }
  keep <- !is.na(vals) & vals != 0
  if (!all(keep)) {
    vals <- vals[keep]
    poss <- poss[keep]
  }
  o <- order(poss)
  vals <- vals[o]
  poss <- as.integer(poss[o])
  new("sparse_numeric", value = as.numeric(vals), pos = poss, length = as.integer(len))
}

## Methods for arithmetic

setMethod("sparse_add", c("sparse_numeric", "sparse_numeric"),
          function(x, y, ...) {
            if (x@length != y@length) stop("vectors must have same length")
            # if one is all-zero
            if (length(x@pos) == 0L) return(y)
            if (length(y@pos) == 0L) return(x)
            
            # merge positions
            allpos <- sort(unique(c(x@pos, y@pos)))
            # initialize values to zero
            vals_x <- numeric(length(allpos))
            vals_y <- numeric(length(allpos))
            # map values
            match_x <- match(allpos, x@pos, nomatch = 0L)
            match_y <- match(allpos, y@pos, nomatch = 0L)
            vals_x[match_x != 0L] <- x@value[match_x[match_x != 0L]]
            vals_y[match_y != 0L] <- y@value[match_y[match_y != 0L]]
            res_vals <- vals_x + vals_y
            .sanitize_sparse(res_vals, as.integer(allpos), x@length)
          })

setMethod("sparse_sub", c("sparse_numeric", "sparse_numeric"),
          function(x, y, ...) {
            if (x@length != y@length) stop("vectors must have same length")
            if (length(x@pos) == 0L && length(y@pos) == 0L) {
              return(new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length))
            }
            if (length(x@pos) == 0L) {
              # -y
              neg_vals <- -y@value
              return(new("sparse_numeric", value = neg_vals, pos = y@pos, length = x@length))
            }
            if (length(y@pos) == 0L) return(x)
            
            allpos <- sort(unique(c(x@pos, y@pos)))
            vals_x <- numeric(length(allpos))
            vals_y <- numeric(length(allpos))
            match_x <- match(allpos, x@pos, nomatch = 0L)
            match_y <- match(allpos, y@pos, nomatch = 0L)
            vals_x[match_x != 0L] <- x@value[match_x[match_x != 0L]]
            vals_y[match_y != 0L] <- y@value[match_y != 0L]
            res_vals <- vals_x - vals_y
            .sanitize_sparse(res_vals, as.integer(allpos), x@length)
          })

setMethod("sparse_mult", c("sparse_numeric", "sparse_numeric"),
          function(x, y, ...) {
            # Multiplication only matters for overlapping indices
            common <- intersect(x@pos, y@pos)
            vals <- x@value[match(common, x@pos)] * y@value[match(common, y@pos)]
            new("sparse_numeric", value = numeric(0), pos = integer(0), length = x@length)
            .sanitize_sparse(vals, as.integer(common), x@length)
          })

setMethod("sparse_crossprod", c("sparse_numeric", "sparse_numeric"),
          function(x, y, ...) {
            # dot product
            common <- intersect(x@pos, y@pos)
            sum(x@value[match(common, x@pos)] * y@value[match(common, y@pos)])
          })

setMethod("sparse_norm", "sparse_numeric", function(x, ...) {
  if (length(x@value) == 0L) return(0)
  sqrt(sum(x@value^2))
})


## Link operators to generics
setMethod("+", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_add(e1, e2))
setMethod("-", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_sub(e1, e2))
setMethod("*", c("sparse_numeric", "sparse_numeric"),
          function(e1, e2) sparse_mult(e1, e2))


## show method
setMethod("show", "sparse_numeric", function(object) {
  nnz <- length(object@pos)
  cat("An object of class 'sparse_numeric'\n")
  cat("  length:", as.integer(object@length), "\n")
  cat("  non-zero entries:", nnz, "\n")
  if (nnz > 0L) {
    show_n <- min(10L, nnz)
    df <- data.frame(pos = object@pos[seq_len(show_n)],
                     value = object@value[seq_len(show_n)])
    print(df, row.names = FALSE)
    if (nnz > show_n) cat("  ... (", nnz - show_n, "more non-zero entries)\n", sep = "")
  } else {
    cat("  (all zeros)\n")
  }
  invisible(NULL)
})


## plot method for two sparse_numeric objects
setMethod("plot", signature(x = "sparse_numeric", y = "sparse_numeric"),
          function(x, y, xlab = "Position", ylab = "Value", main = "Non-zero overlap", ...) {
            if (x@length != y@length) stop("vectors must have same length")
            allpos <- sort(unique(c(x@pos, y@pos)))
            rngx <- range(c(x@value, y@value), finite = TRUE)
            if (length(rngx) == 0L || any(is.infinite(rngx))) rngx <- c(-1, 1)
            plot(NA, xlim = range(allpos, na.rm = TRUE), ylim = rngx,
                 xlab = xlab, ylab = ylab, main = main, ...)
            if (length(x@pos) > 0L) points(x@pos, x@value, pch = 17)
            if (length(y@pos) > 0L) points(y@pos, y@value, pch = 19)
            common <- intersect(x@pos, y@pos)
            if (length(common) > 0L) {
              ix <- match(common, x@pos)
              iy <- match(common, y@pos)
              for (k in seq_along(common)) {
                segments(x0 = common[k], y0 = x@value[ix[k]],
                         x1 = common[k], y1 = y@value[iy[k]])
              }
            }
            legend("topright", legend = c("x non-zero", "y non-zero"),
                   pch = c(17, 19), bty = "n")
            invisible(NULL)
          })
