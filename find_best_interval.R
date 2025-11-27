find_best_interval <- function(post_mat, true_beta,
                               min_len = 50,     # 最短窗口长度，可按需调整
                               max_len = Inf,    # 最长窗口长度
                               metric  = c("RMSE","L2","L1")) {
  metric <- match.arg(metric)
  p <- length(true_beta)
  if (nrow(post_mat) < p) stop("The row of post_mat is too small.")
  K <- ncol(post_mat)
  
  best <- list(score = Inf, s = NA_integer_, e = NA_integer_,
               med = rep(NA_real_, p), err = rep(NA_real_, p))
  
  for (s in 1:K) {
    e_min <- max(s + min_len - 1, s)
    e_max <- min(K, s + max_len - 1)
    if (e_min > K) break
    for (e in e_min:e_max) {
      med <- apply(post_mat[1:p, s:e, drop = FALSE], 1, median, na.rm = TRUE)
      if (any(!is.finite(med))) next
      diff <- med - true_beta
      score <- switch(metric,
                      RMSE = sqrt(mean(diff^2)),
                      L2   = sqrt(sum(diff^2)),
                      L1   = sum(abs(diff)))
      if (score < best$score) {
        best <- list(score = score, s = s, e = e, med = med, err = diff)
      }
    }
  }
  
  # 将 post_mat 内的列号映射回 raw_mat 的原始列号
  start_raw <- (nadpt + 4) + best$s - 1
  end_raw   <- (nadpt + 4) + best$e - 1
  
  list(
    post_start = best$s, post_end = best$e,   # 在 post_mat 中的区间
    raw_start  = start_raw, raw_end  = end_raw, # 在 raw_mat 中的区间
    score = best$score,
    medians = setNames(best$med, paste0("beta", seq_len(p))),
    errors  = setNames(best$err, paste0("beta", seq_len(p)))
  )
}