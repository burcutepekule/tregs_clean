steady_state_idx = function(x, k = 50, tail_frac = 0.25, min_run = 20) {
  
  n = length(x)
  tail_n = ceiling(n * tail_frac)
  
  tail_mean = mean(tail(x, tail_n))
  tail_sd = sd(tail(x, tail_n))
  
  signal_range = diff(range(x))
  
  if (signal_range < 1e-10) return(1L)
  
  tol_mean = max(0.5 * tail_sd, 0.05 * abs(tail_mean), 0.01 * signal_range)
  tol_sd = max(0.5 * tail_sd, 0.05 * abs(tail_mean), 0.01 * signal_range)
  
  m = zoo::rollapply(x, k, mean, align = "right", fill = NA)
  s = zoo::rollapply(x, k, sd, align = "right", fill = NA)
  sl = zoo::rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  
  mean_ok = abs(m - tail_mean) < tol_mean
  sd_ok = abs(s - tail_sd) < tol_sd
  
  slope_tol = 0.01 * signal_range / k
  slope_ok = abs(sl) < slope_tol
  
  close_enough = abs(m - tail_mean) < 0.1 * signal_range
  
  # Relaxed SD check: within 2x of tail_sd OR sd is small relative to signal
  sd_relaxed = s < max(2 * tail_sd, 0.05 * signal_range)
  
  ok = (mean_ok & sd_ok) | (slope_ok & sd_relaxed & close_enough)
  ok[is.na(ok)] = FALSE
  
  runs = rle(ok)
  run_ends = cumsum(runs$lengths)
  run_starts = c(1, head(run_ends, -1) + 1)
  
  good_runs = which(runs$values & runs$lengths >= min_run)
  
  if (length(good_runs) == 0) return(NA_integer_)
  
  run_starts[good_runs[1]]
}

steady_state_idx_slow = function(x, k = 50, tail_frac = 0.25, min_run = 20) {
  
  n = length(x)
  tail_n = ceiling(n * tail_frac)
  
  tail_mean = mean(tail(x, tail_n))
  tail_sd = sd(tail(x, tail_n))
  
  # Signal range for scaling tolerances
  signal_range = diff(range(x))
  
  tol_mean = max(0.5 * tail_sd, 0.05 * abs(tail_mean))
  tol_sd = max(0.5 * tail_sd, 0.05 * abs(tail_mean))
  
  m = zoo::rollapply(x, k, mean, align = "right", fill = NA)
  s = zoo::rollapply(x, k, sd, align = "right", fill = NA)
  sl = zoo::rollapply(x, k, function(v) coef(lm(v ~ seq_along(v)))[2], align = "right", fill = NA)
  
  mean_ok = abs(m - tail_mean) < tol_mean
  sd_ok = abs(s - tail_sd) < tol_sd
  
  # Slope tolerance - signal has stopped changing meaningfully
  slope_tol = 0.01 * signal_range / k
  slope_ok = abs(sl) < slope_tol
  
  # Relaxed mean criterion - within 10% of signal range from tail mean
  close_enough = abs(m - tail_mean) < 0.1 * signal_range
  
  ok = (mean_ok & sd_ok) | (slope_ok & sd_ok & close_enough)
  ok[is.na(ok)] = FALSE
  
  runs = rle(ok)
  run_ends = cumsum(runs$lengths)
  run_starts = c(1, head(run_ends, -1) + 1)
  
  good_runs = which(runs$values & runs$lengths >= min_run)
  
  if (length(good_runs) == 0) return(NA_integer_)
  
  run_starts[good_runs[1]]
}

compute_oscillation_metrics = function(x) {
  x = na.omit(x)
  if (length(x) < 20) {
    return(list(
      acf_peak_lag = NA,
      acf_peak_value = NA,
      acf_regularity = NA,
      spectral_concentration = NA,
      spectral_sharpness = NA
    ))
  }
  
  # === ACF metrics ===
  acf_result = acf(x, lag.max = floor(length(x)/2), plot = FALSE)
  acf_vals = as.numeric(acf_result$acf[-1])
  
  # Find all positive peaks in ACF
  peaks_acf = which(diff(sign(diff(acf_vals))) == -2) + 1
  positive_peaks = peaks_acf[acf_vals[peaks_acf] > 0]
  
  if (length(positive_peaks) > 0) {
    first_peak_idx = positive_peaks[1]
    acf_peak_lag = first_peak_idx
    acf_peak_value = acf_vals[first_peak_idx]
  } else {
    acf_peak_lag = NA
    acf_peak_value = 0
  }
  
  # === ACF regularity ===
  if (length(positive_peaks) >= 3) {
    acf_intervals = diff(positive_peaks)
    acf_regularity = sd(acf_intervals) / mean(acf_intervals)
  } else {
    acf_regularity = NA
  }
  
  # === Spectral metrics ===
  spec_result = spectrum(x, plot = FALSE, detrend = TRUE, taper = 0.1)
  power = spec_result$spec
  
  spectral_concentration = max(power) / sum(power)
  
  max_idx = which.max(power)
  window_start = max(1, max_idx - 2)
  window_end = min(length(power), max_idx + 2)
  neighborhood = power[window_start:window_end]
  neighborhood_mean = mean(neighborhood[neighborhood != max(power)])
  if (is.na(neighborhood_mean) || neighborhood_mean == 0) {
    spectral_sharpness = NA
  } else {
    spectral_sharpness = max(power) / neighborhood_mean
  }
  
  return(list(
    acf_peak_lag = acf_peak_lag,
    acf_peak_value = acf_peak_value,
    acf_regularity = acf_regularity,
    spectral_concentration = spectral_concentration,
    spectral_sharpness = spectral_sharpness
  ))
}

# Function to calculate Cohen's d
cohens_d = function(x, y) {
  nx = length(x)
  ny = length(y)
  mx = mean(x, na.rm = TRUE)
  my = mean(y, na.rm = TRUE)
  sx = sd(x, na.rm = TRUE)
  sy = sd(y, na.rm = TRUE)
  
  # Pooled standard deviation
  pooled_sd = sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  
  # Cohen's d
  d = (my - mx) / pooled_sd
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(0)  # or NA, depending on how you want to interpret it
  }
  return(d)
}

calculate_js_divergence_simple = function(vec1, vec2, n_bins = 50) {
  # Create matching histograms with the same breaks
  breaks = seq(min(c(vec1, vec2)),
                max(c(vec1, vec2)),
                length.out = n_bins + 1)
  
  hist1 = hist(vec1, breaks = breaks, plot = FALSE)
  hist2 = hist(vec2, breaks = breaks, plot = FALSE)
  
  # Convert counts to probabilities
  p = hist1$counts / sum(hist1$counts)
  q = hist2$counts / sum(hist2$counts)
  
  # Calculate JS divergence
  js_div = JSD(rbind(p, q), unit = "log2")
  
  return(js_div)
}


calculate_js_divergence = function(vec1, vec2, n_bins = NULL, method = "sturges") {

  # Check for empty or invalid input vectors
  if (length(vec1) == 0 || length(vec2) == 0) {
    warning("Empty input vector detected, returning JS divergence of NA")
    return(c(NA, 0))
  }

  # Remove NA values
  vec1 = vec1[!is.na(vec1)]
  vec2 = vec2[!is.na(vec2)]

  # Check again after removing NAs
  if (length(vec1) == 0 || length(vec2) == 0) {
    warning("No valid data after removing NAs, returning JS divergence of NA")
    return(c(NA, 0))
  }

  # If n_bins not specified, calculate it
  if (is.null(n_bins)) {
    n1 = length(vec1)
    n2 = length(vec2)
    n = n1 + n2  # Total sample size
    
    n_bins = switch(method,
                     # Sturges' rule (default in R's hist)
                     "sturges" = ceiling(log2(n) + 1),
                     
                     # Scott's rule (assumes normal-like distribution)
                     "scott" = {
                       sd_val = sd(c(vec1, vec2))
                       if (sd_val == 0 || is.na(sd_val)) {
                         ceiling(log2(n) + 1)  # Fallback to Sturges
                       } else {
                         h = 3.5 * sd_val / (n^(1/3))
                         bins = ceiling((max(c(vec1, vec2)) - min(c(vec1, vec2))) / h)
                         if (is.na(bins) || bins <= 0 || is.infinite(bins)) {
                           ceiling(log2(n) + 1)  # Fallback
                         } else {
                           bins
                         }
                       }
                     },
                     
                     # Freedman-Diaconis rule (robust to outliers)
                     "fd" = {
                       iqr_val = IQR(c(vec1, vec2))
                       if (iqr_val == 0 || is.na(iqr_val)) {
                         ceiling(log2(n) + 1)  # Fallback to Sturges if IQR is zero - this happens when distributions are very very different
                       } else {
                         h = 2 * iqr_val / (n^(1/3))
                         bins = ceiling((max(c(vec1, vec2)) - min(c(vec1, vec2))) / h)
                         if (is.na(bins) || bins <= 0 || is.infinite(bins)) {
                           ceiling(log2(n) + 1)  # Fallback
                         } else {
                           bins
                         }
                       }
                     },
                     
                     # Square root rule
                     "sqrt" = ceiling(sqrt(n)),
                     
                     # Rice rule
                     "rice" = ceiling(2 * n^(1/3)),
                     
                     # Default
                     ceiling(log2(n) + 1)
    )
    
    # Ensure reasonable bounds
    n_bins = max(10, min(n_bins, 200))
  }

  # Check if the data has zero or near-zero range
  data_min = min(c(vec1, vec2))
  data_max = max(c(vec1, vec2))
  data_range = data_max - data_min

  # If range is zero or extremely small, the distributions are essentially identical
  # Return JS divergence of 0
  if (data_range == 0 || is.na(data_range) || data_range < .Machine$double.eps * 100) {
    return(c(0, n_bins))
  }

  # Create matching histograms with the same breaks
  breaks = seq(data_min, data_max, length.out = n_bins + 1)

  hist1 = hist(vec1, breaks = breaks, plot = FALSE)
  hist2 = hist(vec2, breaks = breaks, plot = FALSE)

  # Check if histograms are valid (non-empty counts)
  if (sum(hist1$counts) == 0 || sum(hist2$counts) == 0) {
    warning("Empty histogram detected, returning JS divergence of 0")
    return(c(0, n_bins))
  }

  # Convert counts to probabilities
  p = hist1$counts / sum(hist1$counts)
  q = hist2$counts / sum(hist2$counts)

  # Check for NaN or invalid probabilities
  if (any(is.na(p)) || any(is.na(q)) || any(is.nan(p)) || any(is.nan(q))) {
    warning("Invalid probability distribution detected, returning JS divergence of 0")
    return(c(0, n_bins))
  }

  # Calculate JS divergence
  js_div = suppressMessages(JSD(rbind(p, q), unit = "log2"))

  # Check if JS divergence is valid
  if (is.na(js_div) || is.nan(js_div)) {
    warning("Invalid JS divergence result, returning 0")
    return(c(0, n_bins))
  }

  return(c(js_div, n_bins))
}

### note that the tolerance values are adjusted for epithelial score range 
steady_state_idx_old = function(x, k = 20, tail_frac = 0.25,
                            tol_abs = 0.05*(150-25), # 2.5 units
                            tol_sd  = 0.02*(150-25),# 1.25 units
                            tol_slope = 0.005*(150-25))  # 0.125/step
{
  n      = length(x)
  tail_n = ceiling(n*tail_frac)
  x_asym = mean(tail(x, tail_n))
  
  m    = rollapply(x, k, mean, align = "right", fill = NA)
  sd   = rollapply(x, k, sd,   align = "right", fill = NA)
  sl   = rollapply(x, k, function(v) mean(diff(v)), align = "right", fill = NA)
  cand = which(abs(m-x_asym) <= tol_abs & sd <= tol_sd & abs(sl) <= tol_slope)
  
  if (length(cand) == 0) return(NA_integer_)
  cand[1]
}