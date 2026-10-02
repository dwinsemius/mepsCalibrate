extract_cross_partials = function(gam_model, interaction_name, grid_size = 50) {
  smooth_objs = gam_model$smooth
  ti_idx = which(sapply(smooth_objs, function(x) x$label == interaction_name))
  if (length(ti_idx) == 0) stop("Specified interaction term not found in model.")
  
  v_names = smooth_objs[[ti_idx]]$term
  u_seq = seq(0.05, 0.95, length.out = grid_size)
  v_seq = seq(0.05, 0.95, length.out = grid_size)
  eval_grid = expand.grid(u = u_seq, v = v_seq)
  names(eval_grid) = v_names
  
  h_u = 1e-4
  h_v = 1e-4
  
  grid_pp = eval_grid; grid_pp[[v_names[1]]] = grid_pp[[v_names[1]]] + h_u; grid_pp[[v_names[2]]] = grid_pp[[v_names[2]]] + h_v
  grid_pm = eval_grid; grid_pm[[v_names[1]]] = grid_pm[[v_names[1]]] + h_u; grid_pm[[v_names[2]]] = grid_pm[[v_names[2]]] - h_v
  grid_mp = eval_grid; grid_mp[[v_names[1]]] = grid_mp[[v_names[1]]] - h_u; grid_mp[[v_names[2]]] = grid_mp[[v_names[2]]] + h_v
  grid_mm = eval_grid; grid_mm[[v_names[1]]] = grid_mm[[v_names[1]]] - h_u; grid_mm[[v_names[2]]] = grid_mm[[v_names[2]]] - h_v
  
  all_vars = names(gam_model$model)
  extra_vars = setdiff(all_vars, c(v_names, as.character(formula(gam_model))))
  for (ev in extra_vars) {
    if(ev %in% names(eval_grid)) next
    val = if(is.numeric(gam_model$model[[ev]])) median(gam_model$model[[ev]], na.rm=TRUE) else gam_model$model[[ev]]
    grid_pp[[ev]] = grid_pm[[ev]] = grid_mp[[ev]] = grid_mm[[ev]] = val
  }
  
  X_pp = predict(gam_model, newdata = grid_pp, type = "lpmatrix")
  X_pm = predict(gam_model, newdata = grid_pm, type = "lpmatrix")
  X_mp = predict(gam_model, newdata = grid_mp, type = "lpmatrix")
  X_mm = predict(gam_model, newdata = grid_mm, type = "lpmatrix")
  
  first_param = smooth_objs[[ti_idx]]$first.para
  last_param = smooth_objs[[ti_idx]]$last.para
  col_indices = first_param:last_param
  
  X_diff = (X_pp[, col_indices] - X_pm[, col_indices] - X_mp[, col_indices] + X_mm[, col_indices]) / (4 * h_u * h_v)
  
  beta = coef(gam_model)[col_indices]
  V_beta = vcov(gam_model)[col_indices, col_indices]
  
  deriv_values = as.vector(X_diff %*% beta)
  deriv_se = sqrt(rowSums((X_diff %*% V_beta) * X_diff))
  
  output = data.frame(
    u_coord = eval_grid[[v_names[1]]],
    v_coord = eval_grid[[v_names[2]]],
    derivative = deriv_values,
    se = deriv_se,
    t_stat = deriv_values / deriv_se
  )
  return(output)
}
