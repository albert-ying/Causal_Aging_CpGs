require(tidyverse)
require(tidymodels)
require(magrittr)

# new_data = drop_list[[2]]
# coeff_tbl = step.1.feature.list[[2]]
predict_by_coeff = function(new_data, coeff_tbl) {
  intercept = coeff_tbl %>%
    filter(term == "(Intercept)") %>%
    pull(estimate)
  features = coeff_tbl %>%
    filter(term != "(Intercept)") %>%
    pull(term)
  shared_feature = colnames(new_data) %>%
    intersect(features) %>%
    sort()
  estimate = coeff_tbl %>%
    filter(term != "(Intercept)") %>%
    filter(term %in% shared_feature) %>%
    arrange(term) %>%
    pull(estimate) %>%
    as.matrix(1)
  new_mat = new_data[,all_of(shared_feature)] %>%
    as.matrix
  pred = new_mat %*% estimate + intercept
  result = tibble(.pred = as.numeric(pred))
  return(result)
}

## First elastic net model: set function
## function for plot prediction plot with metrics

fit2plot = function(predi, titles) {
  mae = rmse = rsq = pr = pp = rho = NA
  tryCatch({
    mae = predi %>%
      mae(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    rmse = predi %>%
      yardstick::rmse(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    mse = rmse**2 |> round(2)
    rsq = predi %>%
      rsq(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    pearson = predi %$%
      cor.test(age, .pred, method = "pearson")
    spearman = predi %$%
      cor.test(age, .pred, method = "spearman")
    rho = spearman$estimate %>% round(2)
    pr = pearson$estimate %>% round(2)
    pp = pearson$p.value %>% signif()
  }, error = function(e){})
  plot = ggplot(predi, aes(x = age, y = .pred)) |>
    geom_point(color = "black", fill = NA, alpha = 1, size = 2.1, show.legend = F) |>
    base_mode() |>
    geom_smooth(method = "lm", alpha = 0.1) |>
    geom_point(color = "white", alpha = 1, size = 1.6, show.legend = F, pch = 19) +
    geom_point(aes(color = abs(age - .pred)), alpha = 0.2, size = 1.6, show.legend = F, pch = 19) +
    labs(
      title = titles,
      subtitle = glue::glue("MAE = {mae}, R^2 = {rsq}<br><sub>MSE = {mse}, Spearman's Rho = {rho}</sub>"),
      x = "Age",
      y = "Predicted Age"
    ) +
    theme(legend.position = "none")
  return(plot)
}

fit2tbl = function(predi, titles) {
  mae = rmse = rsq = pr = pp = NA
  tryCatch({
    mae = predi %>%
      mae(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    rmse = predi %>%
      yardstick::rmse(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    mse = rmse**2 |> round(2)
    rsq = predi %>%
      rsq(truth = age, .pred) %>%
      pull(.estimate) %>%
      round(2)
    pearson = predi %$%
      cor.test(age, .pred, method = "pearson")
    pr = pearson$estimate %>% round(2)
    pp = pearson$p.value %>% signif()
    spearman = predi %$%
      cor.test(age, .pred, method = "spearman")
    rho = spearman$estimate %>% round(2)
  }, error = function(e){})
  tbl = tibble(
    mae = mae,
    rmse = rmse,
    rsq = rsq,
    pearson_r = pr,
    pearson_p = pp,
    method = titles,
    rho = rho,
    mse = mse
  )
  return(tbl)
}


remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}
