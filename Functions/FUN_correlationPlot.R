correlate_plot = function(cor_output, order = "absolute", ...){
  require(corrplot)
  require(tidyverse)

  full_tb = cor_output |>
    mutate(p = p.adjust(p, method = "bonferroni")) |>
    bind_rows(select(cor_output, Parameter1 = Parameter2, Parameter2 = Parameter1, everything())) |>
    select(Parameter1, Parameter2, rho, p) |>
    mutate(Parameter1 = str_replace_all(Parameter1, "(-|_)", " ")) |>
    mutate(Parameter2 = str_replace_all(Parameter2, "(-|_)", " "))

  rho_matrix = full_tb |>
    select(Parameter1, Parameter2, rho) |>
    pivot_wider(names_from = Parameter2, values_from = rho) |>
    column_to_rownames("Parameter1") |>
    as.matrix()

  rho_matrix[is.na(rho_matrix)] = 1
  rho_matrix[rho_matrix > 1] = 1
  rho_matrix[rho_matrix < -1] = -1
  rho_matrix = rho_matrix[sort(rownames(rho_matrix)), sort(colnames(rho_matrix))]

  p_matrix = full_tb |>
    select(Parameter1, Parameter2, p) |>
    pivot_wider(names_from = Parameter2, values_from = p) |>
    column_to_rownames("Parameter1") |>
    as.matrix()
  p_matrix[is.na(p_matrix)] = 1

  p_matrix = p_matrix[sort(rownames(p_matrix)),sort(colnames(p_matrix))]

  if (order == "absolute") {
    matrix = abs(rho_matrix)
  } else {
    matrix = rho_matrix
  }
  
  ord = corrMatOrder(matrix, order="hclust")
  rho_matrix = rho_matrix[ord, ord]
  p_matrix = p_matrix[ord, ord]

  pdf(...)
  corrplot(
    rho_matrix,
    # order = "hclust",
    order = "original",
    method = "square",
    outline = F,
    addgrid.col = "darkgray",
    type = "full",
    tl.col = "black",
    tl.pos = "tp",
    diag = T,
    tl.cex = 2,
    cl.cex = 2,
    p.mat = p_matrix,
    insig = "label_sig",
    sig.level = c(0.0005, 0.005, .05),
    pch.cex = 1.6,
    pch.col = "black",
    na.label = "square",
    na.label.col = "white",
    col = colorRampPalette(c("#ec4c2f", "white", "#4891a1"))(200)
  )
  dev.off()
}
