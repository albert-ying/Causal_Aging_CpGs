pcorrmap = function(df, ...) {
  require(pheatmap)
  require(RColorBrewer)
  edge = 1
  npal = 90
  myBreaks = c(
    seq(-edge, 0, length.out = ceiling(npal / 2) + 1),
    seq(edge / npal, edge, length.out = floor(npal / 2))
  )
  pheatmap(df,
    breaks = myBreaks,
    cellwidth = 10,
    cellheight = 10,
    border_color = NA,
    fontsize_number = 3,
    fontsize = 10,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(90),
    angle_col = 90,
    ...
  )
}
