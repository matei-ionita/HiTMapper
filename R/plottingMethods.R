#' @title plotMapper
#' @description Generates plots of the network color-coded by
#' each of the specified markers.
#' @param mapper Existing mapper object.
#' @param markers A subset of the marker names.
#' @param community A factor giving community membership of the
#' nodes. If present, an additional plot color-coded by community
#' membership is generated.
#' @param path The path where plots are saved.
#' @param device The device used for image encoding. Supports
#' png for bitmaps, pdf for vector graphics.
#' @export
plotMapper <- function(mapper, markers=colnames(mapper$nodeStats$q50),
                       community=NULL, path = "", device="png") {

  if (device!="png" & device!="pdf")
    stop("Supported devices are png and pdf.")

  layout <- create_layout(mapper$gr, layout="fr")
  size <- sapply(mapper$nodes, length)

  for (marker in markers) {
    message(marker)
    g <- ggraph(layout) +
      geom_edge_link(aes(alpha = weight)) +
      geom_node_point(aes(size=size, color=mapper$nodeStats$q50[,marker])) +
      scale_color_gradient(low = "black", high = "red", name = marker) +
      scale_size(range=c(2,6)) +
      scale_edge_alpha(guide="none") +
      theme_graph(base_family = "sans")

    if (device == "png") {
      png(paste0(path, marker, ".png"), width=1200, height=1000)
      plot(g)
      dev.off()
    } else {
      pdf(paste0(path, marker, ".pdf"), width=12, height=10)
      plot(g)
      dev.off()
    }
  }

  if (!is.null(community)) {
    g <- ggraph(layout) +
      geom_edge_link(aes(alpha = weight)) +
      geom_node_point(aes(color=community, size=size)) +
      scale_edge_alpha(guide="none") +
      theme_graph(base_family = "sans") +
      theme(text=element_text(size = 18))

    if (device == "png") {
      png(paste0(path, "Community.png"), width=1200, height=1000)
      plot(g)
      dev.off()
    } else {
      pdf(paste0(path, "Community.pdf"), width=12, height=10)
      plot(g)
      dev.off()
    }
  }
}


#' @title plotMapperInteractive
#' @description Generates an interactive plot of the network and
#' saves it as an html file.
#' @param mapper Existing mapper object.
#' @param community A factor giving community membership of the
#' nodes.
#' @param color A factor or numeric vector which is mapped to the
#' color of the nodes.
#' @param labels A character vector describing each node, displayed
#' when hovering over the node.
#' @param path The path where plots are saved.
#' @param title The name of the output html file.
#' @import ggiraph
#' @export
plotMapperInteractive <- function(mapper, community, color=community,
                                  labels=paste0("Community: ", community, "\n"),
                                  path="", title="HiTMapper") {
  size <- sapply(mapper$nodes, length)
  layout <- create_layout(mapper$gr, layout="fr")
  css <- "background-color:gray;color:white;padding:10px;border-radius:5px;"

  g <- ggraph(layout) +
    geom_edge_link(aes(alpha=weight)) +
    geom_point_interactive(mapping = aes(x = x, y = y, size=size,
                                         color = color,
                                         data_id = community,
                                         tooltip = labels)) +
    scale_size(range = c(0.5,3)) +
    scale_edge_alpha(guide="none") +
    guides(color = FALSE) +
    theme_graph(base_family = "sans")

  myplot <- girafe(ggobj = g,
                   options = list(
                     opts_tooltip(css=css),
                     opts_hover(css="stroke:black;"),
                     opts_hover_inv(css="opacity:0.3;")
                   ))
  outName <- paste0(path, title, ".html")
  htmlwidgets::saveWidget(myplot, outName)
}




