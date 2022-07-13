#' @title plot_mapper
#' @description Generates plots of the network color-coded by
#' each of the specified markers.
#' @param mapper Existing mapper object.
#' @param markers A subset of the marker names.
#' @param path The path where plots are saved.
#' @param device The device used for image encoding. Supports
#' png for bitmaps, pdf for vector graphics.
#' @export
plot_mapper <- function(mapper, markers=colnames(mapper$node_stats$q50),
                       path = NULL, device="png", remove_smaller=0) {

  if (device!="png" & device!="pdf")
    stop("Supported devices are png and pdf.")

  size <- sapply(mapper$nodes, length)
  keep <- which(size > remove_smaller)
  size <- size[keep]
  community <- mapper$community[keep]

  gr <- mapper$gr %>%
    delete_vertices(setdiff(seq_along(mapper$nodes), keep))

  layout <- create_layout(gr, layout="fr")

  if (!is.null(mapper$community))
    comm <- data.frame(x=layout$x,
                       y=layout$y,
                       c=community) %>%
    group_by(c) %>%
    summarise(x=median(x), y=median(y))

  for (marker in markers) {
    message(marker)
    g <- ggraph(layout) +
      geom_edge_link(aes(alpha = weight)) +
      geom_node_point(shape=21, aes(size = size,
                                    fill = mapper$node_stats$q50[keep,marker])) +
      scale_fill_gradient2(low = "white", mid="white",
                           high = "red",name = marker) +
      scale_size(range=c(1,6), name="count") +
      scale_edge_alpha(guide="none") +
      theme_graph(base_family = "sans")

    if (!is.null(mapper$community))
      g <- g + geom_label(data=comm, aes(x=x,y=y,label=c),
                          alpha=0.5, inherit.aes = FALSE)

    if (is.null(path)) {
      plot(g)
    } else {
      if (device == "png") {
        png(paste0(path, marker, ".png"), width=1200, height=1000, res=120)
        plot(g)
        dev.off()
      } else {
        pdf(paste0(path, marker, ".pdf"), width=12, height=10)
        plot(g)
        dev.off()
      }
    }
  }

  # if (!is.null(mapper$community)) {
  #   g <- ggraph(layout) +
  #     geom_edge_link(aes(alpha = weight)) +
  #     geom_node_point(aes(color=community, size=size)) +
  #     scale_color_discrete(name="community") +
  #     scale_edge_alpha(guide="none") +
  #     scale_size(range=c(1,6), name="count") +
  #     theme_graph(base_family = "sans") +
  #     theme(text=element_text(size = 18))
  #
  #   if (is.null(path)) {
  #     plot(g)
  #   } else {
  #     if (device == "png") {
  #       png(paste0(path, "Community.png"), width=1200, height=1000)
  #       plot(g)
  #       dev.off()
  #     } else {
  #       pdf(paste0(path, "Community.pdf"), width=12, height=10)
  #       plot(g)
  #       dev.off()
  #     }
  #   }
  # }
}


#' @title plot_mapper_interactive
#' @description Generates an interactive plot of the network and
#' saves it as an html file.
#' @param mapper Existing mapper object.
#' @param color A factor or numeric vector which is mapped to the
#' color of the nodes.
#' @param labels A character vector describing each node, displayed
#' when hovering over the node.
#' @param path The path where plots are saved.
#' @param title The name of the output html file.
#' @import ggiraph
#' @export
plot_mapper_interactive <- function(mapper, color=mapper$community,
                                  labels=paste0("Community: ", mapper$community, "\n"),
                                  path="", title="HiTMapper") {
  size <- sapply(mapper$nodes, length)
  layout <- create_layout(mapper$gr, layout="fr")
  css <- "background-color:gray;color:white;padding:10px;border-radius:5px;"

  g <- ggraph(layout) +
    geom_edge_link(aes(alpha=weight)) +
    geom_point_interactive(mapping = aes(x = x, y = y, size=size,
                                         color = color,
                                         data_id = mapper$community,
                                         tooltip = labels)) +
    scale_size(range = c(0.5,3)) +
    scale_edge_alpha(guide="none") +
    guides(color = "none") +
    theme_graph(base_family = "sans")

  myplot <- girafe(ggobj = g,
                   options = list(
                     opts_tooltip(css=css),
                     opts_hover(css="stroke:black;"),
                     opts_hover_inv(css="opacity:0.3;")
                   ))
  out_name <- paste0(path, title, ".html")
  htmlwidgets::saveWidget(myplot, out_name)
}


