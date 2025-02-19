library(ggplot2)
library(devEMF)

save_plot <- function(filename, plot, plot_width,plot_height, scale_f = scalefactor,dpi_val = 600) {
  # Save as EMF
  ggsave(
    filename = paste0(filename, ".emf"),
    plot = plot,
    dpi = dpi_val,
    units = "in",
    width = plot_width * scale_f,
    height = plot_height * scale_f,
    device = {function(filename, ...) devEMF::emf(file = filename, ...)}
  )
  
  # Save as JPEG
  ggsave(
    filename = paste0(filename, ".jpeg"),
    plot = plot,
    dpi = dpi_val,
    units = "in",
    width = plot_width * scale_f,
    height = plot_height * scale_f,
    device = "jpeg"
  )

  # ggsave(
  #   filename = paste0(filename, ".png"),
  #   plot = plot,
  #   dpi = dpi_val,
  #   units = "in",
  #   width = plot_width * scale_f,
  #   height = plot_height * scale_f,
  #   device = "png",
  #   bg = "transparent"
  # )
}

# Example usage
#Figure5 <- ggplot(...)  # Replace ... with your actual ggplot code
#save_plot_as_jpeg_and_emf("Figure5_v3_dish5_plate01_region02_PPRphysiology_from_tracking_O-phys_paper", Figure5)
