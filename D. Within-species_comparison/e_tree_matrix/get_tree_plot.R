suppressPackageStartupMessages({
  library(treeio)
  library(ggtree)
  library(ggplot2)
})

# ==========================================
# 1. INPUT: Path to your Newick file
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript get_tree_plot.R <newick_file1> [<newick_file2> ...]", call. = FALSE)
}

# Process each provided Newick file and produce one plot per file
for (newick_file in args) {
  if (!nzchar(newick_file)) next

  if (!file.exists(newick_file)) {
    warning(paste("File not found, skipping:", newick_file))
    next
  }

# ==========================================
# 2. LOAD & PLOT THE TREE
# ==========================================
# Read the Newick file
  tree <- read.newick(newick_file)

# Create a beautiful left-to-right tree
beautiful_tree <- ggtree(tree, layout = "rectangular", linewidth = 0.5, color = "#2c3e50") +
  
  # Expand limits so long text labels do not get cut off at the edges
  hexpand(0.3, direction = 1) + 
  
  # Add an elegant, minimalist title
    labs(title = "Phylogenetic Tree Visualization",
       subtitle = paste("Generated from:", basename(newick_file))) +

  # Add and style tip labels (organism names)
  geom_tiplab(size = 3.5, color = "#34495e", fontface = "italic", offset = 0.01) +
  
  # Clean up theme elements
  theme_tree2() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "#2c3e50"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7f8c8d"),
    legend.position = "none"
  )

# Display the tree in your RStudio plot viewer
  # Display the tree (useful in interactive sessions)
  print(beautiful_tree)

  # Save to high-quality PNG file named after the input Newick
  out_fname_png <- paste0(tools::file_path_sans_ext(basename(newick_file)), "_tree.png")
  ggsave(
    filename = out_fname_png,
    plot = beautiful_tree,
    width = 10,
    height = 10,
    dpi = 300
  )
  
  # Save to high-quality PDF file named after the input Newick
  out_fname_pdf <- paste0(tools::file_path_sans_ext(basename(newick_file)), "_tree.pdf")
  ggsave(
    filename = out_fname_pdf,
    plot = beautiful_tree,
    width = 10,
    height = 10,
    dpi = 300
  )
}

