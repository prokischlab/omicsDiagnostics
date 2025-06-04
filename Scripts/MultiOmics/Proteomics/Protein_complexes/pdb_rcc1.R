
# Load required libraries
library(bio3d)
library(rgl)
library(dplyr)

# Step 1: Download the PDB structure (6G2J)
pdb <- read.pdb("6G2J")

# Step 2: Sample UniProt-PDB mapping (you will use your own data)
# Example: protein names mapped to chain IDs in 6G2J PDB file
uniprot_pdb_mapping <- data.frame(
  protein_name = c("NDUFS1", "NDUFV2", "NDUFB8", "ND1", "ND2", "NDUFA9"),
  uniprot_id = c("P22303", "P19404", "P19408", "P03886", "P03891", "O95182"),
  pdb_id = "6G2J",  # PDB ID for Complex I
  chain_id = c("G", "H", "A", "M", "K", "L")  # Corresponding chain IDs
)

# Step 3: Load your fold-change data
df_proteins <- data.frame(
  protein_name = c("NDUFS1", "NDUFV2", "NDUFB8", "ND1", "ND2", "NDUFA9"),
  fold_change = c(-2.5, 1.3, -0.7, 0.5, -1.2, 2.4)  # Example fold changes
)

# Step 4: Merge the fold change data with the UniProt-PDB mapping
df_merged <- merge(df_proteins, uniprot_pdb_mapping, by = "protein_name")

# Step 5: Create a color gradient based on fold change values
colors <- colorRampPalette(c("blue", "white", "red"))(100)  # Blue for downregulation, red for upregulation

# Normalize fold changes to scale between 1 and 100
fold_change_scaled <- as.numeric(cut(df_merged$fold_change, breaks = 100))

# Function to map fold change to color
get_color <- function(fold_change) {
  colors[as.numeric(cut(fold_change, breaks = 100))]
}

# Step 6: Assign colors to chains based on fold changes
chain_colors <- rep("gray", length(pdb$atom$resno))  # Initialize with neutral gray color

# Color the chains based on the fold change values
for (i in 1:nrow(df_merged)) {
  chain_id <- df_merged$chain_id[i]  # Get the chain ID from mapping
  fold_change_value <- df_merged$fold_change[i]  # Get the fold change value
  chain_color <- get_color(fold_change_value)  # Get the color for the fold change
  
  # Apply the color to the atoms of the corresponding chain
  chain_atoms <- pdb$atom$chain == chain_id
  chain_colors[chain_atoms] <- chain_color
}

# Step 7: Visualize the structure using rgl (interactive 3D plot)
# Convert the chain colors to a format rgl understands
chain_colors_rgb <- col2rgb(chain_colors)/255

# Plot the protein structure with colored chains
plot3d(pdb, type = "spheres", col = chain_colors, main = "Complex I Structure with Fold Changes")

# Step 8: Add a color legend (optional)
legend3d("topright", legend = c("Downregulated", "No Change", "Upregulated"),
         fill = c("blue", "white", "red"), cex = 0.8)
