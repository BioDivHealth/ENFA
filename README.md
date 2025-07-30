# Environmental Niche Factor Analysis (ENFA)

**Environmental Niche Factor Analysis (ENFA)** is a powerful statistical method used in ecology and conservation biology to model species distributions using **presence-only data**. Unlike traditional methods that require both presence and absence data, ENFA leverages environmental variables to characterize the ecological niche of a species and assess habitat suitability across a landscape.

ENFA is particularly valuable in conservation contexts where absence data is unreliable or unavailable—such as for rare, elusive, or poorly surveyed species. It helps researchers and conservationists:

- **Identify critical habitats** for endangered species.
- **Assess the impact of environmental change** on species distributions.
- **Support reserve design and spatial planning** by highlighting areas of high suitability.
- **Compare niche characteristics** across species or populations.

By calculating the **marginality** (how different the species' niche is from the average environment) and **specialization** (how restricted the species is to specific conditions), ENFA provides a nuanced view of species–environment relationships. The resulting **suitability maps** can guide conservation priorities and inform ecological research.

## How ENFA Works

ENFA operates by comparing the environmental conditions at species presence locations to the overall environmental conditions in the study area. It calculates:

- **Marginality**: The difference between the mean environmental conditions of the species and the global mean.
- **Specialization**: The ratio of variance in the global environment to the variance in the species environment.
- **Mahalanobis distance**: A multivariate measure of distance from each location to the species niche centroid.

In this implementation, **suitability** is derived using the **inverse of the chi-square distribution** applied to the Mahalanobis distance matrix. This transformation allows for a probabilistic interpretation of habitat suitability across the landscape.

## Workflow Summary

### 1. Load Species Spatial Data

- **Presence records**: Loaded from CSV files.
- **Species range**: Loaded from IUCN shapefiles.
- **Study area**: Defined by buffering the species range and selecting intersecting countries.

### 2. Load Environmental Data

- Uses WorldClim bioclimatic variables at 2.5 arc-minute resolution.
- Cropped to the study area.

### 3. Clean Records

- Filters presence points to those within the IUCN range.

### 4. Prepare Data for ENFA

- Extracts environmental values for presence and background locations.
- Constructs a presence/background index.

### 5. Run ENFA

- Calculates niche centroid, marginality, specialization, and suitability.
- Outputs raster layers for visualization.

### 6. Visualize Results

- Plots:
  - Mahalanobis distance to niche centroid
  - Marginality
  - Specialization
  - Suitability

## Main Outputs

- **ENFA raster layers**:
  - `Mahalanobis_dist`
  - `Marginality`
  - `Specificity`
  - `Suitability`
- **Niche centroid coordinates**
- **Marginality loadings**

## Dependencies

- `sf`, `terra`, `geodata`, `viridis`, `dplyr`, `ggplot2`
- Custom functions: `ENFA_function`, `plot_enfa`

## How to Run

1. Clone the repository.
2. Place species records and range shapefiles in `./Data/Sp_info/Records` and `./Data/Sp_info/Ranges`.
3. Run the R script or R Markdown file.
4. Outputs will be saved in `./Data/Env_vars` and visualized in the notebook.

## Author

**Gonzalo Albaladejo Robles**  
Bioinformatician, Wellcome Henry Dale Fellow
