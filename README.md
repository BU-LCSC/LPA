# LPA
Landsat Phenology Algorithm

The scripts provided here allow users to apply the Landsat Phenology Algorithm across Landsat ARD tiles. Here is a brief description of each script:

01_untar_ARD.R - Untars Landsat ARD files and generate directories for further processing 

02_topocorr_v3.R - Topographic correction functions 

03_calc_vi.R - Generates EVI2 GeoTIFF for each image in the ARD stack and applies topographic correction

04_landsat_pheno_pixel.R - Landsat Phenology algorithm

05_calc_phenology.R - Calculates phenometrics using the Landsat Phenology algorithm 

06_load_phenology.R - Ggenerates maps of annual/long-term average spring and autumn phenology 
for each Landsat ARD tile