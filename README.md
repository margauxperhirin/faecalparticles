# faecalparticles

Scripts and raw data for the manuscript "Could faecal pellets be identified from *in situ* images in the Arctic?"

Authors : Margaux Perhirin<sup>1</sup>*, Laure Vilgrain<sup>2</sup>, Marc Picheral<sup>3</sup>, Frédéric Maps<sup>4</sup>, Sakina-Dorothée Ayata<sup>1,5</sup> 

<sup>1</sup> Sorbonne Université, CNRS, IRD, MNHN, Laboratoire d’Océanographie et du Climat : Expérimentation & Approches Numériques,  LOCEAN-IPSL, Paris, France; 

<sup>2</sup> Memorial University of Newfoundland, Centre for Fisheries Ecosystems Research (CFER), Canada;

<sup>3</sup> Sorbonne Université, CNRS, Laboratoire d’Océanographie de Villefranche (LOV), Villefranche-sur-Mer, France; 

<sup>4</sup> Takuvik Joint International Laboratory Université Laval-CNRS, Département de Biologie and Québec-Océan, Université Laval, Québec, Canada; 

<sup>5</sup> Institut Universitaire de France, Paris, France.

* corresponding author margaux.perhirin@locean.ipsl.fr 

## How does it work?

On this GitHub you can find :
- Data : all the raw and transformed data needed for the scripts to be run
    - ecotaxa_export_149_20231221_1029.tsv (associated images can be found on EcoTaxa web application https://ecotaxa.obs-vlfr.fr/, project ‘UVP5hd GreenEdge 2016 [149]’, access authorisation required)
    - taxon_ecotaxa.csv
    - zoo_env_greenedge_final.csv
    - Arctic_Driftingtraps25m.csv
    - faecalpellets_inddata.csv
    - faecalpellets_incubations_mm.csv
    - poopingcopepods_mm.csv
    - env_stations.csv.gz
    - export_detailed_20240404_15_40_PAR_Aggregated.tsv
    - faecalpellets_copepods-sizereview.csv
    - faecalpellets_zoo-sizereview
- Scripts : all the R scripts to run the analyses presented in the paper
    - 1.fromEcoTaxafile.R : to clean the file from EcoTaxa, to save the morphological features, to trim the images and to save them
    - 2.watervolume.R : to get the watervolume sampled at each station per 10-m bins
    - 3.detritus_descriptorsandadditionalvariables.R : to BoxCox-transformed the morphological descriptors and to compute additional variables
    - 4.detritus_morphospace&kmeans.R : to build the morphological space and to find the best clusters
    - 5.faecalparticles_sizecomparison.R : to compare lengths distribution
    - 6.allparticles_sizedistribution.R : to compare the coverage of each group of faecal pellets/particles with UVP5 size range, and data from literature


Please cite the manuscript if using this code, partially or in its totality.
