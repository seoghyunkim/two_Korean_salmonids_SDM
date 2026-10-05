================

-   [Article information](#article-information)
-   [File descriptions](#file-descriptions)
-   [License](#license)
-   [Session information](#session-information)


# Article information
Title: Species-specific responses to climate change in two salmonids in South Korea” for consideration in Ecology of Freshwater Fish

Author: Jae-woo Joo, Jeongwon Lee,  Myeong-Hun Ko, Seoghyun Kim

Corresponding author email: seoghyunkim@kangwon.ac.kr


# Description of the data and file structure
- The dataset includes 9 files: 2 R code files for each salmonid species, and 6 GeoPackage files.

1. Files
The GeoPackage files consist of 2 basin boundary files (Basin_boundary.gpkg and Standard_basin_boundary.gpkg) and occurrence segment file (O_masou_presence_stream_segments.gpkg).
Streamnetwork file can be downloaded in following link: https://drive.google.com/drive/folders/1rt-sx3vcr8tYVxwhLbxf1gFIJ3IvxvZ6 (Seoghyun Kim's Google Drive)
Occurrence coordinates of endangered species are not publicly shared due to the data sensitivity and species protection concerns.

3. Code/software
Ensemble species distribution models were developed for each species using six algorithms: a maximum entropy (MaxEnt), random forest (RF), gradient boosting machine (GBM), artificial neural network (ANN), extreme gradient boosting (XGBoost), and light gradient boosting machine (LightGBM).
All analyses were conducted in R version 4.4.1 (R Core Team 2024). MaxEnt models were tuned and fitted using R package ENMeval (version 2.0.5.2) and dismo (version 1.3-16). RF, GBM, and ANN models were fitted using R package caret (version 6.0-94), ranger (version 0.16.0), gbm (version 2.2.2), and nnet (version 7.3-19). XGBoost and LightGBM models were fitted using R package xgboost (version 3.2.1.1) and lightgbm (version 4.6.0), respectively. An area under the curve (AUC) and maximized the sum of sensitivity and specificity (maxSSS) metrics were calculated using pROC (version 1.18.5) and ecospat (version 4.1.0).

4. Access information
Occurrence data were compiled from the National Ecosystem Survey datasets provided by the National Institute of Ecology and from the National Aquatic Ecological Monitoring Program provided through the Water Environment Information System.
