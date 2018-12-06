# scRNA-seq_workflow

## How to use the scripts in this repository: 
Assuming you have docker installed already on your computer. Follow these steps to run this. Open up command line (Terminal or what have you). 

### 1. Build a docker image with the Dockerfile and create a container
``` bash
$ docker build -< Dockerfile -t <DESIRED_IMAGE_TAG_HERE>
```
### 2. Run a container with this image
```bash
$ docker run -it --rm --mount type=volume,dst=/home/rstudio/kitematic,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=<PUT_DESIRED_LOCAL_DIRECTORY_PATH_HERE> -e PASSWORD=<DESIRED_PASSWORD_HERE> -p 8787:8787 <SAME_DESIRED_IMAGE_TAG_AS_ABOVE_HERE>
```
### 3. Go to the command line in your container: 
Run this so you can find out what your container id is.
```bash
$ docker ps
$ docker exec -it abcde12345 bash
```
### 4. Run the pipeline to process the data with salmon and create a data matrix.RDS file 
$ cd <PATH_TO_THE_CLONED_REPOSITORY>
$ bash run_pipeline.sh

### 5. Open up Rstudio and run Rmd with the post processing analyses you wish to run
Go to your internet browser and enter: `localhost:8787`
In rstudio, open up one of these and follow along: 

- *seurat_data_processing.Rmd* - basic Rmd taking you through the beginning steps of seurat pipeline QC
- *tSNE_and_PCA.Rmd* - Runs tSNE and PCA and labels it with metadata to do initial clustering analysis
- *asap_data_prep.Rmd* - Makes data into a format that ASAP will take. [ASAP](https://ASAP.epfl.ch) has a gui with a pipeline with a lot of different options. 
