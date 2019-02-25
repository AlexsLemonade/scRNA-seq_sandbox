# scRNA-seq_workflow

## About the scripts in this repository:

### 1) Pre-processing pipeline:
For starting from scratch with fastq files, an experiment ID or set of SRA's you would like to process into a counts gene matrix file.    
### A) For Smart-seq2 data or other single cell data that is available on SRA and not droplet or tag-based.   
### B) For 10X Genomics data (or Drop-seq data).  
    
### 2) Post-processing pipeline:
For starting from counts gene matrix data file that you would like to normalize and do further analyses on.

## Step 0: Set up docker image
For either pipeline, first clone the repository and then set up the docker image to work from:
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
Run the first line so you can find out what your container id is.
```bash
$ docker ps
```
It will be something like "a1b23c45" (a jumble of lower case letters and numbers). And you'll put that here:
```bash
$ docker exec -it <CONAINER_ID> bash
```

## 1) Pre-processing pipeline:
For starting from scratch with an experiment ID or set of SRA's you would like to process into a counts gene matrix file.

### A) For Smart-seq2 data or other single cell data that is available on SRA and not droplet or tag-based.

#### Step A1. Open `run_pre-processing_pipeline.sh` and change the variables to the dataset you
are working with OR keep this the same and follow this example's dataset.
```
# Change your directory name, GEO ID, and SRP here. Then run the script.
dir=darmanis_data
GSE=GSE84465
SRP=SRP079058
label=darmanis
```

#### Step A2. Run the pipeline to process the data with salmon and tximport to create a gene matrix file
Depending on how many samples are in the dataset this will take an hour or days (if you have thousands of samples)

#### Step A3. Open up Rstudio and prep the data for post-processing
To open Rstudio in docker, go to your internet browser and enter: `localhost:8787`
Follow the example in `darmanis_data_prep.Rmd` to set up data.

```bash
$ cd <PATH_TO_THE_CLONED_REPOSITORY>  
$ bash run_pre-processing_pipeline.sh
```
### B) For 10X Genomics data (or Drop-seq data)

#### Step B1. Open `run_tag-based_pre-processing_pipeline.sh` and change the url and variables to the dataset you
are working with OR keep this the same and follow this example's dataset.
```
# Change your directory name, and label here.
dir=pbmc_data
label=pbmc_1k_v2
```

Change this line to the url of the fastq files for the dataset you want to work with. Or keep as is and follow the example
```
cd ${dir}
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_fastqs.tar
tar -xvf ${label}.tar
```

#### Step B2. Run the pipeline to process the data with Alevin to create a gene matrix file
Depending on how many samples are in the dataset this will take a few hours or so.


## 2) Post-processing pipeline:
For starting from counts gene matrix data file that you would like to normalize and do further analyses on.

### Step 1.  Open up Rstudio and prep the data for post-processing
To open Rstudio in docker, go to your internet browser and enter: `localhost:8787`
Follow the example in `darmanis_data_prep.Rmd` to set up data.

### Step 2.  Open `run_post-processing_pipeline.sh` and change the variables to the dataset you
are working with OR keep this the same and follow this example's dataset.
```bash
dir=darmanis_data
label=darmanis
```

### Step 3. Run the pipeline
```bash
$ cd <PATH_TO_THE_CLONED_REPOSITORY>  
$ bash run_post-processing_pipeline.sh
```
