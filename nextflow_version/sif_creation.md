# You'll need a docker login and a AWS login.

Follow this to make a EC2 instance.
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html

### Set these settings for the EC2 instance.
- Make with ubuntu
- Free tier 64bit (x86)
- t3.micro
- Generated pem
- 50gb gp2 storage

Change permissions for pem key-pair
```bash
chmod 400 ~/Downloads/newkey.pem
```

ssh -i <path to pem> ubuntu@<public ip4>
```bash
ssh -i ~/Downloads/newkey.pem ubuntu@34.222.152.147
```


## Install singularity, go, and docker. Use docker to pull the ubuntu image
https://docs.sylabs.io/guides/latest/user-guide/quick_start.html
```bash
#install basics ##
# Ensure repositories are up-to-date
sudo apt-get update
# Install debian packages for dependencies
sudo apt-get install -y \
   autoconf \
   automake \
   cryptsetup \
   git \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev \
   make

## Install go ##
export VERSION=1.21.0 OS=linux ARCH=amd64 && \
  wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
  sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
  rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc && \
  source ~/.bashrc

## Install singularity ##
export VERSION=4.1.0 && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}

./mconfig && \
    make -C builddir && \
    sudo make -C builddir install



## Install docker and use to pull the ubuntu image. ##
#Convert ubuntu image to sif file for building our own SIFs
cd 
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh ./get-docker.sh
sudo chmod 666 /var/run/docker.sock

#Pull ubuntu image to act as our base ##
#docker pull ubuntu:latest
#docker create --name ubuntu -p 80:80 ubuntu:latest
#docker images #take IMAGE ID 
#sudo docker save fd1d8f58e8ae -o ubuntu.tar #save image as tar
#sudo singularity build ubuntu.sif docker-archive://ubuntu.tar #convert tar to sif to act as local bootstrap

```

## Multiome BC Container Building

```bash
sudo singularity build --sandbox bc_multiome_test/ /home/ubuntu/ubuntu.sif
sudo singularity shell --writable bc_multiome_test/

#test r installs and stuff (just go line by line of def below)

```

bc_multiome.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/miniconda3/bin:$PATH
    export PYTHONPATH=/miniconda3/lib/python3.9/:$PYTHONPATH

    # activate conda environment
    conda init
    source activate base;
    
%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake build-essential bzip2 wget git default-jre unzip

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/miniconda3/bin:$PATH"
	#build full conda environment for sif
	conda install -y -c conda-forge mamba 
	conda config --add channels bioconda
	conda config --add channels conda-forge

	#install python libraries
	pip install MACS2 #
	pip install scrublet #
	pip install scipy #
	pip install matplotlib #
	pip install numpy #
	pip install pandas #

	# denotes installed in sandbox without error
	#install additional tools
	mamba install -y -f bwa #
	mamba install -y -f samtools # 
	mamba install -y -f bedtools #
	mamba install -y -f fastqc #
	mamba install -y -f multiqc #
	mamba install -y -f conda-forge::parallel #

	#install R packages
	mamba install -y -f r-base=4.2 #
	mamba install -y -f r-devtools #
	mamba install -y -f r-biocmanager=1.30.19 #
	mamba install -y -f r-rlang #
	mamba install -y -f r-ggplot2 # 
	mamba install -y -f bioconda::bioconductor-dirichletmultinomial #
	mamba install -y -f conda-forge::r-igraph #
	mamba install -y -f conda-forge::r-rjags #
	mamba install -y -f conda-forge::r-leiden #
	mamba install -y -f conda-forge::r-hdf5r #
	mamba install -y -f conda-forge::r-rmpfr #
	mamba install -y -f conda-forge::r-ggraph #
	mamba install -y -f conda-forge::r-nloptr #
	mamba install -y -f conda-forge::r-jomo #

	#R utility libraries
	R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("Matrix", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	#Bioconductor packages through conda
	mamba install -y -f bioconda::bioconductor-biocparallel
	mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38
	mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86
	mamba install -y -f bioconda::bioconductor-jaspar2020
	mamba install -y -f bioconda::bioconductor-org.hs.eg.db
	mamba install -y -f bioconda::bioconductor-tfbstools
	mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
	mamba install -y -f bioconda::bioconductor-universalmotif
	mamba install -y -f bioconda::bioconductor-chromvar
	mamba install -y -f bioconda::bioconductor-motifmatchr
	mamba install -y -f bioconda::bioconductor-decoupler
	mamba install -y -f bioconda::bioconductor-scran
	mamba install -y -f bioconda::bioconductor-infercnv
	mamba install -y -f bioconda::bioconductor-complexheatmap

	#R --slave -e 'BiocManager::install("BiocParallel")' #
	#R --slave -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")' #
	#R --slave -e 'BiocManager::install("EnsDb.Hsapiens.v86")' #
	#R --slave -e 'BiocManager::install("JASPAR2020")' #
	#R --slave -e 'BiocManager::install("org.Hs.eg.db")' #
	#R --slave -e 'BiocManager::install("TFBSTools")' #
	#R --slave -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")' #
	#R --slave -e 'BiocManager::install("universalmotif")' #
	#R --slave -e 'BiocManager::install("chromVAR")' #
	#R --slave -e 'BiocManager::install("motifmatchr")' #
	#R --slave -e 'BiocManager::install("infercnv")' #
	#R --slave -e 'BiocManager::install("ComplexHeatmap")' #
	#R --slave -e 'BiocManager::install("decoupleR")' #
	#R --slave -e 'BiocManager::install("scran")' #

	#install cistopic
	R --slave -e 'devtools::install_github("aertslab/RcisTarget")' #
	R --slave -e 'devtools::install_github("aertslab/AUCell")' #
	R --slave -e 'install.packages("lda", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("doSNOW", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("DT", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("feather", repos="http://cran.us.r-project.org")' #
	cd #
	wget https://github.com/aertslab/cisTopic/archive/refs/tags/v2.1.0.tar.gz #
	R --slave -e 'install.packages("v2.1.0.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

	#Funner stuff!
	R --slave -e 'install.packages("rliger", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("SoupX", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("Signac", repos="http://cran.us.r-project.org")' #
	R --slave -e 'remotes::install_github("satijalab/seurat-wrappers")' #
	R --slave -e 'devtools::install_github("JuliusCampbell/TITAN")' #
	R --slave -e 'devtools::install_github("caleblareau/BuenColors")' #
	R --slave -e 'devtools::install_github("buenrostrolab/FigR")' #
	R --slave -e 'devtools::install_github("quadbio/Pando")' #
	R --slave -e 'install.packages(c("DescTools", "reshape2", "ggridges", "mice"), repos="http://cran.us.r-project.org")' #
	R --slave -e 'devtools::install_github("SydneyBioX/scDC")' 

	#TO ADD??
	#R --slave -e 'devtools::install_github("navinlabcode/copykat")'
	#R --slave -e 'remotes::install_github("akdess/CaSpER")' #https://rpubs.com/akdes/673120
	#copyscAT
	#AneuFinder


```

```bash
sudo singularity build multiome_bc.sif bc_multiome.def
sudo singularity shell multiome_bc.sif
```

## Pull Images
Use sftp to get images off cluster

```bash
sftp -i ~/Downloads/newkey.pem ubuntu@34.222.152.147
get *sif
```
Now test on exacloud/seadragon
```bash
sftp mulqueen@acc.ohsu.edu
put multiome_bc.sif

ssh mulqueen@acc.ohsu.edu
ssh exacloud.ohsu.edu
module load singularity
singularity shell multiome_bc.sif
cp ~/multiome_bc.sif /home/groups/CEDAR/mulqueen/bc_multiome

```