# You'll need a docker login and a AWS login.

Follow this to make a EC2 instance.
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html

### Set these settings for the EC2 instance.
- Make with ubuntu
- Free tier 64bit (x86)
- t3.xlarge
- Generated pem
- 100gb gp2 storage

Change permissions for pem key-pair
```bash
chmod 400 ~/Downloads/newkey2.pem
```

ssh -i <path to pem> ubuntu@<public ip4>
```bash
ssh -i ~/Downloads/newkey2.pem ubuntu@34.220.46.119
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
sudo singularity build --sandbox bc_multiome_test/ docker://ubuntu:latest
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
    export PATH=/opt/miniconda3/bin:$PATH
    export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda install -y -c conda-forge mamba 
	conda config --add channels bioconda
	conda config --add channels conda-forge

	#install python libraries
	pip install macs3 #
	pip install scrublet # 
	pip install scipy #
	pip install matplotlib # 
	pip install numpy #
	pip install pandas #
	pip install h5py #
	pip install tables #
	pip install annoy==1.15 #
	
	# denotes installed in sandbox without error
	#install additional tools
	mamba install -y -f bioconda::bwa #
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f bioconda::fastqc #
	mamba install -y -f bioconda::multiqc #
	mamba install -y -f anaconda::graphviz #
	mamba install -y -f conda-forge::parallel #
	conda install -y -f conda-forge::ncurses #

	#install R packages
	conda install -y -f r-base=4.2 #
	mamba install -y -f conda-forge::r-devtools #
	mamba install -y -f conda-forge::r-biocmanager=1.30.19 #
	mamba install -y -f conda-forge::r-rlang #
	mamba install -y -f conda-forge::r-ggplot2 #
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
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	#Bioconductor packages through conda
	mamba install -y -f bioconda::bioconductor-biocparallel #
	mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
	mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #
	mamba install -y -f bioconda::bioconductor-jaspar2020 #
	mamba install -y -f bioconda::bioconductor-org.hs.eg.db #
	mamba install -y -f bioconda::bioconductor-tfbstools #
	mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene #
	mamba install -y -f bioconda::bioconductor-universalmotif #
	mamba install -y -f bioconda::bioconductor-chromvar #
	mamba install -y -f bioconda::bioconductor-motifmatchr #
	mamba install -y -f bioconda::bioconductor-decoupler #
	mamba install -y -f bioconda::bioconductor-scran #
	mamba install -y -f bioconda::bioconductor-infercnv #
	mamba install -y -f bioconda::bioconductor-complexheatmap #
	mamba install -y -f bioconda::bioconductor-biovizbase #
	mamba install -y -f bioconductor-glmgampoi #
'GenomeInfoDb', 'GenomicRanges', 'IRanges', 'Rsamtools', 'S4Vectors', 'BiocGenerics' 
	#install cistopic
	R --slave -e 'devtools::install_github("aertslab/RcisTarget")' #
	R --slave -e 'install.packages("lda", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("doSNOW", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("DT", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("feather", repos="http://cran.us.r-project.org")' #
	cd #
	wget https://github.com/aertslab/cisTopic/archive/refs/tags/v2.1.0.tar.gz 
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
	R --slave -e 'devtools::install_github("SydneyBioX/scDC")' #

	#Correct matrix version
	mamba install -y -f conda-forge::r-matrix=1.6_1 # set this version
	mamba install -y -f r::r-irlba # set this version
	R --slave -e 'install.packages("Matrix", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
	R --slave -e 'install.packages("irlba", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
	R --slave -e 'oo <- options(repos = "https://cran.r-project.org/");
		tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]];
		install.packages("lme4", type = "source",repos = "http://cran.us.r-project.org");
		options(oo)'





#Changelog v0.2:
#added changed r irlba and rmatrix installation due to an error in newer matrix installs. 
#added genefu for pseudobulk pam50 analysis

%labels
    Author Ryan Mulqueen
    Version v0.2
    MyLabel Multiome Breast Cancer Processing


```

```bash
sudo singularity build multiome_bc.sif multiome_bc.def
sudo singularity shell multiome_bc.sif
```

## Pull Images
Use sftp to get images off cluster
Now test on exacloud
```bash
sftp mulqueen@acc.ohsu.edu
cd /home/groups/CEDAR/mulqueen/bc_multiome
put multiome_bc.sif
#also using Demuxafy.sif for scrublet
cd /home/groups/CEDAR/mulqueen/bc_multiome
wget -O Demuxafy.sif 'https://www.dropbox.com/scl/fi/g0cuyjwomdavom6u6kb2v/Demuxafy.sif?rlkey=xfey1agg371jo4lubsljfavkh&'

```



## Building scrublet image
multiome_scrub.def

```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/opt/miniconda3/bin:$PATH
	echo "set enable-bracketed-paste off" >> ~/.inputrc
%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda remove menuinst
	conda install -y python=3.7

	#install python libraries
	pip install scipy #
	pip install matplotlib # 
	pip install numpy==1.19.4  #
	pip install pandas #
	pip install h5py #
	pip install tables #
	pip install argparse #
	pip install scrublet # 
	pip install annoy==1.15.2 #

%labels
    Author Ryan Mulqueen
    Version v0.1
    MyLabel Multiome Breast Cancer Processing Scrublet

```
```bash
sudo singularity build multiome_scrub.sif multiome_scrub.def

sudo singularity build --sandbox scrub docker://ubuntu:latest
sudo singularity shell --writable scrub/
```


## Building Harmony Image

```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/opt/miniconda3/bin:$PATH
    export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda install -y -c conda-forge mamba 
	conda config --add channels bioconda
	conda config --add channels conda-forge

	#install R packages
	mamba install -y -f conda-forge::parallel #
	conda install -y -f conda-forge::ncurses #
	conda install -y -f r-base #
	mamba install -y -f conda-forge::r-devtools #
	mamba install -y -f conda-forge::r-biocmanager#
	mamba install -y -f conda-forge::r-rlang #
	mamba install -y -f conda-forge::r-ggplot2 #
	mamba install -y -f bioconda::bioconductor-genomeinfodb
	mamba install -y -f bioconda::bioconductor-genomicranges
	mamba install -y -f bioconda::bioconductor-iranges
	mamba install -y -f bioconda::bioconductor-rsamtools
	mamba install -y -f bioconda::bioconductor-s4vectors
	mamba install -y -f bioconda::bioconductor-biocgenerics

	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
	R --slave -e 'devtools::install_github("stuart-lab/signac", "develop")' #
	R --slave -e 'install.packages("harmony",repos = "http://cran.us.r-project.org")'
	R --slave -e 'install.packages("ggalluvial",repos = "http://cran.us.r-project.org")'

	#Correct matrix version
	mamba install -y -f conda-forge::r-matrix=1.6_1 # set this version
	mamba install -y -f r::r-irlba # set this version
	R --slave -e 'install.packages("Matrix", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
	R --slave -e 'install.packages("irlba", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
	R --slave -e 'install.packages("SeuratObject", type = "source",repos="http://cran.us.r-project.org")' #reinstall from source
	R --slave -e 'oo <- options(repos = "https://cran.r-project.org/");
		tools::package_dependencies("Matrix", which = "LinkingTo", reverse = TRUE)[[1L]];
		install.packages("lme4", type = "source",repos = "http://cran.us.r-project.org");
		options(oo)'

%labels
    Author Ryan Mulqueen
    Version v0.1
    MyLabel Multiome Breast Cancer Processing Harmony

```
```bash
sudo singularity build harmony.sif harmony.def

sudo singularity shell harmony.sif
```



## Building NMF Image

```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/opt/miniconda3/bin:$PATH
    export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda install -y -c conda-forge mamba 
	conda config --add channels bioconda
	conda config --add channels conda-forge

	#install python libraries
	pip install macs3 #
	pip install scrublet # 
	pip install scipy #
	pip install matplotlib # 
	pip install numpy #
	pip install pandas #
	pip install h5py #
	pip install tables #
	pip install annoy==1.15 #
	
	# denotes installed in sandbox without error
	#install additional tools
	mamba install -y -f bioconda::bwa #
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f bioconda::fastqc #
	mamba install -y -f bioconda::multiqc #
	mamba install -y -f anaconda::graphviz #
	mamba install -y -f conda-forge::parallel #
	conda install -y -f conda-forge::ncurses #

	#install R packages
	conda install -y -f r-base=4.2 #
	mamba install -y -f conda-forge::r-devtools #
	mamba install -y -f conda-forge::r-biocmanager=1.30.19 #
	mamba install -y -f conda-forge::r-rlang #
	mamba install -y -f conda-forge::r-ggplot2 #
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
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	#Bioconductor packages through conda
	mamba install -y -f bioconda::bioconductor-biocparallel #
	mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
	mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #
	mamba install -y -f bioconda::bioconductor-jaspar2020 #
	mamba install -y -f bioconda::bioconductor-org.hs.eg.db #
	mamba install -y -f bioconda::bioconductor-tfbstools #
	mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene #
	mamba install -y -f bioconda::bioconductor-universalmotif #
	mamba install -y -f bioconda::bioconductor-chromvar #
	mamba install -y -f bioconda::bioconductor-motifmatchr #
	mamba install -y -f bioconda::bioconductor-decoupler #
	mamba install -y -f bioconda::bioconductor-scran #
	mamba install -y -f bioconda::bioconductor-complexheatmap #
	mamba install -y -f bioconda::bioconductor-biovizbase #
	mamba install -y -f bioconductor-glmgampoi #
	mamba install -y -f bioconductor-ucell #
	mamba install -y -f bioconductor-fgsea
	
	#install cistopic
	R --slave -e 'devtools::install_github("aertslab/RcisTarget")' #
	R --slave -e 'install.packages("lda", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("doSNOW", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("DT", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("feather", repos="http://cran.us.r-project.org")' #
	cd #
	wget https://github.com/aertslab/cisTopic/archive/refs/tags/v2.1.0.tar.gz 
	R --slave -e 'install.packages("v2.1.0.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

	#Funner stuff!
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("Signac", repos="http://cran.us.r-project.org")' #
	R --slave -e 'remotes::install_github("satijalab/seurat-wrappers")' #
	R --slave -e 'devtools::install_github("JuliusCampbell/TITAN")' #
	R --slave -e 'install.packages("GeneNMF", type = "source",repos="http://cran.us.r-project.org")' 
	R --slave -e 'install.packages("RcppML", type = "source",repos="http://cran.us.r-project.org")' 
	R --slave -e 'install.packages("msigdbr", type = "source",repos="http://cran.us.r-project.org")' 

#Changelog v0.2:
#added changed r irlba and rmatrix installation due to an error in newer matrix installs. 
#added genefu for pseudobulk pam50 analysis

%labels
    Author Ryan Mulqueen
    Version v0.2
    MyLabel Multiome Breast Cancer Processing with GeneNMF

```

```bash
singularity build --fakeroot multiome_nmf.sif multiome_nmf.def


sftp mulqueen@acc.ohsu.edu
cd /home/groups/CEDAR/mulqueen/bc_multiome
put multiome_scrub.sif

singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome multiome_scrub.sif
cd /home/groups/CEDAR/mulqueen/bc_multiome/cellranger_data/third_round/IDC_4/outs
#test run

```


# Making a container for deeptools and SCENIC


scenic.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/opt/miniconda3/bin:$PATH
    export PYTHONPATH=/opt/miniconda3/lib/python3.11/:$PYTHONPATH

%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	nano \
	bedtools

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda install -y python==3.11.8
	conda config --add channels bioconda
	conda config --add channels conda-forge
	conda install -y -c bioconda pysam

	#install deeptools
	pip install deeptools

	#install scenicplus
	git clone https://github.com/aertslab/scenicplus
	cd scenicplus
	pip install .

	pip install h5py numpy argparse pybedtools pandas scipy tables
	pip install numpy scipy cython numba matplotlib scikit-learn h5py click
	pip install velocyto

%files
	/home/ubuntu/Mallet-202108 /container_mallet

%labels
    Author Ryan Mulqueen
    Version v0.2
    MyLabel SCENIC+ and DeepTools Env
```

```bash
#predownload mallet to include in container
wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
tar -xf Mallet-202108-bin.tar.gz
sudo singularity build scenicplus.sif scenic.def
```