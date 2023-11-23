# cgevirus (k-mer gene typer)

# Installation

The following Conda channels are required:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

Mamba can be installed with:

`conda install -c conda-forge mamba`

For a fast install of cgevirus, use mamba:
`mamba install -c genomicepidemiology cgevirus`


# Datebase

Download the cge_db database:

`wget https://cge.food.dtu.dk/services/MINTyper/cge_db.tar.gz`

`tar -xvzf cge_db.tar.gz`

`sudo mkdir -777 /opt/cge`

`mv cge_db /opt/cge/.`


# Usage

Standard Usage:
`cgevirus -i <input_file> -o <output_file>`

If you have a folder of many fastq.gz files:
`cgevirus -f <input_folder> -name <sample_name> -o <output_folder>`

Help Message:
`cgevirus -h`
