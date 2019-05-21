#!/usr/bin/bash

#
#   This section checks if channels are installed and if not
#   install them
#

channels=$(conda config --get channels)

if [[ $channels != *"--add channels 'conda-forge'"* ]]; then
  conda config --add channels 'conda-forge'
fi

if [[ $channels != *"--add channels 'rdkit'"* ]]; then
  conda config --add channels 'rdkit'
fi

{
    #
    #  Try to install the environment with the yml file
    #  if env already exists will fail
    #
    conda env create -f environment.yml --name chemdb
} || {
    #
    #  ligmol env does exists, so lets activate it
    #  and install some stuff
    #
    conda activate chemdb
    {
        #
        #  Install via `conda` directly.
        #  This will fail to install all
        #  dependencies. If one fails,
        #  all dependencies will fail to install.
        #
        conda install --yes --file requirements.txt
    } || {
        #
        #  To go around issue above, one can
        #  iterate over all lines in the
        #  requirements.txt file.
        #
        while read requirement; do conda install --yes $requirement; done < requirements.txt
    }
}
