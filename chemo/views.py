from django.shortcuts import render
from django.views import generic
from django.conf import settings

# RDkit imports
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors

# Model import
from chemo.models import *

# General imports
import os
import pandas as pd

# Create your views here.

def loadSDF(sdfPath):

    # usar esse
    suppl = Chem.SDMolSupplier(sdfPath)
    for mol in suppl:
        try:
            # Smiles data
            smiles = Chem.MolToSmiles(mol)

            # Formal charge
            fCharge = Chem.GetFormalCharge(mol)

            # Mol weight
            molMass = Descriptors.ExactMolWt(mol)

            # clogp
            clogp = Crippen.MolLogP(mol)

            # tpsa
            tpsa = rdMolDescriptors.CalcTPSA(mol)

            # rotable
            rotable = rdMolDescriptors.CalcNumRotatableBonds(mol)

            print(mol.GetNumAtoms())
        except:
            print(-1)

def updateCountries():
    # Load countries.csv
    countries = pd.read_csv(os.path.join(settings.FILES_DIR, 'database/countries/countries.csv'))

    # For each record in csv
    for index, row in countries.iterrows():
        # Add it to database if doesn't exists yet
        obj, created = Countries.objects.get_or_create(
            name=row['Country'].strip(),
            continentName=row['Continent'].strip(),
        )

class IndexView(generic.ListView):
    """
    Class to work with index.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        loadSDF(os.path.join(settings.FILES_DIR, 'ApprovedDrugs2015.sdf'))

        # Render page index.html within the request and variables
        return render(request, 'chemo/index.html', {'countries': updateCountries()})