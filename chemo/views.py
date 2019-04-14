from django.shortcuts import render
from django.views import generic
from django.conf import settings

# RDkit imports
from rdkit import Chem

# Model import
from chemo.models import *

# General imports
import os
import pandas as pd

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

# Create your views here.

def loadSDF(sdfPath):

    # usar esse
    suppl = Chem.SDMolSupplier(sdfPath)
    for mol in suppl:
        try:
            print(mol.GetNumAtoms())
        except:
            print(-1)


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