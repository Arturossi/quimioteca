from django.conf import settings
from django.core.paginator import Paginator
from django.shortcuts import render, redirect
from django.views import generic

# RDkit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import inchi
from rdkit.Chem import rdMolDescriptors

# Model import
from chemo.models import *

# General imports
import os
import json
import math
import logging

import pandas as pd

# Logger to log
logger = logging.getLogger(__name__)

# Create your views here.

# Function to load and parse .sdf files into database [UNFINISHED]
def loadSDF(sdfPath):

    # Create images
    generateImages(sdfPath)
    
    # Create a molecule supplier
    suppl = Chem.SDMolSupplier(sdfPath)
    
    # Filter empty entries
    sdf = [x for x in suppl if x is not None]
    
    # For each molecule in supplier
    for mol in sdf:
        data = {}
        
        try:
            data['fCharge'] = Chem.GetFormalCharge(mol)
        except:
            data['fCharge'] = None
            
        try:
            data['name'] = mol.GetProp('_Name') # Get molecule name
        except:
            data['name'] = 'unkown'
            
        try:
            # pka
            data['molMass'] = Descriptors.ExactMolWt(mol) # Mol weight ???
            # monoisotripic mass
            data['cLogP'] = Crippen.MolLogP(mol) # clogp ???
            data['tpsa'] = rdMolDescriptors.CalcTPSA(mol)
            data['hbond'] = rdMolDescriptors.CalcNumHBA(mol)
            data['hbondDonnors'] = rdMolDescriptors.CalcNumHBD(mol)
            data['rotable'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
            data['smiles'] = Chem.MolToSmiles(mol)
            data['InChI'] = inchi.MolToInchi(mol)
            data['inchiKey'] = inchi.MolToInchiKey(mol)
            data['numAtoms'] = mol.GetNumAtoms()
            data['Molecule Name'] = mol.GetProp('_Molecule Name')
            if Entry.objects.filter(inChIKey=linchiKey).exists():
                if not Entry.objects.filter(provider=lprovider).exists():
                    feedDatabase(data)
                else:
                    continue
                    
            else:
                feedDatabase(data)
                
        except:
            print("Molecule not processed")
            continue
            

def feedDatabase(data):
    try: # Try to parse (this avoid those badly filled entries e.g.: half empties)
        obj, created = Compounds.objects.get_or_create(
                moleculeName=data['name'].strip(),
                totalMolweight=0.0,#float(data['Total Molweight']),
                cLogP=float(data['cLogP']),
                cLogS=0.0,#float(data['cLogS']),
                hAcceptors=0,#int(data['H-Acceptors']),
                hDonors=int(data['hbondDonnors']),
                totalSurfaceArea=0.0,#float(data['Total Surface Area']),
                polarSurfaceArea=0.0,#float(data['Polar Surface Area']),
                mutagenic='none',#data['Mutagenic'].strip(),
                tumorigenic='none',#data['Tumorigenic'].strip(),
                irritant='none',#data['Irritant'].strip(),
                nonHAtoms=0,#int(data['Non-H Atoms']),
                stereoCenters=0,#int(data['Stereo Centers']),
                rotatableBonds=int(data['rotable']),
                smiles=data['smiles'].strip(),
                inChI=data['InChI'].strip(),
                inChIKey=data['inchiKey'].strip(),
                provider = data['provider'].strip(),
                numAtoms = int(data['numAtoms']),
                molname = data['Molecule Name'].strip(),
            ) # Try to fetch elements, if fails insert elements into database (make database unique)
    except:
        logger.warn("The molecule " + data['name'].strip() +
                    " is experiencig some problems, skipping it.") # Show a warning


def updateCountries():
    '''
        Update the Countries table in database from file in 'database/countries/countries.csv' path
    '''

    # Load countries.csv
    countries = pd.read_csv(os.path.join(
        settings.FILES_DIR, 'database/countries/countries.csv'))

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

        #loadSDF(os.path.join(settings.FILES_DIR, 'ApprovedDrugs2015.sdf'))

        # Render page index.html within the request and variables
        return render(request, 'chemo/index.html', {'countries': updateCountries()})


class compoundView(generic.ListView):
    """
    Class to work with compound.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        return render(request, 'chemo/compound.html')


class aboutView(generic.ListView):
    """
    Class to work with about.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        return render(request, 'chemo/about.html')


class contactUsView(generic.ListView):
    """
    Class to work with contactUs.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        return render(request, 'chemo/contactUs.html')


class loginView(generic.ListView):
    """
    Class to work with login.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        return render(request, 'chemo/login.html')
    
class cadastroMolView(generic.ListView):
    """
    Class to work with login.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        return render(request, 'chemo/cadastroMoleculas.html')


class quimiotecaDatabaseView(generic.ListView):
    """
    Class to work with quimiotecaDatabase.html template
    """

    def get(self, request, **kwargs):
        """
        Get function to the class
        """

        mols = Compounds.objects.all() # Read database

        count = mols.count() # Count number of elements on it

        if request.method == 'GET': # If there is any GET response
            if 'elements' in request.GET: # If in GET a 'elements' field has been passed
                try:
                    elements = int(request.GET.get('elements')) # Try to parse its value to int (it comes as string, and if it fails, it is an invalid value)

                    if not elements or elements not in [10, 25, 50, 100]: # If elements are not one of those in list (avoid the get injection)
                        elements = 10 # Set elements to default value
                except:
                    elements = 10 # Set elements to default value
            else:
                elements = 10 # Set elements to default value

            if 'page' in request.GET: # If in GET a 'page' field has been passed
                try: 
                    page = int(request.GET.get('page')) # Try to parse its value to int (it comes as string, and if it fails, it is an invalid value)

                    if not page or page < 1 or page > (int(count) / elements + 1): # If page are not accepted range (avoid the get injection or type error)
                        page = 1 # Set page to default value
                except:
                    page = 1 # Set page to default value
            else:
                page = 1 # Set page to default value
        else:
            page = 1 # Set page to default value
            elements = 10 # Set elements to default value

        paginator = Paginator(mols, elements) # Make the pages
        mols = paginator.page(page) # Set the page
        endpage = math.ceil(count/elements) # Calculate how many pages we have


        nextpage = page + 1 # Value to be put in next '>' button

        previouspage = page - 1 # Value to be put in previous '<' button

        starting = (page-1)*elements + 1 # Starting element in table

        ending = page*elements # Ending element in table

        if ending > count: # Check if the last element is higher than the size of database
            ending = count # Correct the issue


        variables = {
            'mols': mols,
            'page': page,
            'next': nextpage,
            'previous': previouspage,
            'count': count,
            'starting': starting,
            'ending': ending,
            'endpage': endpage,
            'elements': elements,
        } # Set all variables

        return render(request, 'chemo/quimiotecaDatabase.html', variables)
