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
'''def loadSDF(sdfPath):

    # usar esse
    suppl = Chem.SDMolSupplier(sdfPath)
    for mol in suppl:
        try:
            # pka

            # Formal charge
            fCharge = Chem.GetFormalCharge(mol)
        except:
        
        try:
            fCharge = NULL
            # Mol weight ???
            molMass = Descriptors.ExactMolWt(mol)

            # monoisotripic mass


            # clogp ???
            clogp = Crippen.MolLogP(mol)

            # tpsa
            tpsa = rdMolDescriptors.CalcTPSA(mol)

            # hbond
            hbond = rdMolDescriptors.CalcNumHBA(mol)

            # hbondDonnord
            hbondDonnors = rdMolDescriptors.CalcNumHBD(mol)

            # rotable
            rotable = rdMolDescriptors.CalcNumRotatableBonds(mol)

            # Smiles data
            smiles = Chem.MolToSmiles(mol)

            # Inchi
            inchi = inchi.MolToInchi(mol)

            # Inchi Key
            inchiKey = inchi.MolToInchiKey(mol)

            print(mol.GetNumAtoms())
        except:
            print(-1)'''

def generateImages(path):
    '''
        From a .sdf file, generate images and save them in data/molImages/ being filename the name of the molecule and the format of the image .png

        VARIABLES:
            - path [string]: Path to .sdf file
    '''
    
    sdf = Chem.SDMolSupplier(path) # Load sdf file

    ms = [x for x in sdf if x is not None] # Filter empty entries

    # Not sure about the need of this piece of code VVVV

    for m in ms: # For each element in list
        tmp = AllChem.Compute2DCoords(m) # Compute its coordinates

    # Not sure about the need of this piece of code ^^^^

    for i in range(len(ms)): # Iterate over all elements
        name = ms[i].GetProp('_Name') # Get molecule name
        Draw.MolToFile(ms[i], 
            os.path.join(settings.FILES_DIR, f'molImages/{name}.png'),
            size=(300,300),
            kekulize=True, 
            wedgeBonds=True,
            fitImage=True) # Save it

        Draw.MolToFile(ms[i], 
            os.path.join(settings.FILES_DIR, f'molThumbs/{name}.png'),
            size=(150,150),
            kekulize=True,
            wedgeBonds=True,
            fitImage=True) # Save it
    return

'''def csv2json(request, path='GREENIDGE_APAGAR', separator='\t'):
    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf"))

    df = pd.read_csv(os.path.join(settings.FILES_DIR,
                                  str(path) + ".txt"), sep=separator)
    df2 = df.drop(['Structure [idcode]', 'Unnamed: 18'], axis=1)
    df2.to_json(os.path.join(settings.FILES_DIR, str(path) + ".json"))

    return render(request, 'chemo/quimiotecaDatabase.html')


def readjson(path='GREENIDGE_APAGAR'):
    fullPath = os.path.join(settings.FILES_DIR, str(path) + ".json")
    data = json.load(fullPath)

    return data'''


def feedDatabase(request, path='GREENIDGE_APAGAR', separator='\t'):
    '''
        Feeds the postgreesql database with data from provided path as input.

        VARIABLES:
            - path [string]: Path to file to be read with pandas (format of file should be csv like).
            - separator [string]: Key to perform the split in pandas library.
    '''

    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf")) # Create images for all elements

    df = pd.read_csv(
            os.path.join(settings.FILES_DIR,
            str(path) + ".txt"),
            dtype=str,
            sep=separator
        ) # Parse all elements into a dataframe with pandas
    df = df.drop(
            ['Structure [idcode]', 'Unnamed: 18'],
            axis=1
        ) # Remove useless columns

    for index, row in df.iterrows(): # For each element in dataframe
        try: # Try to parse (this avoid those badly filled entries e.g.: half empties)
            obj, created = Compounds.objects.get_or_create(
                moleculeName=row['Molecule Name'].strip(),
                totalMolweight=float(row['Total Molweight']),
                cLogP=float(row['cLogP']),
                cLogS=float(row['cLogS']),
                hAcceptors=int(row['H-Acceptors']),
                hDonors=int(row['H-Donors']),
                totalSurfaceArea=float(row['Total Surface Area']),
                polarSurfaceArea=float(row['Polar Surface Area']),
                mutagenic=row['Mutagenic'].strip(),
                tumorigenic=row['Tumorigenic'].strip(),
                irritant=row['Irritant'].strip(),
                nonHAtoms=int(row['Non-H Atoms']),
                stereoCenters=int(row['Stereo Centers']),
                rotatableBonds=int(row['Rotatable Bonds']),
                smiles=row['Smiles'].strip(),
                inChI=row['InChI'].strip(),
                inChIKey=row['InChI-Key'].strip(),
            ) # Try to fetch elements, if fails insert elements into database (make database unique)
        except:
            logger.warn("The molecule " + row['Molecule Name'].strip() +
                        " is experiencig some problems, skipping it.") # Show a warning

# region otherdatabase (to fill the other models) [UNFINISHED]
    '''database = ['Molecule Name', 'Total Molweight', 'cLogP', 'cLogS', 'H-Acceptors', 'H-Donors', 'Total Surface Area',
              'Polar Surface Area', 'Mutagenic', 'Tumorigenic', 'Irritant', 'Non-H Atoms', 'Stereo Centers',
              'Rotatable Bonds', 'Smiles', 'InChI', 'InChI-Key']
    
    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf"))
    
    df = pd.read_table(os.path.join(settings.FILES_DIR, str(path) + ".txt"), sep=separator)
    df = df.drop(['Structure [idcode]', 'Unnamed: 18'], axis=1)
    
    for index, row in df.iterrows():
        #for idx in index:
        #    if idx in database:
        objMol, createdMol = Molecules.objects.get_or_create(
            name=row['Molecule Name'].strip(),
        )

        objProp, createdProp = Properties.objects.get_or_create(
            pka = 0.0,#row[''],
            charge = 0.0,#row[''],
            molarMass = 0.0,#row[''],
            monoIsotropicMass = 0.0,#$row[''],
            clogp = row['cLogP'],
            tpsa = 0.0,#row[''],
            lipinski = 0,#row[''],
            hBond = row['H-Acceptors'],
            hBondDonnors = row['H-Donors'],
            rotatableBonds = row['Rotatable Bonds'],
            moleculeID = row['Molecule Name'].strip(),
        )
        
        objStru, createdStru = Structures.objects.get_or_create(
                smiles = row['Smiles'].strip(),
                inChi = row['InChI'].strip(),
                inChiKey = row['InChI-Key'].strip(),
        )

        objFig, createdFig = Figures.objects.get_or_create(
            molecule = 'xxx',#row[''].strip(),
            surface = 'xxx',#row[''].strip(),
            pka = 'xxx',#row[''],
        )'''
# endregion
    # render(request, 'chemo/quimiotecaDatabase.html')
    return redirect('/quimiotecaDatabase')


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
