from django.shortcuts import render, redirect
from django.views import generic
from django.conf import settings

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
import logging

import pandas as pd

# Logger to log
logger = logging.getLogger(__name__)

# Create your views here.

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

#a - Desenho (PNG)
#b - <Molecule Name>
#c - <Total Molweight>
#d - <cLogP>
#e - <cLogS>
#f - <Polar Surface Area>

def generateImages(path):
    sdf = Chem.SDMolSupplier(path)#'GREENIDGE_APAGAR.sdf')
    ms = [x for x in sdf if x is not None]

    for m in ms:
        tmp=AllChem.Compute2DCoords(m)

    for i in range(len(ms)) :
        name=ms[i].GetProp('_Name')
        Draw.MolToFile(ms[i], os.path.join(settings.FILES_DIR, f'molImages/{name}.png'))
    
    return

def csv2json(request, path='GREENIDGE_APAGAR', separator='\t'):
    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf"))
    
    df = pd.read_csv(os.path.join(settings.FILES_DIR, str(path) + ".txt"), sep=separator)
    df2 = df.drop(['Structure [idcode]', 'Unnamed: 18'], axis=1)
    df2.to_json(os.path.join(settings.FILES_DIR, str(path) + ".json"))

    return render(request, 'chemo/quimiotecaDatabase.html')

def readjson(path='GREENIDGE_APAGAR'):
    fullPath = os.path.join(settings.FILES_DIR, str(path) + ".json")
    data = json.load(fullPath)

    return data

def feedDatabase(request, path='GREENIDGE_APAGAR', separator='\t'):
    '''
        Feeds the postgreesql database with data from provided path as input.

        VARIABLES:
            - path [string]: Path to file to be read with pandas (format of file should be csv like).
            - separator [string]: Key to perform the split in pandas library.
    '''
    
    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf"))

    df = pd.read_csv(os.path.join(settings.FILES_DIR, str(path) + ".txt"), dtype=str, sep=separator)
    df = df.drop(['Structure [idcode]', 'Unnamed: 18'], axis=1)
    
    for index, row in df.iterrows():
        try:
            obj, created = Compounds.objects.get_or_create(
                moleculeName = row['Molecule Name'].strip(),
                totalMolweight = float(row['Total Molweight']),
                cLogP = float(row['cLogP']),
                cLogS = float(row['cLogS']),
                hAcceptors = int(row['H-Acceptors']),
                hDonors = int(row['H-Donors']),
                totalSurfaceArea = float(row['Total Surface Area']),
                polarSurfaceArea = float(row['Polar Surface Area']),
                mutagenic = row['Mutagenic'].strip(),
                tumorigenic = row['Tumorigenic'].strip(),
                irritant = row['Irritant'].strip(),
                nonHAtoms = int(row['Non-H Atoms']),
                stereoCenters = int(row['Stereo Centers']),
                rotatableBonds = int(row['Rotatable Bonds']),
                smiles = row['Smiles'].strip(),
                inChI = row['InChI'].strip(),
                inChIKey = row['InChI-Key'].strip(),
            )
        except:
            logger.warn("The molecule " + row['Molecule Name'].strip() + " is experiencig some problems, skipping it.")

#region otherdatabase
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
#endregion
    return redirect('/quimiotecaDatabase')#render(request, 'chemo/quimiotecaDatabase.html')

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

        #loadSDF(os.path.join(settings.FILES_DIR, 'ApprovedDrugs2015.sdf'))

        # Render page index.html within the request and variables
        return render(request, 'chemo/index.html', {'countries': updateCountries()})

class compoundView(generic.ListView):
    def get(self, request, **kwargs):
        return render(request, 'chemo/compound.html')

class aboutView(generic.ListView):
    def get(self, request, **kwargs):
        return render(request, 'chemo/about.html')

class contactUsView(generic.ListView):
    def get(self, request, **kwargs):
        return render(request, 'chemo/contactUs.html')

class loginView(generic.ListView):
    def get(self, request, **kwargs):
        return render(request, 'chemo/login.html')

class quimiotecaDatabaseView(generic.ListView):
    def get(self, request, **kwargs):
        mols = Compounds.objects.all()

        variables = {
            'mols': mols,
        }

        #return render(request, 'chemo/quimiotecaDatabase.html', variables)
        return render(request, 'chemo/quimiotecaDatabase.html')