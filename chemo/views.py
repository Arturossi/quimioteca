from django.shortcuts import render
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
import pandas as pd

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

def feedDatabase(path):
    database = ['Molecule Name', 'Total Molweight', 'cLogP', 'cLogS', 'H-Acceptors', 'H-Donors', 'Total Surface Area',
              'Polar Surface', 'Area', 'Mutagenic', 'Tumorigenic', 'Irritant', 'Non-H Atoms', 'Stereo Centers',
              'Rotatable Bonds', 'Smiles', 'InChI', 'InChI-Key']
    
    generateImages(os.path.join(settings.FILES_DIR, 'GREENIDGE_APAGAR.sdf'))
    
    df = pd.read_table(os.path.join(settings.FILES_DIR, 'GREENIDGE_APAGAR.txt'), sep='\t')
    
    for index, row in df.iterrows():
        #for idx in index:
        #    if idx in database:
        objMol, createdMol = Molecules.objects.get_or_create(
            name=row['Molecule Name'].strip(),
        )

        objProp, createdProp = Properties.objects.get_or_create(
            pka = row[''].strip(),
            charge = row[''].strip(),
            molarMass = row[''].strip(),
            monoIsotropicMass = row[''].strip(),
            clogp = row['cLogP'].strip(),
            tpsa = row[''].strip(),
            lipinski = row[''].strip(),
            hBond = row['H-Acceptors'].strip(),
            hBondDonnors = row['H-Donors'].strip(),
            rotatableBonds = row['Rotatable Bonds'].strip(),
        )

        objStru, createdStru = Structures.objects.get_or_create(
                smiles = row['Smiles'].strip(),
                inChi = row['InChI'].strip(),
                inChiKey = row['InChI-Key'].strip(),
        )

        objFig, createdFig = Figures.objects.get_or_create(
            molecule = row[''].strip(),
            surface = row[''].strip(),
            pka = row[''].strip(),
        )

    return

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