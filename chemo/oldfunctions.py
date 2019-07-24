### old Functions

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

def csv2json(request, path='GREENIDGE_APAGAR', separator='\t'):
    generateImages(os.path.join(settings.FILES_DIR, str(path) + ".sdf"))

    df = pd.read_csv(os.path.join(settings.FILES_DIR,
                                  str(path) + ".txt"), sep=separator)
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
    #return redirect('/quimiotecaDatabase')
