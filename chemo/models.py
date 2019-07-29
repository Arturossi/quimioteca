from django.db import models

# Create your models here.

class Countries(models.Model):
    countryID = models.AutoField(primary_key=True)
    name = models.CharField(max_length=60) # Largest country name should be "The United Kingdom of Great Britain and Northern Ireland" (56 chars), 60 just in case
    continentName = models.CharField(max_length=20) # Largest continent nome should be either North or South America (13 chars), 20 just in case

    # Metadata
    class Meta:
        ordering = ('countryID', 'continentName', 'name' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Users(models.Model):
    GENDERS = (
        ("Male", "Male"),
        ("Female", "Female"),
        ("Other", "Other"),
        ("Not specified", "Not specified"),
    )

    userID        = models.AutoField(primary_key=True)
    fullName      = models.CharField(max_length=400)
    email         = models.CharField(max_length=256, unique=True) # Biggest email should be 254 (I've put +2 just to be power of 2)
    gender        = models.CharField(max_length=15, choices=GENDERS)
    birthDate     = models.DateField()
    userCreatedAt = models.DateField()
    countryID     = models.ForeignKey(Countries, on_delete=models.CASCADE)

    # Metadata
    class Meta:
        ordering = ('userID', 'fullName', 'email', 'gender', 'birthDate', 'userCreatedAt', 'countryID' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Orders(models.Model):
    orderID        = models.AutoField(primary_key=True)
    userID         = models.ForeignKey(Users, on_delete=models.CASCADE)
    status         = models.CharField(max_length=100)
    orderCreatedAt = models.DateField()

    # Metadata
    class Meta:
        ordering = ('orderID', 'userID', 'status', 'orderCreatedAt' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class OrderItems(models.Model):
    orderItemID = models.AutoField(primary_key=True)
    productID   = models.ForeignKey(Users, on_delete=models.CASCADE)
    quantitity  = models.PositiveIntegerField()

    # Metadata
    class Meta:
        ordering = ('orderItemID', 'productID', 'quantitity' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Merchants(models.Model):
    merchantID        = models.AutoField(primary_key=True)
    merchantName      = models.CharField(max_length=1000)
    countryID         = models.ForeignKey(Countries, on_delete=models.CASCADE)
    merchantCreatedAt = models.DateField()
    userID            = models.ForeignKey(Users, on_delete=models.CASCADE)

    # Metadata
    class Meta:
        ordering = ('merchantID', 'merchantName', 'countryID', 'merchantCreatedAt', 'userID' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Products(models.Model):
    productID        = models.AutoField(primary_key=True)
    merchantID       = models.ForeignKey(Merchants, on_delete=models.CASCADE)
    name             = models.CharField(max_length=250)
    price            = models.FloatField()
    status           = models.CharField(max_length=100)
    productCreatedAt = models.DateField()

    # Metadata
    class Meta:
        ordering = ('productID', 'merchantID', 'price', 'price', 'productCreatedAt', 'name' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Compounds(models.Model):
    moleculeName      = models.CharField(max_length=100) # 2i4x
    totalMolweight    = models.FloatField() # 728.794
    cLogP             = models.FloatField() # 2.8511
    cLogS             = models.FloatField() # -5.343
    hAcceptors        = models.PositiveIntegerField() # 14
    hDonors           = models.PositiveIntegerField() # 2
    totalSurfaceArea  = models.FloatField() # 534.30
    polarSurfaceArea  = models.FloatField() # 186.58
    mutagenic         = models.CharField(max_length=10) # none
    tumorigenic       = models.CharField(max_length=10) # none
    irritant          = models.CharField(max_length=10) # none
    nonHAtoms         = models.PositiveIntegerField() # 49
    stereoCenters     = models.PositiveIntegerField() # 5
    rotatableBonds    = models.PositiveIntegerField() # 19
    smiles            = models.CharField(max_length=5000) # CCOP(COc1ccc(C[C@@H]([C@@H](CN(CC(C)C)S(c(cc2)ccc2OC)(=O)=O)O)NC(O[C@@H]2[C@H](CCO3)[C@H]3OC2)=O)cc1)(OCC)=O
    inChI             = models.CharField(max_length=5000) # InChI=1S/C33H49N2O12PS/c1-6-45-48(38,46-7-2)22-44-26-10-8-24(9-11-26)18-29(34-33(37)47-31-21-43-32-28(31)16-17-42-32)30(36)20-35(19-23(3)4)49(39,40)27-14-12-25(41-5)13-15-27/h8-15,23,28-32,36H,6-7,16-22H2,1-5H3,(H,34,37)/t28-,29-,30+,31-,32-/m0/s1
    inChIKey          = models.CharField(max_length=1000) # FCLYPCIMVVLLRN-LUKCZKMGSA-N
    provider          = models.CharField(max_length=100, default='unkown') # nome do lab
    numAtoms          = models.IntegerField() # 50

class Molecules(models.Model):
    moleculeID = models.AutoField(primary_key=True)
    name = models.CharField(max_length=250)

    # Metadata
    class Meta:
        ordering = ('moleculeID', 'name') 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Properties(models.Model):
    propertyID          = models.AutoField(primary_key=True)
    moleculeID          = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    pka                 = models.FloatField()
    charge              = models.FloatField()
    molarMass           = models.FloatField()
    monoIsotropicMass   = models.FloatField()
    clogp               = models.FloatField()
    tpsa                = models.FloatField()
    lipinski            = models.SmallIntegerField()
    hBond               = models.SmallIntegerField()
    hBondDonnors        = models.SmallIntegerField()
    rotatableBonds      = models.SmallIntegerField()

    # Metadata
    class Meta:
        ordering = ('propertyID', 'moleculeID', 'pka', 'charge', 'molarMass', 'monoIsotropicMass', 'clogp', 'tpsa', 'lipinski', 'hBond', 'hBondDonnors', 'rotatableBonds' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Structures(models.Model):
    structureID = models.AutoField(primary_key=True)
    moleculeID  = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    smiles      = models.CharField(max_length=5000)
    inChi       = models.CharField(max_length=5000)
    inChiKey    = models.CharField(max_length=14)

    # Metadata
    class Meta:
        ordering = ('structureID', 'moleculeID', 'smiles', 'inChi', 'inChiKey' )
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Figures(models.Model):
    figureID    = models.AutoField(primary_key=True)
    moleculeID  = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    molecule    = models.CharField(max_length=5000)
    surface     = models.CharField(max_length=5000)
    pka         = models.FloatField()

    # Metadata
    class Meta:
        ordering = ('figureID', 'moleculeID', 'molecule', 'surface', 'pka' )
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())
