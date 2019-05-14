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

    userID = models.AutoField(primary_key=True)
    fullName = models.CharField(max_length=400)
    email = models.CharField(max_length=256, unique=True) # Biggest email should be 254 (I've put +2 just to be power of 2)
    gender = models.CharField(max_length=15, choices=GENDERS)
    birthDate = models.DateField()
    userCreatedAt = models.DateField()
    countryID = models.ForeignKey(Countries, on_delete=models.CASCADE)

    # Metadata
    class Meta:
        ordering = ('userID', 'fullName', 'email', 'gender', 'birthDate', 'userCreatedAt', 'countryID' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Orders(models.Model):
    orderID = models.AutoField(primary_key=True)
    userID = models.ForeignKey(Users, on_delete=models.CASCADE)
    status = models.CharField(max_length=100)
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
    productID = models.ForeignKey(Users, on_delete=models.CASCADE)
    quantitity = models.PositiveIntegerField()

    # Metadata
    class Meta:
        ordering = ('orderItemID', 'productID', 'quantitity' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Merchants(models.Model):
    merchantID = models.AutoField(primary_key=True)
    merchantName = models.CharField(max_length=1000)
    countryID = models.ForeignKey(Countries, on_delete=models.CASCADE)
    merchantCreatedAt = models.DateField()
    userID = models.ForeignKey(Users, on_delete=models.CASCADE)

    # Metadata
    class Meta:
        ordering = ('merchantID', 'merchantName', 'countryID', 'merchantCreatedAt', 'userID' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Products(models.Model):
    productID = models.AutoField(primary_key=True)
    merchantID = models.ForeignKey(Merchants, on_delete=models.CASCADE)
    name = models.CharField(max_length=250)
    price = models.FloatField()
    status = models.CharField(max_length=100)
    productCreatedAt = models.DateField()

    # Metadata
    class Meta:
        ordering = ('productID', 'merchantID', 'price', 'price', 'productCreatedAt', 'name' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

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
    propertyID = models.AutoField(primary_key=True)
    moleculeID = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    pka = models.FloatField()
    charge = models.FloatField()
    molarMass = models.FloatField()
    monoIsotropicMass = models.FloatField()
    clogp = models.FloatField()
    tpsa = models.FloatField()
    lipinski = models.SmallIntegerField()
    hBond = models.SmallIntegerField()
    hBondDonnors = models.SmallIntegerField()
    rotatableBonds = models.SmallIntegerField()

    # Metadata
    class Meta:
        ordering = ('propertyID', 'moleculeID', 'pka', 'charge', 'molarMass', 'monoIsotropicMass', 'clogp', 'tpsa', 'lipinski', 'hBond', 'hBondDonnors', 'rotatableBonds' ) 
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Structures(models.Model):
    structureID = models.AutoField(primary_key=True)
    moleculeID = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=5000)
    inChi = models.CharField(max_length=5000)
    inChiKey = models.CharField(max_length=14)

    # Metadata
    class Meta:
        ordering = ('structureID', 'moleculeID', 'smiles', 'inChi', 'inChiKey' )
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())

class Figures(models.Model):
    figureID = models.AutoField(primary_key=True)
    moleculeID = models.ForeignKey(Molecules, on_delete=models.CASCADE)
    molecule = models.CharField(max_length=5000)
    surface = models.CharField(max_length=5000)
    pka = models.FloatField()

    # Metadata
    class Meta:
        ordering = ('figureID', 'moleculeID', 'molecule', 'surface', 'pka' )
         # helps in alphabetical listing. Sould be a tuple

    def __str__(self):
        # It is not a nice print but it is sort of readable (SORT OF)
        return ', '.join('{}{}'.format(key, val) for key, val in self.__dict__.items())