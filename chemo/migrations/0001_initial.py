# Generated by Django 2.2 on 2019-05-30 20:47

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Compounds',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('moleculeName', models.CharField(max_length=10)),
                ('totalMolweight', models.FloatField()),
                ('cLogP', models.FloatField()),
                ('cLogS', models.FloatField()),
                ('hAcceptors', models.PositiveIntegerField()),
                ('hDonors', models.PositiveIntegerField()),
                ('totalSurfaceArea', models.FloatField()),
                ('polarSurfaceArea', models.FloatField()),
                ('mutagenic', models.CharField(max_length=10)),
                ('tumorigenic', models.CharField(max_length=10)),
                ('irritant', models.CharField(max_length=10)),
                ('nonHAtoms', models.PositiveIntegerField()),
                ('stereoCenters', models.PositiveIntegerField()),
                ('rotatableBonds', models.PositiveIntegerField()),
                ('smiles', models.CharField(max_length=5000)),
                ('inChI', models.CharField(max_length=5000)),
                ('inChIKey', models.CharField(max_length=1000)),
            ],
        ),
        migrations.CreateModel(
            name='Countries',
            fields=[
                ('countryID', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=60)),
                ('continentName', models.CharField(max_length=20)),
            ],
            options={
                'ordering': ('countryID', 'continentName', 'name'),
            },
        ),
        migrations.CreateModel(
            name='Merchants',
            fields=[
                ('merchantID', models.AutoField(primary_key=True, serialize=False)),
                ('merchantName', models.CharField(max_length=1000)),
                ('merchantCreatedAt', models.DateField()),
                ('countryID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Countries')),
            ],
            options={
                'ordering': ('merchantID', 'merchantName', 'countryID', 'merchantCreatedAt', 'userID'),
            },
        ),
        migrations.CreateModel(
            name='Molecules',
            fields=[
                ('moleculeID', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=250)),
            ],
            options={
                'ordering': ('moleculeID', 'name'),
            },
        ),
        migrations.CreateModel(
            name='Users',
            fields=[
                ('userID', models.AutoField(primary_key=True, serialize=False)),
                ('fullName', models.CharField(max_length=400)),
                ('email', models.CharField(max_length=256, unique=True)),
                ('gender', models.CharField(choices=[('Male', 'Male'), ('Female', 'Female'), ('Other', 'Other'), ('Not specified', 'Not specified')], max_length=15)),
                ('birthDate', models.DateField()),
                ('userCreatedAt', models.DateField()),
                ('countryID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Countries')),
            ],
            options={
                'ordering': ('userID', 'fullName', 'email', 'gender', 'birthDate', 'userCreatedAt', 'countryID'),
            },
        ),
        migrations.CreateModel(
            name='Structures',
            fields=[
                ('structureID', models.AutoField(primary_key=True, serialize=False)),
                ('smiles', models.CharField(max_length=5000)),
                ('inChi', models.CharField(max_length=5000)),
                ('inChiKey', models.CharField(max_length=14)),
                ('moleculeID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Molecules')),
            ],
            options={
                'ordering': ('structureID', 'moleculeID', 'smiles', 'inChi', 'inChiKey'),
            },
        ),
        migrations.CreateModel(
            name='Properties',
            fields=[
                ('propertyID', models.AutoField(primary_key=True, serialize=False)),
                ('pka', models.FloatField()),
                ('charge', models.FloatField()),
                ('molarMass', models.FloatField()),
                ('monoIsotropicMass', models.FloatField()),
                ('clogp', models.FloatField()),
                ('tpsa', models.FloatField()),
                ('lipinski', models.SmallIntegerField()),
                ('hBond', models.SmallIntegerField()),
                ('hBondDonnors', models.SmallIntegerField()),
                ('rotatableBonds', models.SmallIntegerField()),
                ('moleculeID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Molecules')),
            ],
            options={
                'ordering': ('propertyID', 'moleculeID', 'pka', 'charge', 'molarMass', 'monoIsotropicMass', 'clogp', 'tpsa', 'lipinski', 'hBond', 'hBondDonnors', 'rotatableBonds'),
            },
        ),
        migrations.CreateModel(
            name='Products',
            fields=[
                ('productID', models.AutoField(primary_key=True, serialize=False)),
                ('name', models.CharField(max_length=250)),
                ('price', models.FloatField()),
                ('status', models.CharField(max_length=100)),
                ('productCreatedAt', models.DateField()),
                ('merchantID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Merchants')),
            ],
            options={
                'ordering': ('productID', 'merchantID', 'price', 'price', 'productCreatedAt', 'name'),
            },
        ),
        migrations.CreateModel(
            name='Orders',
            fields=[
                ('orderID', models.AutoField(primary_key=True, serialize=False)),
                ('status', models.CharField(max_length=100)),
                ('orderCreatedAt', models.DateField()),
                ('userID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Users')),
            ],
            options={
                'ordering': ('orderID', 'userID', 'status', 'orderCreatedAt'),
            },
        ),
        migrations.CreateModel(
            name='OrderItems',
            fields=[
                ('orderItemID', models.AutoField(primary_key=True, serialize=False)),
                ('quantitity', models.PositiveIntegerField()),
                ('productID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Users')),
            ],
            options={
                'ordering': ('orderItemID', 'productID', 'quantitity'),
            },
        ),
        migrations.AddField(
            model_name='merchants',
            name='userID',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Users'),
        ),
        migrations.CreateModel(
            name='Figures',
            fields=[
                ('figureID', models.AutoField(primary_key=True, serialize=False)),
                ('molecule', models.CharField(max_length=5000)),
                ('surface', models.CharField(max_length=5000)),
                ('pka', models.FloatField()),
                ('moleculeID', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemo.Molecules')),
            ],
            options={
                'ordering': ('figureID', 'moleculeID', 'molecule', 'surface', 'pka'),
            },
        ),
    ]
