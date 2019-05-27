# Generated by Django 2.2 on 2019-05-27 19:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemo', '0001_initial'),
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
    ]