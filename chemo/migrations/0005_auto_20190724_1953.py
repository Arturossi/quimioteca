# Generated by Django 2.2.3 on 2019-07-24 19:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemo', '0004_compounds_molname'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='compounds',
            name='molName',
        ),
        migrations.AlterField(
            model_name='compounds',
            name='moleculeName',
            field=models.CharField(max_length=100),
        ),
    ]