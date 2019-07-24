import re

from django import forms

from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _

### CUSTOM VALIDATORS ###

def validate_path(value):
    if not re.match(r'^[a-zA-Z0-9+_.]*$', value):
        raise forms.ValidationError('Invalid value')

class HistForm(forms.Form):
    mincutoff = forms.FloatField(label='min cutoff', required=False)
    maxcutoff = forms.FloatField(label='max cutoff', required=False)

class CheckForm(forms.Form):
    choices = forms.MultipleChoiceField(
        widget = forms.CheckboxSelectMultiple(),
        choices = list(map(list, zip(*[list(range(0, 100000)),list(range(0, 100000))])))
    )

class TwodmapFormLine(forms.Form):
    types2d = forms.ChoiceField(
        choices = [("simple", "Simple plot"), ("runningAvg", "Running Average plot")],
        initial = 'simple',
        widget = forms.RadioSelect(),
        required = False
    )
    mincutoff2d = forms.FloatField(label='min cutoff: ', required=False, initial = 0)
    maxcutoff2d = forms.FloatField(label='max cutoff: ', required=False)

class TwodmapFormHeat(forms.Form):
    mincutoffheat2d = forms.FloatField(label='min cutoff: ', required=False, initial = 0)
    maxcutoffheat2d = forms.FloatField(label='max cutoff: ', required=False)

class TwodmapFormDistrib(forms.Form):
    types2ddistrib = forms.ChoiceField(
        choices = [("strip", "Strip plot"), ("box", "Box plot"), ("violin", "Violin plot")],
        initial = 'strip',
        widget = forms.RadioSelect(),
        required = False
    )
    mincutoffdistrib2d = forms.FloatField(label='min cutoff: ', required=False, initial = 0)
    maxcutoffdistrib2d = forms.FloatField(label='max cutoff: ', required=False)

class TwodmapFormFacet(forms.Form):
    types2dfacet = forms.ChoiceField(
        choices = [("fg", "Facet Grids"), ("fgrolling", "Facet Grids Rolling"), ("fgdistplot", "Facet Grid Dist Plot"), ("fgseparate", "Facet Grids Separate")],
        initial = 'fg',
        widget = forms.RadioSelect(),
        required = False
    )
    mincutofffacet2d = forms.FloatField(label='min cutoff: ', required=False, initial = 0)
    maxcutofffacet2d = forms.FloatField(label='max cutoff: ', required=False)

class FilesSubfiles(forms.Form):
    # choices = forms.MultipleChoiceField(
    #     widget  = forms.CheckboxSelectMultiple,
    #     validators=[validate_path]
    # )
    choices = forms.CharField(
        widget=forms.SelectMultiple
    )
