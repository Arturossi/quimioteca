from django import forms

### CUSTOM VALIDATORS ###

class UploadFileForm(forms.Form):
    title = forms.CharField(max_length=50)
    file = forms.FileField()
