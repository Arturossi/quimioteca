from django.urls import path

from .views import *

app_name = 'chemo'

# Patterns which leads to a page
urlpatterns = [
    path('', IndexView.as_view(), name='home'),
    path('index', IndexView.as_view(), name='index'),
    path('compound', compoundView.as_view(), name='compound'),
    path('about', aboutView.as_view(), name='abnout'),
    path('contactUs', contactUsView.as_view(), name='contactUs'),
    path('login', loginView.as_view(), name='login'),
    path('quimiotecaDatabase', quimiotecaDatabaseView.as_view(),
         name='quimiotecaDatabase'),
    path('cadastroMoleculas', cadastroMolView.as_view(), name='cadastroMoleculas'),
    path('sucesso', sucessoView.as_view(), name='sucesso'),
]

noviewsurls = [
    path('updateDatabase', feedDatabase, name='updateDatabase'),
    path('uploadMolecule', uploadMolecule, name='uploadMolecule'),
]

urlpatterns = urlpatterns + noviewsurls
