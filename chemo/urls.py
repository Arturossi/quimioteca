from django.urls import path

from . import views

app_name = 'chemo'

# Patterns which leads to a page
urlpatterns = [
    path('', views.IndexView.as_view(), name='home'),
    path('index', views.IndexView.as_view(), name='index'),
    path('compound', views.compoundView.as_view(), name='compound'),
    path('about', views.aboutView.as_view(), name='abnout'),
    path('contactUs', views.contactUsView.as_view(), name='contactUs'),
    path('login', views.loginView.as_view(), name='login'),
    path('quimiotecaDatabase', views.quimiotecaDatabaseView.as_view(),
         name='quimiotecaDatabase'),
    path('cadastroMoleculas', views.cadastroMolView.as_view(), name='cadastroMoleculas'),
    path('sucesso', views.sucessoView.as_view(), name='sucesso'),
]

noviewsurls = [
    path('updateDatabase', views.feedDatabase, name='updateDatabase'),
    path('uploadMolecule', views.uploadMolecule, name='uploadMolecule'),
]

urlpatterns = urlpatterns + noviewsurls
