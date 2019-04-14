from django.urls import path

from . import views

app_name = 'chemo'

# Patterns which leads to a page
urlpatterns = [
    path('', views.IndexView.as_view(), name='test'),
]