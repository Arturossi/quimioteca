from django.contrib import admin

from .models import *

# Register your models here.

admin.site.register(Countries)
admin.site.register(Users)
admin.site.register(Orders)
admin.site.register(OrderItems)
admin.site.register(Merchants)
admin.site.register(Products)

admin.site.register(Molecules)
admin.site.Properties(Products)
admin.site.Structures(Products)
admin.site.Figures(Products)