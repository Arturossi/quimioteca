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
admin.site.register(Properties)
admin.site.register(Structures)
admin.site.register(Figures)

admin.site.register(Compounds)