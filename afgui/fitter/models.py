from __future__ import unicode_literals

from django.db import models

# Create your models here.

class KfitsSession(models.Model):
    sid = models.CharField(max_length=32, primary_key=True)
    timestamp = models.DateTimeField(auto_now=True)
    saved_output = models.TextField(null=True)
    kinetic_params = models.TextField(null=True)
