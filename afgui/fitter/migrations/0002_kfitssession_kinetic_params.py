# -*- coding: utf-8 -*-
# Generated by Django 1.10.1 on 2017-08-02 13:27
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('fitter', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='kfitssession',
            name='kinetic_params',
            field=models.TextField(null=True),
        ),
    ]
