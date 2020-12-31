from django.db import models
from django.forms import ModelForm
from math import pi

class Input(models.Model):
    A = models.IntegerField(
        verbose_name=' maxstep ', default= 4000)
    B = models.IntegerField(
        verbose_name=' grid ', default= 71.0)
    C = models.FloatField(
        verbose_name=' Length ', default= 1.0)
    D = models.FloatField(
        verbose_name=' Uwall', default= 2.0)
    E = models.FloatField(
        verbose_name=' nu ', default= 0.05)
    F = models.FloatField(
        verbose_name=' deltat ', default= 0.0003)

class InputForm(ModelForm):
    class Meta:
        model = Input
        fields = '__all__'
# Create your models here.
