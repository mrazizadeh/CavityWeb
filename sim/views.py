from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponse
from .models import InputForm
from .compute import compute_CFD
import os

def index(request):
    os.chdir(os.path.dirname(__file__))
    result = None
    if request.method == 'POST':
        form = InputForm(request.POST)
        if form.is_valid():
            form2 = form.save(commit=False)
            result = compute_CFD (form2.A, form2.B, form2.C, form2.D, form2.E, form2.F)
            result = result.replace('static/', '')
    else:
        form = InputForm()
    context =  {'form': form,
             'result': result,
             }
    return render(request,'sim.html', context
            )
# Create your views here.
'''
'''