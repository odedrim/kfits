from django import forms

class UploadFileForm(forms.Form):
    fdata = forms.FileField()
