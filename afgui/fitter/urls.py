from django.conf import urls
from . import views

urlpatterns = [
    urls.url(r'^$', views.index, name='index'),
    urls.url(r'^index.html?$', views.index, name='index'),
    urls.url(r'^test.html?$', views.test, name='test'),
    urls.url(r'^bootstrap/.*$', views.bootstrap, name='bootstrap'),
    urls.url(r'^backend$', views.backend, name='fbackend'),
    urls.url(r'^upload_text.html?$', views.upload_text, name='upload_text'),
]
