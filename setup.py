from setuptools import setup

setup(name='psquery',
      version='0.6.0',
      url='http://github.com/realfastvla/psquery',
      packages=['psquery'],
      install_requies=[
          'astropy',
          'extinction',
          'sedpy',
          'dustmaps',
          'mastcasjobs',
          'astroquery',
          'pyvo'
      ],
      zip_safe=False)
