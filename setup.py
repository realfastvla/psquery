from setuptools import setup

setup(name='psquery-astro',
      version='0.6.1',
      url='http://github.com/realfastvla/psquery',
      packages=['psquery'],
      author="Casey Law",
      author_email="caseyjlaw@gmail.com",
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
