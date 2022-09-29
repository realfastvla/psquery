from setuptools import setup

setup(name='psquery-astro',
      version='0.6.8',
      url='http://github.com/realfastvla/psquery',
      packages=['psquery'],
      author="Casey Law",
      author_email="caseyjlaw@gmail.com",
      install_requires=[
          'astropy',
          'astroquery',
          'pandas',
          'extinction',
          'sedpy',
          'dustmaps',
          'mastcasjobs',
          'pyvo',
          'noaodatalab',
          'cfod'
      ],
      zip_safe=False)
