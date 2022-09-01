# psquery
Query astronomical catalogs with Python APIs or libraries. Also some support for optical SED modeling.

The best supported use case is a cone search in radio and optical catalogs. See `psquery_cone_examples.ipynb` for sample of cone search queries.

## Dependencies
- astropy
- astroquery
- urllib
- requests
- pandas

Optional (for specific data sets):
- pyvo
- noaodatalab
- cfod
- mastcasjobs: `pip install git+https://github.com/rlwastro/mastcasjobs@master`

Optional (for SED modeling):
- sedpy
- extinction
- dustmaps

## Installation

`python setup.py install`

mastcasjobs may require some hand holding. sedpy requires lots of hand holding (only needed for optical SED modeling).