# psquery
Query catalogs with MAST API. Basic query is a cone search with arguments of (RA, Dec, radius). 

Currently, set ups as two modules for PS1 and 2MASS queries.

Returns either astropy coordinate objects or PS1 query result.

See notebook for example execution.

## Dependencies
- astropy
- urllib
- requests
- pandas

## Installation

`python setup.py install`
