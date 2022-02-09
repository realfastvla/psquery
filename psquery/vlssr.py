from astroquery import heasarc
from astropy import coordinates, units
import click

heq = heasarc.Heasarc()

@click.command()
@click.argument('source')
@click.option('--radius', default=0.5, help='radius in degrees')
def vlssr(source, radius):
    """ Source can be an radec string (no comma, hmsdms or deg) or source name as string.
    """

    try:
        tab = heq.query_region(source, mission='VLSSr', radius=0.5*units.deg)
        print(tab)
    except TypeError:
        print("No sources found")

if __name__ == '__main__':
    vlssr()
