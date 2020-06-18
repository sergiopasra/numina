
from ..datamodel import QueryAttribute, KeyDefinition

import numina.util.convert as conv
from ..datamodel import get_imgid

_gtc_mappings = {
    'instrument': 'INSTRUME',
    'object': 'OBJECT',
    'observation_date': ('DATE-OBS', 0, conv.convert_date),
    'uuid': 'uuid',
    'type': 'numtype',
    'mode': 'obsmode',
    'exptime': 'EXPTIME',
    'darktime': 'DARKTIME',
    'quality_control': ('NUMRQC', 0, conv.convert_qc),
    'insmode': ('INSMODE', 'undefined'),
    'imgid': get_imgid,
    'insconf': 'INSCONF',
    'blckuuid': lambda x: '1',
    'insconf_uuid': 'insconf', # Alias
    'block_uuid': 'blckuuid', # Alias
}


class FITSKeyExtractor(object):
    """Extract values from FITS images"""
    def __init__(self, values):
        self.map = {}
        for key, entry in values.items():
            if isinstance(entry, KeyDefinition):
                newval = entry
            elif isinstance(entry, tuple):
                if len(entry) == 3:
                    keyname = entry[0]
                    hduname = entry[1]
                    convert = entry[2]
                    default = None
                elif len(entry) == 2:
                    keyname = entry[0]
                    default = entry[1]
                    hduname = 0
                    convert = None
                else:
                    raise ValueError

                newval = KeyDefinition(
                    keyname,
                    ext=hduname,
                    convert=convert,
                    default=default
                )
            elif isinstance(entry, str):
                newval = KeyDefinition(
                    entry
                )
            else:
                newval = entry

            self.map[key] = newval

    def extract(self, hdulist, value):
        extractor = self.map[value]
        return extractor(hdulist)

def gather_info(dframe):
    """Obtain a summary of information about the image."""
    with dframe.open() as hdulist:
        info = gather_info_hdu(hdulist)
    return info


def gather_info_dframe(self, img):
    """Obtain a summary of information about the image."""
    return gather_info(img)


_meta_dinfo_headers = [
        'instrument',
        'object',
        'observation_date',
        'uuid',
        'type',
        'mode',
        'exptime',
        'darktime',
        'insconf',
        'blckuuid',
        'quality_control',
        'block_uuid',  # Alias
        'insconf_uuid',  # Alias
        'imgid',
    ]


def gather_info_hdu(hdulist):
    """Obtain a summary of information about the image."""
    values = {}
    values['n_ext'] = len(hdulist)
    extnames = [hdu.header.get('extname', '') for hdu in hdulist[1:]]
    values['name_ext'] = ['PRIMARY'] + extnames

    fits_extractor = FITSKeyExtractor(values)
    for key in _meta_dinfo_headers:
        values[key] = fits_extractor.extract(hdulist, key)

    return values
