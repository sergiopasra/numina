#
# Copyright 2018 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""User command line interface for Montecarlo simulations."""

import os

from .clirundal import mode_run_common


def complete_config(config):
    """Complete config with default values"""

    if not config.has_section('run'):
        config.add_section('run')

    values = {
        'basedir': os.getcwd(),
        'task_control': 'control.yaml',
    }

    for k, v in values.items():
        if not config.has_option('run', k):
            config.set('run', k, v)

    return config

def register(subparsers, config):

    complete_config(config)

    task_control_base = config.get('run', 'task_control')

    parser_run = subparsers.add_parser(
        'mc',
        help='run MC simulations'
        )

    parser_run.set_defaults(command=mode_run_mc)

    parser_run.add_argument(
        '-r', '--repeat', dest='repeat', default=10,
        help='repeat times',
        )

    parser_run.add_argument(
        'file', help='base file for simulation',
        )

    return parser_run


def mode_run_mc(args, extra_args):
    import numina.drps

    print('MC')
    print(args.file)

    import astropy.io.fits as fits
    with fits.open(args.file) as ffile:
        rv_fits(ffile)


def rv_fits(fitsfile):
    import numina.drps
    import astropy.io.fits as fits

    drpsys = numina.drps.get_system_drps()

    instrument = fitsfile[0].header['INSTRUME']
    dm = drpsys.drps[instrument].datamodel

    if hasattr(dm, 'rvs'):
        for idx, s in enumerate(dm.rvs(fitsfile, size=1)):
            print(s)
            s.writeto('/tmp/dum-{}.fits'.format(idx), overwrite=True)
            diff = s[0].data.astype('float') - fitsfile[0].data.astype('float')
            fits.writeto('/tmp/diff-{}.fits'.format(idx), diff, overwrite=True)
    else:
        print('datamodel does not have random variate sampling')
