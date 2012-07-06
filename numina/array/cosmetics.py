#
# Copyright 2008-2012 Universidad Complutense de Madrid
# 
# This file is part of Numina
# 
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

import logging

import numpy
import scipy.stats
import scipy.ndimage

from numina.array.blocks import max_blk_coverage, blk_nd_short

# Values stored in integer masks
PIXEL_HOT = 2
PIXEL_DEAD = 1
PIXEL_VALID = 0

#
HIGH_SIGMA = 200
LOW_SIGMA = -200 

_logger = logging.getLogger('numina.cosmetics')

def update_mask(mask, gmask, newmask, value):
    f1_mask = mask[gmask]
    f1_mask[newmask] = value
    mask[gmask] = f1_mask
    gmask = mask == PIXEL_VALID
    smask = mask != PIXEL_VALID
    return mask, gmask, smask

def cosmetics(flat1, flat2, mask, lowercut=4.0, uppercut=4.0, 
                nsig=2.0, posthook=None):
    '''Find cosmetic defects in a detector using two flat field images.
    
    Two arrays representing flat fields of different exposure times are
    required. Cosmetic defects are selected as points that deviate
    significantly of the expected normal distribution of pixels in
    the ratio between `flat2` and `flat1`.
    
    The median of the ratio array is computed and subtracted to it.
    
    The standard deviation of the distribution of pixels is computed
    obtaining the percentiles nearest the pixel values corresponding to
    `nsig` in the normal CDF. The standar deviation is then the distance
    between the pixel values divided by two times `nsig`.
    The ratio image is then normalized with this standard deviation.
    
    The values in the ratio above `uppercut` are flagged as hot pixels,
    and those below '-lowercut` are flagged as dead pixels in the output mask. 
    
    :parameter flat1: an array representing a flat illuminated exposure.
    :parameter flat2: an array representing a flat illuminated exposure.
    :parameter mask: an integer array representing initial mask.
    :parameter lowercut: values bellow this sigma level are flagged as dead pixels.
    :parameter uppercut: values above this sigma level are flagged as hot pixels.
    :parameter nsig: level to estimate the standar deviation.
    :returns: the normalized ratio of the flats, the updated mask and standard deviation
    
    '''

    # check inputs

    for vname in ['lowercut', 'uppercut', 'nsig']:
        val = locals()[vname]
        if val <= 0:
            raise ValueError('%s must be > 0' % vname)

    ratio = numpy.zeros_like(flat1)
    invalid = numpy.zeros_like(flat1)

    gmask = mask == PIXEL_VALID
    _logger.info('valid points in input mask %d', numpy.count_nonzero(gmask))
    smask = mask != PIXEL_VALID
    _logger.info('invalid points in input mask %d', numpy.count_nonzero(smask))

    # check if there are zeros in flat1 and flat2
    zero_mask = numpy.logical_or(flat1[gmask] <= 0, flat2[gmask] <= 0)
    
    # if there is something in zero mask
    # we update the mask
    if numpy.any(zero_mask):
        mask, gmask, smask = update_mask(mask, gmask, zero_mask, PIXEL_DEAD)
        invalid[mask == PIXEL_DEAD] = LOW_SIGMA

    # ratio of flats
    ratio[gmask] = flat2[gmask] / flat1[gmask]
    ratio[smask] = invalid[smask]

    _logger.info('computing median using all the image')
    ratio_med = numpy.median(ratio[gmask])

    # subtracting the median map
    ratio[gmask] -= ratio_med

    # Quantiles that contain nsig sigma in normal distribution
    _logger.debug('estimating sigma using nsig=%f', nsig)
    qns = 100 * scipy.stats.norm.cdf(nsig)
    pns = 100 - qns
    _logger.debug('percentiles at nsig=%f', nsig)
    _logger.debug('low %f%% high %f%%', pns, qns)

    valid = ratio[gmask]
    ls = scipy.stats.scoreatpercentile(valid.flat, pns)
    hs = scipy.stats.scoreatpercentile(valid.flat, qns)

    _logger.debug('score at percentiles')
    _logger.debug('low %f high %f', ls, hs)

    # sigma estimation
    sig = (hs - ls) / (2 * nsig)
    _logger.info('sigma estimation is %f ', sig)

    # normalized points
    ratio_base = ratio[gmask]
    ratio[gmask] /= sig

    f1_ratio = ratio[gmask]
    f1_mask = mask[gmask]
    f1_mask[f1_ratio >= uppercut] = PIXEL_HOT
    f1_mask[f1_ratio <= -lowercut] = PIXEL_DEAD
    mask[gmask] = f1_mask

    if posthook is not None:
        posthook(ratio_base, 0.0, sig)

    return ratio, mask, sig

# IRAF task
def ccdmask(flat1, flat2=None, mask=None, ncmed=7, nlmed=7, 
                ncsig=15, nlsig=15,
                lowercut=6.0, uppercut=6.0, 
                nsig=1.0, posthook=None):
    '''Find cosmetic defects in a detector using two flat field images.
    
    Two arrays representing flat fields of different exposure times are
    required. Cosmetic defects are selected as points that deviate
    significantly of the expected normal distribution of pixels in
    the ratio between `flat2` and `flat1`.
    
    The median of the ratio array is computed in boxes of (nlmed, ncmed) 
    and subtracted to it.
    
    The standard deviation of the distribution of pixels is computed
    in boxes of (nlsig, ncsig) pixels. It is
    obtained computing the percentiles nearest to the pixel values 
    corresponding to`nsig` in the normal CDF. 
    The standard deviation is then the distance
    between the pixel values divided by two times `nsig`.
    The ratio image is then normalized with this standard deviation.
    
    The values in the ratio above `uppercut` are flagged as hot pixels,
    and those below '-lowercut` are flagged as dead pixels in the output mask. 
    
    :parameter flat1: an array representing a flat illuminated exposure.
    :parameter flat2: an array representing a flat illuminated exposure.
    :parameter mask: an integer array representing initial mask.
    :parameter ncmed: number of columns used to compute the median
    :parameter nlmed: number of lines used to compute the median
    :parameter ncsig: number of columns used to estimate the standard deviation
    :parameter nlsig: number of lines used to estimate the standard deviation
    :parameter lowercut: values below this sigma level are flagged as dead pixels.
    :parameter uppercut: values above this sigma level are flagged as hot pixels.
    :parameter nsig: level to estimate the standard deviation.
    :returns: the normalized ratio of the flats, the updated mask and standard deviation
    
    .. note::
    
        This function is based on the description of the task
        ccdmask of IRAF
    
    .. seealso::
    
        :py:func:`cosmetics`
            Operates much like this function but computes
            median and sigma in the hole image instead of in boxes          
    
    
    '''
    # check inputs
    for vname in ['lowercut', 'uppercut', 'nsig', 'ncsig', 'nlsig',
                    'ncmed', 'nlmed']:
        val = locals()[vname]
        if val <= 0:
            raise ValueError('%s must be > 0' % vname)

    if ncsig * nlsig < 100:
        raise ValueError('ncsig * nlsig >= 100 required to have enough points for statistics')

    if flat2 is None:
        # we have to swap flat1 and flat2, and
        # make flat1 an array of 1s
        flat1, flat2 = flat2, flat1
        flat1 = numpy.ones_like(flat2)

    if mask is None:
        mask = numpy.zeros_like(flat1, dtype='int')

    ratio = numpy.zeros_like(flat1)
    invalid = numpy.zeros_like(flat1)
    invalid[mask == PIXEL_HOT] = HIGH_SIGMA
    invalid[mask == PIXEL_DEAD] = LOW_SIGMA

    gmask = mask == PIXEL_VALID
    _logger.info('valid points in input mask %d', numpy.count_nonzero(gmask))
    smask = mask != PIXEL_VALID
    _logger.info('invalid points in input mask %d', numpy.count_nonzero(smask))

    # check if there are zeros in flat1 and flat2
    zero_mask = numpy.logical_or(flat1[gmask] <= 0, flat2[gmask] <= 0)

    # if there is something in zero mask
    # we update the mask
    if numpy.any(zero_mask):
        mask, gmask, smask = update_mask(mask, gmask, zero_mask, PIXEL_DEAD)
        invalid[mask == PIXEL_DEAD] = LOW_SIGMA
  
    # ratio of flats
    ratio[gmask] = flat2[gmask] / flat1[gmask]
    ratio[smask] = invalid[smask]

    _logger.info('computing median in boxes of %r', (nlmed, ncmed))
    ratio_med = scipy.ndimage.filters.median_filter(ratio, size=(nlmed, ncmed))

    # subtracting the median map
    ratio[gmask] -= ratio_med[gmask]

    # Quantiles that contain nsig sigma in normal distribution
    qns = 100 * scipy.stats.norm.cdf(nsig)
    pns = 100 - qns
    _logger.info('percentiles at nsig=%f', nsig)
    _logger.info('low %f%% high %f%%', pns, qns)

    # in several blocks of shape nlsig, ncsig
    # we estimate sigma
    sigma = numpy.zeros_like(ratio)
    medianbox = (nlsig, ncsig)
    mshape = max_blk_coverage(blk=medianbox, shape=ratio.shape)

    _logger.info('estimating sigma in boxes of %r', medianbox)
    _logger.info('shape covered by boxes is  %r', mshape)
    for blk in blk_nd_short(blk=medianbox, shape=ratio.shape):
        # mask for this region
        m = mask[blk] == PIXEL_VALID
        valid_points = numpy.ravel(ratio[blk][m])
        ls = scipy.stats.scoreatpercentile(valid_points, pns)
        hs = scipy.stats.scoreatpercentile(valid_points, qns)

        _logger.debug('score at percentiles')
        _logger.debug('low %f high %f', ls, hs)

        # sigma estimation
        sig = (hs - ls) / (2 * nsig)
        _logger.debug('sigma estimation is %f ', sig)

        # normalized points
        sigma[blk] = sig

    # fill regions of sigma not computed
    fill0 = ratio.shape[0] - mshape[0]    
    fill1 = ratio.shape[1] - mshape[1]    
    if fill0 > 0:
        _logger.info('filling %d rows in sigma image', fill0)
        sigma[:, mshape[0]:] = sigma[:,mshape[0] - fill0: mshape[0]]

    if fill1 > 0:
        _logger.info('filling %d columns in sigma image', fill1)
        sigma[mshape[1]:, :] = sigma[mshape[1] - fill1: mshape[1], :]

    invalid_sigma = sigma <= 0.0
    if numpy.any(invalid_sigma):
        _logger.info('updating mask with points where sigma <=0')
        mask, gmask, smask = update_mask(mask, gmask, invalid_sigma, PIXEL_HOT)
        invalid[mask == PIXEL_HOT] = HIGH_SIGMA

    ratio[gmask] /= sigma[gmask]

    f1_ratio = ratio[gmask]
    f1_mask = mask[gmask]
    f1_mask[f1_ratio >= uppercut] = PIXEL_HOT
    f1_mask[f1_ratio <= -lowercut] = PIXEL_DEAD
    mask[gmask] = f1_mask

    return ratio, mask, sigma
