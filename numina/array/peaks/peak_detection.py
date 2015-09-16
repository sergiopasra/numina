import numpy

def find_peaks_index(data, window=5, threshold=0.0):
    from scipy.ndimage.filters import generic_filter
    from ._kernels import kernel_peak_function

    kernel_peak2 = kernel_peak_function(threshold)

    out = generic_filter(data, kernel_peak2, window)
    result, =  numpy.nonzero(out)
    return result


# Auxiliary kernels
def generate_kernel(window):
    from numpy.linalg import inv

    xm = numpy.linspace(-1, 1, window)
    xv = numpy.vander(xm, N=3, increasing=True)
    return numpy.dot(inv(numpy.dot(xv.T, xv)), xv.T)


def refine_peaks(data, ipeaks, window):

    step = window // 2

    winoff = numpy.arange(-step, step+1)
    peakwin = ipeaks[:, numpy.newaxis] + winoff

    ycols = data[peakwin]
    ww = generate_kernel(window)

    coff2 = numpy.dot(ww, ycols.T)

    uc = -0.5 * coff2[1] / coff2[2]
    yc = coff2[0] + uc * (coff2[1] + coff2[2] * uc)
    xc = ipeaks + 0.5 * (window-1) * uc

    return xc, yc
