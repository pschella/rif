"""This is a Python module for working with data from the Radio Interferometer.

It has two low level functions `nsamples` and `readdata` to read in the raw data and a series of convenience functions to quickly analyse the data and reduce the total data volume before further processing.

RIF data format
===============

RIF data is stored in a combination of two files:

* A ``.cds`` plain text file with telescope coordinates stored every second and,
* a ``.dat`` binary file with the telescope data stored as ``A0, B0, A1, B1, ...`` where ``A0`` is sample 0 from channel (telescope) A. The data type is an 8 bit (1 byte) integer (or C-type char).

The routines in the ``rif`` module can be used to read and manipulate this data in Python. Alternatively one can use any programming language to read the data directly using the ``rif`` module code (which itself is written in C) as a reference.

API documentation
=================
"""

# Import core functionality
import _rif

# Define which functions to export
__all__ = ['readdata', 'nsamples', 'totalpower', 'dynamicspectrum', 'crosscorrelation']

# Import additional modules
import numpy as np

def readdata(filename, start, n):
    """Read *n* samples from file *filename* starting at sample *start*.

    Returns a two dimensional array of shape=(2, n/2) with the data
    for both channels.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Read 1024 samples for each channel starting at sample 0.
        d = rif.readdata("RIF.dat", 0, 2048)

        # Plot raw time series of:
        plt.plot(d[0]) # channel 0
        plt.plot(d[1]) # channel 1.

        plt.show()

    """
    data = _rif.readdata(filename, start, n)

    return data.reshape((n/2, 2)).transpose()

def readdata_single_channel(filename, start, n):
    """Read *n* samples from file *filename* starting at sample *start*.

    Returns a two dimensional array with the data.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Read 1024 samples for each channel starting at sample 0.
        d = rif.readdata("RIF.dat", 0, 2048)

        # Plot raw time series of:
        plt.plot(d)

        plt.show()

    """
    data = _rif.readdata(filename, start, n)

    return data

def nsamples(filename):
    """Returns the number of samples in file *filename*.

    Example::

        # Import modules
        import rif

        # Print number of samples per channel stored in file RIF.dat
        Ns = rif.nsamples("RIF.dat")

        print "Number of samples per channel =", Ns/2

    """
    return _rif.nsamples(filename)

def totalpower_single_channel(filename, dt):
    """Calculate total power.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a one dimensional array with shape=(Ns/(dt*2*20e6))
    where Ns is the total number of samples in the file.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate total power, averaged over 1s for data in file "RIF.dat"
        P = rif.totalpower_single channel("RIF.dat", 1.0)

        plt.plot(P) # Plot total power

        plt.show()
    
    """
    Ns = nsamples(filename)
    blocksize = int(dt * 20e6)
    nblocks = Ns / blocksize

    P = _rif.totalpower_single_channel(filename, blocksize, nblocks)

    return P

def averagespectrum_single_channel(filename, blocksize=4096):
    """Calculate dynamic spectrum.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *blocksize* number of samples to user for each FFT

    It returns a two dimensional array
    with shape=(blocksize) containing the average
    spectrum.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate average spectrum of data in "RIF.dat"
        # using a blocksize of 4096 samples for the FFT
        S = rif.averagespectrum("RIF.dat", 4096)

        # Plot the average spectrum
        plt.plot(S)
        
        plt.show()
        
    """
    Ns = nsamples(filename)
    navg = int(Ns / blocksize)

    S = _rif.averagespectrum_single_channel(filename, blocksize, navg)

    return S

def totalpower(filename, dt):
    """Calculate total power.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a two dimensional array with shape=(Ns/(dt*2*20e6), 2)
    where Ns is the total number of samples in the file.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate total power, averaged over 1s for data in file "RIF.dat"
        P = rif.totalpower("RIF.dat", 1.0)

        plt.plot(P[0]) # Plot total power of channel 0
        plt.plot(P[1]) # Plot total power of channel 1

        plt.show()
    
    """
    Ns = nsamples(filename)
    blocksize = int(dt * 20e6)
    nblocks = Ns / (2 * blocksize)

    P = _rif.totalpower(filename, blocksize, nblocks)

    return P.reshape((nblocks, 2)).transpose()

def averagespectrum(filename, blocksize=4096):
    """Calculate dynamic spectrum.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *blocksize* number of samples to user for each FFT

    It returns a two dimensional array
    with shape=(2, blocksize) containing the average
    spectrum for each channel.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate average spectrum of data in "RIF.dat"
        # using a blocksize of 4096 samples for the FFT
        S = rif.averagespectrum("RIF.dat", 4096)

        # Plot both spectra
        plt.subplot(2,1,1)
        plt.plot(S[0])
        plt.subplot(2,1,2)
        plt.plot(S{1])
        
        plt.show()
        
    """
    Ns = nsamples(filename)
    nf = blocksize / 2 + 1
    navg = int(Ns / (2*blocksize))

    S = _rif.averagespectrum(filename, blocksize, navg)

    return S.reshape((nf, 2)).transpose()

def dynamicspectrum(filename, dt, blocksize=4096):
    """Calculate dynamic spectrum.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *blocksize* number of samples to user for each FFT
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a three dimensional array
    with shape=(2, Ns/(dt*2*20e6), blocksize) containing the dynamic
    spectrum for each channel.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate dynamic spectrum of data in "RIF.dat" averaging over
        # 1s and using a blocksize of 4096 samples for the FFT
        S = rif.dynamicspectrum("RIF.dat", 1.0, 4096)

        # Plot all spectra in a single graph
        plt.subplot(2,1,1)
        for i in range(S.shape[0]):
            plt.plot(S[1][i])

        # Plot the dynamic spectrum
        plt.subplot(2,1,2)
        plt.imshow(S[1].transpose(), aspect='auto')
        
        plt.show()
        
    """
    Ns = nsamples(filename)
    nf = blocksize / 2 + 1
    navg = int(20e6 * dt / blocksize)
    nblocks = Ns / (2 * navg * blocksize)

    S = _rif.dynamicspectrum(filename, blocksize, navg, nblocks)

    return np.rollaxis(S.reshape((nblocks, nf, 2)), 2)

def crosscorrelation(filename, dt, blocksize=4096):
    """Calculate cross correlation.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *blocksize* number of samples to user for each FFT
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a two dimensional array with shape=(Ns/(dt*2*20e6))
    containing the real (e.g. cosine) and imaginary (e.g. sine) parts
    of the complex normalized cross correlation of the signal of the two
    channels.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate the cross correlation of data in "RIF.dat" averaging
        # over 1s and using a blocksize of 4096 samples for the FFT
        R = rif.crosscorrelation("RIF.dat", 1.0, 4096)

        plt.plot(R[0]) # Plot the real (e.g. cos) component
        plt.plot(R[1]) # Plot the imaginary (e.g. sin) component

        plt.show()
        
    """
    Ns = nsamples(filename)
    navg = int(20e6 * dt / blocksize)
    nblocks = Ns / (2 * navg * blocksize)

    R = _rif.crosscorrelation(filename, blocksize, navg, nblocks)

    return R.reshape((nblocks, 2)).transpose()

def crosscorrelation_parallel(filename, dt, blocksize=4096, nread=1):
    """Calculate cross correlation.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *blocksize* number of samples to user for each FFT
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a two dimensional array with shape=(Ns/(dt*2*20e6))
    containing the real (e.g. cosine) and imaginary (e.g. sine) parts
    of the complex normalized cross correlation of the signal of the two
    channels.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate the cross correlation of data in "RIF.dat" averaging
        # over 1s and using a blocksize of 4096 samples for the FFT
        R = rif.crosscorrelation("RIF.dat", 1.0, 4096)

        plt.plot(R[0]) # Plot the real (e.g. cos) component
        plt.plot(R[1]) # Plot the imaginary (e.g. sin) component

        plt.show()
        
    """
    Ns = nsamples(filename)
    navg = int(20e6 * dt / blocksize)
    nblocks = Ns / (2 * navg * blocksize)

    R = _rif.crosscorrelation_parallel(filename, blocksize, navg, nblocks, nread)

    return R.reshape((nblocks, 2)).transpose()

def frequencies(blocksize=4096):
    """Return frequencies corresponding to blocksize used to calculate FFT.
    """

    return np.linspace(1415.4, 1425.4, blocksize / 2 + 1)

def coordinates(filename, dt=1.0):
    """Returns coordinates RA and DEC interpolated according to the specified timestep.
    
    Returns a 4 dimensional array with the following rows: right ascension of telescope A,
    declination of telescope A, right ascension of telescope B and declination of telescope B.
    
    Example::

        # This example plots the results of a declination scan with telescope A
        
        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate total power, averaged over 0.5s for data in file "RIF.dat"
        P = rif.totalpower_single channel("RIF.dat", 0.5)

        # Calculate corresponding coordinates
        cds = rif.coordinates("RIF.cds", 0.5)
        
        # Declination coordinates for telescope A are stored in the second row
        DEC = cds[1]
        
        # Plot total power as a function of declination
        plt.plot(DEC, P)

        plt.show()
    
    """
    
    # Load datafile with coordinates
    cds = np.loadtxt(filename).transpose()
    
    # Get correct columns
    RA_A = cds[3]
    DEC_A = cds[4]
    RA_B = cds[5]
    DEC_B = cds[6]
    
    # Get number of elements in input and output
    N_in = cds.shape[1]
    N_out = int(np.floor(N_in / dt))
    
    # Get times for input points and centers of averaging steps
    t_in = np.arange(N_in)+0.5
    t_out = np.arange(N_out)*dt
    
    # Return interpolated coordinates
    return np.vstack(( \
      np.interp(t_out, t_in, RA_A),
      np.interp(t_out, t_in, DEC_A),
      np.interp(t_out, t_in, RA_B),
      np.interp(t_out, t_in, DEC_B)))

def sumsquared(filename, dt):
    """Calculate average of added channels squared.

    This function accepts the following parameters:

    *filename* name of the file to read data from
    *dt* time in seconds to average over (each second = 20e6 samples)

    It returns a one dimensional array with shape=(Ns/(dt*2*20e6), )
    where Ns is the total number of samples in the file.

    Example::

        # Import modules
        import rif
        import matplotlib.pyplot as plt

        # Calculate total power, averaged over 1s for data in file "RIF.dat"
        P = rif.sumsquared("RIF.dat", 1.0)

        plt.plot(P) # Plot

        plt.show()
    
    """
    Ns = nsamples(filename)
    blocksize = int(dt * 20e6)
    nblocks = Ns / (2 * blocksize)

    P = _rif.sumsquared(filename, blocksize, nblocks)

    return P
