#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package,top_path)

    config.add_extension('_rif',
                         sources=['_rif.pyf','_rif.c'],
                         libraries=['m', 'fftw3', 'gomp'],
                         extra_compile_args=['-ffast-math', '-fopenmp'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)

