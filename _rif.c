/**************************************************************************
 *  This file is part of the Radio Interferometer DAQ.                    *
 *  Copyright (C) 2010 Pim Schellart <P.Schellart@astro.ru.nl>            *
 *                                                                        *
 *  This library is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        * 
 *  This library is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.  *
 **************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>

/**
 * \brief Read data from RIF binary file.
 * 
 * \param filename name of file to read from
 * \param data array of size n used for output in the form of
 *        (channel 0 sample 0, channel 1 sample 1, channel 0 sample 1, ...)
 * \param start number of sample to start reading from
 * \param n number of samples to read
 */
void readdata(char* filename, char* data, long start, long n)
{
  FILE *fp;

  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  /* Skip to correct position */
  fseek(fp, start * sizeof(char), SEEK_SET);

  /* Read data */
  fread(data, sizeof(char), n, fp);

  /* Close file */
  fclose(fp);
}

/**
 * \brief Determines number of samples in the file
 *
 * \param filename name of file to read from
 *
 * \return number of samples in file (divide by two for number of samples
 *         per channel.
 */
long nsamples(char* filename)
{
  FILE *fp;
  long n;

  /* Open file */
  fp = fopen(filename, "rb");

  /* Find end of file */
  fseek(fp, 0, SEEK_END);
  
  /* Get size */
  n = ftell(fp) / sizeof(char);
  
  /* Close file */
  close(fp);
  
  return n;
}

/**
 * \brief Calculate total power as a function of time (averaged over
 *        specified number of samples.
 *
 * \param filename name of file to read data from
 * \param P array of length nblocks used for output of total power
 * \param blocksize number of samples to average for 1 datapoint in output
 * \param nblocks number of points in output
 */
void totalpower_single_channel(char* filename, double* P, long blocksize, long nblocks)
{
  long i, j;
  long n;
  char *buffer, *b;
  
  FILE *fp;

  /* Allocate memory for buffer */
  n = blocksize;
  buffer = malloc(n * sizeof(char));
  
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the number of samples used for averaging.\n");
    return;
  }

  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  for (i=nblocks; i!=0; --i)
  {
    /* Read data into buffer */
    fread(buffer, sizeof(char), n, fp);
    
    /* Square and sum */
    *P = 0;
    
    b = buffer;
    for (j=blocksize; j!=0; --j)
    {
      *P += ((double) (*b)) *  ((double) (*b));
      ++b;
    }

    /* Mean */
    *P /= blocksize;
    ++P;
  }

  /* Close file */
  fclose(fp);
  
  /* Free memory */
  free(buffer);
}

/**
 * \brief Calculates average spectrum (average FFT for a single channel).
 *
 * \param filename name of file to read data from
 * \param S array of length blocksize for output of dynamic
 *        spectrum.
 * \param blocksize number of samples used for FFT (should be power of 2)
 * \param navg number of blocks to average (after taking absolute value
 *        of transform) for a single spectrum in the output
 * \param nblocks number of points in output
 */
void averagespectrum_single_channel(char* filename, double* S, long blocksize, long navg)
{
  long i, j, k, nf;
  char* buffer;
  double *in;
  fftw_complex *out;
  fftw_plan p;
  FILE *fp;

  /* Number of frequencies in output */
  nf = blocksize / 2 + 1;

  /* Allocate memory for FFT */
  in = (double*) fftw_malloc(sizeof(double) * blocksize);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);

  /* Create plan for FFT */
  p = fftw_plan_dft_r2c_1d(blocksize, in, out, FFTW_ESTIMATE);

  /* Allocate memory for buffer */
  buffer = malloc(blocksize * sizeof(char));
  
  printf("size: %ld\n", sizeof(char));
  printf("size: %ld\n", blocksize);
  
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the blocksize.\n");
    return;
  }
  
  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  /* Clear output array */
  for (i=0; i<nf; i++)
  {
    S[i] = 0;
  }
  
  /* Loop over blocks to average */
  for (j=0; j<navg; j++)
  {
    /* Read data into buffer */
    fread(buffer, sizeof(char), blocksize, fp);

    /* Put data into FFT input arrays */
    for (k=0; k<blocksize; k++)
    {
      in[k] = buffer[k];
    }

    /* Perform FFT */
    fftw_execute(p);

    /* Get data from FFT output arrays */
    for (k=0; k<nf; k++)
    {
      S[k] += cabs(out[k]);
    }
  }

  /* Average */
  for (i=0; i<nf; i++)
  {
    S[i] /= navg;
  }
  
  /* Close file */
  fclose(fp);

  /* Free memory */
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  free(buffer);
}

/**
 * \brief Calculate total power as a function of time (averaged over
 *        specified number of samples.
 *
 * \param filename name of file to read data from
 * \param P array of length nblocks*2 used for output of total power in the
 *        form of (average (block 0 channel 0), average (block 0 channel 1), *        average (block 1 channel 0), ...
 * \param blocksize number of samples to average for 1 datapoint in output
 * \param nblocks number of points in output
 */
void totalpower(char* filename, double* P, long blocksize, long nblocks)
{
  long i, j;
  long n;
  double *P0, *P1;
  char *buffer, *b0, *b1;
  FILE *fp;

  /* Allocate memory for buffer */
  n = blocksize * 2;
  buffer = malloc(n * sizeof(char));
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the number of samples used for averaging.\n");
    return;
  }

  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  P0 = P1 = P; ++P1;
  for (i=nblocks; i!=0; --i)
  {
    /* Feedback */
    printf("processing block %ld of %ld\n", i, nblocks);
    
    /* Read data into buffer */
    fread(buffer, sizeof(char), n, fp);
    
    /* Square and sum */
    *P0 = 0;
    *P1 = 0;
    
    b0 = b1 = buffer; ++b1;
    for (j=blocksize; j!=0; --j)
    {
      *P0 += (double) *b0 * *b0;
      *P1 += (double) *b1 * *b1;
      
      b0 += 2;
      b1 += 2;
    }

    /* Mean */
    *P0 /= blocksize;
    *P1 /= blocksize;
    
    P0 += 2;
    P1 += 2;
  }

  /* Close file */
  fclose(fp);
  
  /* Free memory */
  free(buffer);
}

/**
 * \brief Calculates average spectrum (average FFT for both channels).
 *
 * \param filename name of file to read data from
 * \param S array of length blocksize*2 for output of dynamic
 *        spectrum.
 * \param blocksize number of samples used for FFT (should be power of 2)
 * \param navg number of blocks to average (after taking absolute value
 *        of transform) for a single spectrum in the output
 * \param nblocks number of points in output
 */
void averagespectrum(char* filename, double* S, long blocksize, long navg)
{
  long i, j, k, k0, k1, nf, idx0, idx1, n;
  char* buffer;
  double *in0, *in1;
  fftw_complex *out0, *out1;
  fftw_plan p0, p1;
  FILE *fp;

  /* Number of frequencies in output */
  nf = blocksize / 2 + 1;

  /* Allocate memory for FFT */
  in0 = (double*) fftw_malloc(sizeof(double) * blocksize);
  in1 = (double*) fftw_malloc(sizeof(double) * blocksize);
  out0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);
  out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);

  /* Create plan for FFT */
  p0 = fftw_plan_dft_r2c_1d(blocksize, in0, out0, FFTW_ESTIMATE);
  p1 = fftw_plan_dft_r2c_1d(blocksize, in1, out1, FFTW_ESTIMATE);

  /* Allocate memory for buffer */
  n = blocksize * 2;
  buffer = malloc(n * sizeof(char));
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the blocksize.\n");
    return;
  }
  
  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  /* Clear output array */
  for (i=0; i<nf*2; i++)
  {
    S[i] = 0;
  }
  
  /* Loop over blocks to average */
  for (j=0; j<navg; j++)
  {
    /* Read data into buffer */
    fread(buffer, sizeof(char), n, fp);

    /* Put data into FFT input arrays */
    k0 = 0;
    k1 = 1;
    for (k=0; k<blocksize; k++)
    {
      in0[k] = buffer[k0];
      in1[k] = buffer[k1];
      
      k0 += 2;
      k1 += 2;
    }

    /* Perform FFT */
    fftw_execute(p0);
    fftw_execute(p1);

    /* Get data from FFT output arrays */
    idx0 = 0;
    idx1 = 1;
    for (k=0; k<nf; k++)
    {
      S[idx0] += cabs(out0[k]);
      S[idx1] += cabs(out1[k]);
      
      idx0 += 2;
      idx1 += 2;
    }
  }

  /* Average */
  for (i=0; i<nf*2; i++)
  {
    S[i] = S[i] / navg;
  }
  
  /* Close file */
  fclose(fp);

  /* Free memory */
  fftw_destroy_plan(p0);
  fftw_destroy_plan(p1);
  fftw_free(in0);
  fftw_free(out0);
  fftw_free(in1);
  fftw_free(out1);
  free(buffer);
}

/**
 * \brief Calculates dynamic spectrum (e.g. average FFT for both channels as a function of time).
 *
 * \param filename name of file to read data from
 * \param S array of length nblocks*blocksize*2 for output of dynamic
 *        spectrum.
 * \param blocksize number of samples used for FFT (should be power of 2)
 * \param navg number of blocks to average (after taking absolute value
 *        of transform) for a single spectrum in the output
 * \param nblocks number of points in output
 */
void dynamicspectrum(char* filename, double* S, long blocksize, long navg, long nblocks)
{
  long i, j, k, k0, k1, nf, idx, idx0, idx1, n;
  char* buffer;
  double *in0, *in1;
  fftw_complex *out0, *out1;
  fftw_plan p0, p1;
  FILE *fp;

  /* Number of frequencies in output */
  nf = blocksize / 2 + 1;

  /* Allocate memory for FFT */
  in0 = (double*) fftw_malloc(sizeof(double) * blocksize);
  in1 = (double*) fftw_malloc(sizeof(double) * blocksize);
  out0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);
  out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);

  /* Create plan for FFT */
  p0 = fftw_plan_dft_r2c_1d(blocksize, in0, out0, FFTW_ESTIMATE);
  p1 = fftw_plan_dft_r2c_1d(blocksize, in1, out1, FFTW_ESTIMATE);

  /* Allocate memory for buffer */
  n = blocksize * navg * 2;
  buffer = malloc(n * sizeof(char));
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the number of samples used for averaging.\n");
    return;
  }
  
  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  /* Clear output array */
  for (i=0; i<nblocks*nf*2; i++)
  {
    S[i] = 0;
  }

  /* Loop over output datapoints */
  idx = 0;
  for (i=0; i<nblocks; i++)
  {
    /* Read data into buffer */
    fread(buffer, sizeof(char), n, fp);
    
    /* Loop over blocks to average */
    k0 = 0;
    k1 = 1;
    for (j=0; j<navg; j++)
    {
      /* Put data into FFT input arrays */
      for (k=0; k<blocksize; k++)
      {
        in0[k] = buffer[k0];
        in1[k] = buffer[k1];
        
        k0 += 2;
        k1 += 2;
      }

      /* Perform FFT */
      fftw_execute(p0);
      fftw_execute(p1);

      /* Get data from FFT output arrays */
      idx0 = idx;
      idx1 = idx+1;
      for (k=0; k<nf; k++)
      {
        S[idx0] += cabs(out0[k]);
        S[idx1] += cabs(out1[k]);
        
        idx0 += 2;
        idx1 += 2;
      }
    }
    
    idx += 2*nf;
  }

  /* Average */
  for (i=0; i<nblocks*nf*2; i++)
  {
    S[i] = S[i] / navg;
  }
  
  /* Close file */
  fclose(fp);

  /* Free memory */
  fftw_destroy_plan(p0);
  fftw_destroy_plan(p1);
  fftw_free(in0);
  fftw_free(out0);
  fftw_free(in1);
  fftw_free(out1);
  free(buffer);
}

/**
 * \brief Calculates real (e.g. cos part) and imaginary (e.g. sin part) of
 *        the cross correlation of the signals of the two channels.
 *
 * \param filename name of file to read data from
 * \param R cross correlation of signals (cos, sin, cos, sin, ...)
 * \param blocksize number of samples to use for FFT
 * \param navg number of blocks to average
 * \param nblocks number of points in output
 */
void crosscorrelation(char* filename, double* R, long blocksize, long navg, long nblocks)
{
  long i, j, k, nf, iidx, oidx;
  double c0, c1, re, im, norm, R0, R1;
  complex corr;
  char* buffer;

  double *in0, *in1;
  fftw_complex *out0, *out1;
  fftw_plan p0, p1;

  /* Number of frequencies in output */
  nf = blocksize / 2 + 1;

  /* Allocate memory for FFT */
  in0 = (double*) fftw_malloc(sizeof(double) * blocksize);
  in1 = (double*) fftw_malloc(sizeof(double) * blocksize);
  out0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);
  out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nf);

  /* Create plan for FFT */
  p0 = fftw_plan_dft_r2c_1d(blocksize, in0, out0, FFTW_ESTIMATE);
  p1 = fftw_plan_dft_r2c_1d(blocksize, in1, out1, FFTW_ESTIMATE);

  /* Allocate memory for buffer */
  buffer = malloc(blocksize*navg*2*sizeof(char));
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the number of samples used for averaging.\n");
    return;
  }

  /* Loop over output datapoints */
  for (i=0; i<nblocks; i++)
  {
    /* Data stored cos, sin, cos, sin ... so we need index to pair */
    oidx = 2*i;

    /* Clear output array */
    R[oidx] = 0;
    R[oidx+1] = 0;

    /* Read data into buffer */
    readdata(filename, buffer, i*blocksize*navg*2, blocksize*navg*2);

    /* Loop over blocks to average */
    for (j=0; j<navg; j++)
    {
      /* Put data into FFT input arrays */
      for (k=0; k<blocksize; k++)
      {
        iidx = j*blocksize+2*k;
        in0[k] = buffer[iidx];
        in1[k] = buffer[iidx+1];
      }

      /* Perform FFT */
      fftw_execute(p0);
      fftw_execute(p1);

      /* Loop over frequencies */
      c0 = 0; c1 = 0, R0 = 0; R1 = 0;
      for (k=0; k<nf; k++)
      {
        /* First normalize FFT */
        out0[k] = out0[k] / nf;
        out1[k] = out1[k] / nf;

        /* Then calculate normalization factors */
        re = creal(out0[k]); im = cimag(out0[k]);
        c0 += re * re + im * im;

        re = creal(out1[k]); im = cimag(out1[k]);
        c1 += re * re + im * im;

        /* Calculate cross correlation */
        corr = out0[k] * conj(out1[k]);

        R0 += creal(corr);
        R1 += cimag(corr);
      }

      /* Normalize cross correlation */
      norm = sqrt(c0) * sqrt(c1);

      R[oidx] += R0 / norm;
      R[oidx+1] += R1 / norm;
    }

    /* Average over time */
    R[oidx] /= navg;
    R[oidx+1] /= navg;
  }

  /* Free memory */
  fftw_destroy_plan(p0);
  fftw_destroy_plan(p1);
  fftw_free(in0);
  fftw_free(out0);
  fftw_free(in1);
  fftw_free(out1);
  free(buffer);
}

/**
 * \brief Calculate total power as a function of time (averaged over
 *        specified number of samples.
 *
 * \param filename name of file to read data from
 * \param P array of length nblocks used for output of total power in the
 *        form of (average, average, ...)
 * \param blocksize number of samples to average for 1 datapoint in output
 * \param nblocks number of points in output
 */
void sumsquared(char* filename, double* P, long blocksize, long nblocks)
{
  long i, j;
  long n;
  double tmp;
  char *buffer, *b0, *b1;
  FILE *fp;

  /* Allocate memory for buffer */
  n = blocksize * 2;
  buffer = malloc(n * sizeof(char));
  if (buffer == NULL)
  {
    printf("Cannot allocate memory for read-in buffer, try reducing the number of samples used for averaging.\n");
    return;
  }

  /* Open file */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    printf("Error: could not open file!\n");
    return;
  }

  for (i=nblocks; i!=0; --i)
  {
    /* Feedback */
    printf("processing block %ld of %ld\n", i, nblocks);
    
    /* Read data into buffer */
    fread(buffer, sizeof(char), n, fp);
    
    /* Add, square and sum */
    *P = 0;
    
    b0 = b1 = buffer; ++b1;
    for (j=blocksize; j!=0; --j)
    {
      tmp = ((double) *b0) + ((double) *b1);
      *P += tmp * tmp;
      
      b0 += 2;
      b1 += 2;
    }

    /* Mean */
    *P /= blocksize;
    
    ++P;
  }

  /* Close file */
  fclose(fp);
  
  /* Free memory */
  free(buffer);
}