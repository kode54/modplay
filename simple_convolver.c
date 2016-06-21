/* vim: set et ts=4 sw=4 sts=4 ft=c:
 *
 * Copyright (C) 2016 Christopher Snowhill.  All rights reserved.
 * https://github.com/kode54/fft-resampler
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* Welcome to FFT convolution, with your host, kode54! Today, we'll be following
 * along with some code I mostly picked up from examples I discovered scattered
 * around various parts of the Internet. I chose to use kissfft because it's
 * fairly fast, and also under a less restrictive license than FFTW. */

#include "simple_convolver.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Here comes the magic import header! */

#include "kissfft/kiss_fftr.h"

/* Only the simplest of state information is necessary for this, and it's
 * designed around a simple usage case. As many samples as you can generate,
 * input one at a time, and output pulled one at a time as well. */

typedef struct convolver_state
{
	int fftlen;                             /* size of FFT */
	int stepsize;                           /* size of overlapping steps */
	int buffered_in;                        /* how many input samples buffered */
	int buffered_out;                       /* output samples buffered */
	int inputs;                             /* Input channels */
	int outputs;                            /* Output channels */
	int impulses;				/* Impulse count */
	int channels_per_impulse;               /* Channels per impulse */
	kiss_fftr_cfg cfg_fw, cfg_bw;           /* forward and backwards instances */
	kiss_fft_cpx *f_in, *f_out, **f_ir;     /* output and impulse in frequency domain */
	float *revspace, **outspace, **inspace; /* reverse, output, and input work space */
} convolver_state;

/* Fully opaque convolver state created and returned here, otherwise NULL on
 * failure. Users are welcome to change this to pass in a const pointer to an
 * impulse and its size, which will be copied and no longer needed upon return. */

void * convolver_create(const float * const* impulse, const int impulse_size, int input_channels, int output_channels, int impulse_count, int channels_per_impulse)
{
	convolver_state * state;
	int fftlen, total_channels, i, j, k;

	if ( input_channels > impulse_count )
		return NULL;

	/* Output channels must be an even multiple or factor of input channels */

	if ( (output_channels % input_channels) != 0 &&
	     (input_channels % output_channels) != 0 )
		return NULL;

	/* Impulse count must be an even multiple or factor of input channels */

	if ( (impulse_count % input_channels) != 0 &&
	     (input_channels % impulse_count) != 0 )
		return NULL;

	/* If mixing impulses into an output, output count must be an even
	 * factor of impulses. */

	if ( (output_channels < (impulse_count * channels_per_impulse)) &&
	     ((impulse_count * channels_per_impulse) % output_channels) != 0 )
		return NULL;

	/* Otherwise, if using less impulses than output channels, there must
	 * be an even multiple of output channels per impulse. */

	if ( (output_channels > (impulse_count * channels_per_impulse)) &&
	     (output_channels % (impulse_count * channels_per_impulse)) )
		return NULL;

	state = (convolver_state *) calloc( 1, sizeof( convolver_state ) );

	if ( !state )
		return NULL;

	state->inputs = input_channels;
	state->outputs = output_channels;
	state->impulses = impulse_count;
	state->channels_per_impulse = channels_per_impulse;
	total_channels = impulse_count * channels_per_impulse;

    /* This is bog standard, from the example code I lifted. Mainly needs this
     * calculation for varying impulse sizes */

	fftlen = ( impulse_size * 3 ) / 2 + 10000;
	{
		// round up to a power of two
		int pow = 1;
		while ( fftlen > 2 ) { pow++; fftlen /= 2; }
		fftlen = 2 << pow;
	}

    /* And a bog standard minimum size to work with, I guess. */
    
	if ( fftlen < 32768 )
		fftlen = 32768;

	if ( fftlen > 1024 + impulse_size + 10 )
	{
		int current;

		fftlen = 1024 + impulse_size + 10;

#define TRYFACTOR(n) if ( failed && current*n > fftlen ) { current *= n; failed = 0; }
		
		current = 1;
		while ( current < fftlen )
		{
		    /* These are optimal FFT scale factors that kissfft can handle. */
		    
			int failed = 1;
			TRYFACTOR(2)
			TRYFACTOR(3)
			TRYFACTOR(5)
			TRYFACTOR(6)
			TRYFACTOR(9)
			TRYFACTOR(10)
			TRYFACTOR(12)
			TRYFACTOR(15)
			TRYFACTOR(18)
			TRYFACTOR(20)
			TRYFACTOR(25)
			TRYFACTOR(30)
			TRYFACTOR(36)
			TRYFACTOR(45)
			TRYFACTOR(75)
			if ( failed )
				current *= 4;
		}
#undef TRYFACTOR

		fftlen = current;
	}

	state->fftlen = fftlen;
	state->stepsize = fftlen - impulse_size - 10;
	state->buffered_in = 0;
	state->buffered_out = 0;

    /* Prepare arrays for multiple inputs */
    /* And we use kissfft's aligned malloc functions/macros to allocate these things. */
 
	if ( (state->f_in = (kiss_fft_cpx*) KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen/2+1))) == NULL )
		goto error;

	if ( (state->f_out = (kiss_fft_cpx*) KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen/2+1))) == NULL )
		goto error;

	if ( (state->f_ir = (kiss_fft_cpx**) calloc(sizeof(kiss_fft_cpx*), total_channels)) == NULL )
		goto error;
	for (i = 0; i < total_channels; ++i)
	{
		if ( (state->f_ir[i] = (kiss_fft_cpx*) KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen/2+1))) == NULL )
			goto error;
	}

	if ( (state->revspace = (float *) KISS_FFT_MALLOC(sizeof(float) * fftlen)) == NULL )
		goto error;

	if ( (state->outspace = (float **) calloc(sizeof(float *), output_channels)) == NULL )
		goto error;
	for (i = 0; i < output_channels; ++i)
	{
		if ( (state->outspace[i] = (float *) calloc(sizeof(float), fftlen)) == NULL )
			goto error;
	}

	if ( (state->inspace = (float **) calloc(sizeof(float *), input_channels)) == NULL )
		goto error;
	for (i = 0; i < input_channels; ++i)
	{
		if ( (state->inspace[i] = (float *) calloc(sizeof(float), fftlen)) == NULL )
			goto error;
	}

	if ( (state->cfg_fw = kiss_fftr_alloc( fftlen, 0, NULL, NULL )) == NULL )
		goto error;
	if ( (state->cfg_bw = kiss_fftr_alloc( fftlen, 1, NULL, NULL )) == NULL )
		goto error;

	convolver_restage( state, impulse );

	return state;

error:
	convolver_delete( state );
	return NULL;
}

/* Restage the convolver with a new impulse set, same size/parameters */
void convolver_restage( void * state_, const float * const * impulse )
{
	convolver_state * state = (convolver_state *) state_;

	float * impulse_temp;

	int impulse_count = state->impulses;
	int channels_per_impulse = state->channels_per_impulse;
	int fftlen = state->fftlen;
	int impulse_size = fftlen - state->stepsize - 10;
	int i, j, k;

	kiss_fftr_cfg cfg_fw = state->cfg_fw;
	kiss_fft_cpx ** f_ir = state->f_ir;

    /* Since the FFT requires a full input for every transformaton, we allocate
     * a temporary buffer, which we fill with the impulse, then pad with silence. */
     
	if ( (impulse_temp = (float*) malloc( sizeof(float) * fftlen )) == NULL )
		return;

	memset( impulse_temp + impulse_size, 0, sizeof(float) * (fftlen - impulse_size) );

	for (i = 0; i < impulse_count; ++i)
	{
		for (j = 0; j < channels_per_impulse; ++j)
		{
			for (k = 0; k < impulse_size; ++k)
			{
				impulse_temp[k] = impulse[i][j + k * channels_per_impulse];
			}

    /* Our first actual transformation, which is cached for the life of this convolver. */
    
			kiss_fftr( cfg_fw, impulse_temp, f_ir[i * channels_per_impulse + j] );
		}
	}

	free( impulse_temp );
}

/* Delete our opaque state, by freeing all of its member structures, then the
 * top level structure itself. */

void convolver_delete( void * state_ )
{
	if ( state_ )
	{
		int i, input_channels, output_channels, total_channels;
		convolver_state * state = (convolver_state *) state_;
		input_channels = state->inputs;
		output_channels = state->outputs;
		total_channels = state->impulses * state->channels_per_impulse;
		if ( state->cfg_fw )
			kiss_fftr_free( state->cfg_fw );
		if ( state->cfg_bw )
			kiss_fftr_free( state->cfg_bw );
		if ( state->f_ir )
		{
			for (i = 0; i < total_channels; ++i)
			{
				if ( state->f_ir[i] )
					KISS_FFT_FREE( state->f_ir[i] );
			}
			free( state->f_ir );
		}
		if ( state->f_out )
			KISS_FFT_FREE( state->f_out );
		if ( state->f_in )
			KISS_FFT_FREE( state->f_in );
		if ( state->revspace )
			KISS_FFT_FREE( state->revspace );
		if ( state->outspace )
		{
			for (i = 0; i < output_channels; ++i)
			{
				if ( state->outspace[i] )
					free( state->outspace[i] );
			}
			free( state->outspace );
		}
		if ( state->inspace )
		{
			for (i = 0; i < input_channels; ++i)
			{
				if ( state->inspace[i] )
					free( state->inspace[i] );
			}
			free( state->inspace );
		}
		free( state );
	}
}

/* This resets the state between uses, if you need to restart output on startup. */

void convolver_clear( void * state_ )
{
	if ( state_ )
	{
	    /* Clearing for a new use setup only requires resetting the input and
	     * output buffers, not actually changing any of the FFT state. */
	     
		int i, input_channels, total_channels, fftlen;
		convolver_state * state = (convolver_state *) state_;
		input_channels = state->inputs;
		total_channels = state->impulses * state->channels_per_impulse;
		fftlen = state->fftlen;
		state->buffered_in = 0;
		state->buffered_out = 0;
		for (i = 0; i < input_channels; ++i)
			memset( state->inspace[i], 0, sizeof(float) * fftlen );
		for (i = 0; i < total_channels; ++i)
			memset( state->outspace[i], 0, sizeof(float) * fftlen );
	}
}

/* This returns how many slots are available in the input buffer, before the
 * user may pull output. */

int convolver_get_free_count(void * state_)
{
	if ( state_ )
	{
		convolver_state * state = (convolver_state *) state_;
		return state->stepsize - state->buffered_in;
	}
	return 0;
}

/* This returns how many output samples are currently buffered, or zero if not
 * enough sample data has been buffered. */

int convolver_ready(void * state_)
{
	if ( state_ )
	{
		convolver_state * state = (convolver_state *) state_;
		return state->buffered_out;
	}
	return 0;
}

/* Input sample data is fed in here, one sample at a time. */

void convolver_write(void * state_, const float * input_samples)
{
	if ( state_ )
	{
		int i, j, k, input_channels;
		convolver_state * state = (convolver_state *) state_;
		input_channels = state->inputs;

		for (i = 0; i < input_channels; ++i)
			state->inspace[i][ state->buffered_in ] = input_samples[i];
		++state->buffered_in;
		
		/* And every stepsize samples buffered, it convolves a new block of samples. */
		
		if ( state->buffered_in == state->stepsize )
		{
			int output_channels = state->outputs;
			int impulse_count = state->impulses;
			int channels_per_impulse = state->channels_per_impulse;
			int total_channels = impulse_count * channels_per_impulse;
			int fftlen;
			int index;
			float fftlen_if;
			kiss_fft_cpx *f_in = state->f_in;
			kiss_fft_cpx *f_out = state->f_out;
			float *revspace = state->revspace;
			int inputs_per_impulse = (input_channels + impulse_count - 1) / impulse_count;
			int input_start, input_end;
			int outputs_per_input = (output_channels + input_channels - 1) / input_channels;
			int inputs_per_output = (output_channels + input_channels - 1) / output_channels;

			for ( i = 0; i < impulse_count; ++i )
			{
				int input_index, output_index;
				int output_start, output_end;
				input_start = i * input_channels / impulse_count;

				for (input_index = input_start, input_end = input_start + inputs_per_impulse; input_index < input_end; ++input_index)
				{
            /* First the input samples are transformed to frequency domain, like
             * the cached impulse was in the setup function. */
             
					kiss_fftr( state->cfg_fw, state->inspace[input_index], f_in );

					for ( j = 0; j < channels_per_impulse; ++j )
					{
            /* Then we cross multiply the products of the frequency domain, the
             * real and imaginary values, into output real and imaginary pairs. */
 
						int index = i * channels_per_impulse + j; 
						kiss_fft_cpx *f_ir = state->f_ir[index];
						float *outspace;

						for ( k = 0, fftlen = state->fftlen / 2 + 1; k < fftlen; ++k )
						{
							float re = f_ir[k].r * f_in[k].r - f_ir[k].i * f_in[k].i;
							float im = f_ir[k].i * f_in[k].r + f_ir[k].r * f_in[k].i;
							f_out[k].r = re;
							f_out[k].i = im;
						}

            /* Then we transform back from frequency to time domain. */
            
						kiss_fftri( state->cfg_bw, f_out, revspace );

            /* Then we add the entire revspace block onto our output, dividing
             * each value by the total number of samples in the buffer. Remember,
             * since there is some overlap, this addition step is important. */

						output_start = input_index * output_channels / inputs_per_output;
						output_end = output_start + outputs_per_input;
						for (output_index = output_start; output_index < output_end; ++output_index)
						{
							outspace = state->outspace[output_index];
							for (k = 0, fftlen = state->fftlen, fftlen_if = 1.0f / (float)fftlen; k < fftlen; ++k)
								outspace[k] += revspace[k] * fftlen_if;
						}
					}
				}
			}

            /* Output samples are now buffered and ready for retrieval. */
            
			state->buffered_out = state->stepsize;
			state->buffered_in = 0;
		}
	}
}

/* Call this to retrieve available output samples. */

void convolver_read(void *state_, float * output_samples)
{
	if ( state_ )
	{
		convolver_state * state = (convolver_state *) state_;

        /* We only want to pass in here if we have buffered samples. */

		if ( state->buffered_out )
		{
			int i, output_channels = state->outputs;

			for (i = 0; i < output_channels; ++i)
			{
				float sample = state->outspace[i][state->stepsize - state->buffered_out];

				output_samples[i] = sample;

            /* And if the buffer has run empty again, we want to slide back the buffer
             * contents, eliminating the space used by the block of output samples, and
             * filling in the space at the end with silence. */
			}

			if ( --state->buffered_out == 0 )
			{
				for (i = 0; i < output_channels; ++i)
				{
					float *outspace = state->outspace[i];
					memmove( outspace, outspace + state->stepsize, (state->fftlen - state->stepsize) * sizeof(float) );
					memset( outspace + state->fftlen - state->stepsize, 0, state->stepsize * sizeof(float) );
				}
			}
		}
		else
		{
			memset( output_samples, 0, state->outputs * sizeof(float) );
		}
	}
}
