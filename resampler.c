#include "resampler.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "simple_convolver.h"

enum { FILTER_WIDTH = 8192 };
enum { IIR_ORDER = 52 };

#if 0
static int fEqual(const float b, const float a)
{
    return fabs(a - b) < 1.0e-6;
}

static float sinc(float x)
{
    return fEqual(x, 0.0) ? 1.0 : sin(x * M_PI) / (x * M_PI);
}

static void gen_sinc( int width, int count, float* out )
{
	int half_count = count >> 1;
	int half_width = width >> 1;
	double fw = (double)(half_width);
	double dx = fw / (double)(half_count);
	double dy = 1.0 / (double)(half_count);
	double x = 0.0;
	double y = 0.0;
	int index_out = half_count;
	int index_in = half_count;
	double total = 0.0;
	out[0] = 0.0;
	while ( half_count-- )
	{
#if 0
        	// Blackman
        	float window = 0.42659 - 0.49656 * cos(M_PI + M_PI * y) + 0.076849 * cos(2.0 * M_PI * y);
#elif 1
	        // Nuttal 3 term
	        float window = 0.40897 + 0.5 * cos(M_PI * y) + 0.09103 * cos(2.0 * M_PI * y);
#elif 0
        	// C.R.Helmrich's 2 term window
	        float window = 0.79445 * cos(0.5 * M_PI * y) + 0.20555 * cos(1.5 * M_PI * y);
#elif 0
	        // Lanczos
        	float window = sinc(y);
#endif
		float sinc_ = fabs(x) < fw ? sinc(x) : 0.0;
		float windowed_sinc = sinc_ * window;
		total += windowed_sinc * ((index_in != index_out) ? 2.0 : 1.0);
		out[index_in] = out[index_out] = windowed_sinc;
		++index_out;
		--index_in;
		x += dx;
		y += dy;
	}
	total = 1.0 / total;
	while ( count-- )
	{
		*out++ *= total;
	}
}

static void lowpass_set_rate( float * out, unsigned int width, double ratio )
{
	double const filter = (ratio < 1.0) ? 1.0 : 1.0 / ratio;

	gen_sinc( (int) (width * filter + 1) & ~1, (int) width, out );
}
#elif 0
static void gen_sinc( double rolloff, int width, double offset, double spacing, double scale,
		int count, float* out )
{
	double const maxh = 2048;
	double const step = M_PI / maxh * spacing;
	double const to_w = maxh * 2 / width;
	double const pow_a_n = pow( rolloff, maxh );
	scale /= maxh * 2;
	double angle = (count / 2 - 1 + offset) * -step;

	while ( count-- )
	{
		*out++ = 0;
		double w = angle * to_w;
		if ( fabs( w ) < M_PI )
		{
			double rolloff_cos_a = rolloff * cos( angle );
			double num = 1 - rolloff_cos_a -
					pow_a_n * cos( maxh * angle ) +
					pow_a_n * rolloff * cos( (maxh - 1) * angle );
			double den = 1 - rolloff_cos_a - rolloff_cos_a + rolloff * rolloff;
			double sinc = scale * num / den - scale;

			out [-1] = (float) (cos( w ) * sinc + sinc);
		}
		angle += step;
	}
}

static void lowpass_set_rate( float * out, unsigned int width, double ratio )
{
	double const rolloff = 0.999;
	double const gain = 1.0;

	double const filter = (ratio < 1.0) ? 1.0 : 1.0 / ratio;
	double pos = 0.0;

	if ( filter < 1.0 )
		gen_sinc( rolloff, (int) (width * filter + 1) & ~1, pos, filter,
			(gain * filter), (int) width, out );
	else
	{
		memset(out, 0, sizeof(*out) * width);
		out[width / 2] = 1.0f;
	}
}
#else
typedef struct iir
{
    double cutoff;              //frequency cutoff
    double quality;             //frequency response quality
    double gain;                //peak gain
    double a0, a1, a2, b1, b2;  //coefficients
    double z1, z2;              //second-order IIR
} iir;

static void iir_reset(iir * i, double cutoff, double quality, double gain)
{
    double v, k, q, n;
    
    i->cutoff = cutoff;
    i->quality = quality;
    i->gain = gain;
    
    v = pow(10, fabs(gain) / 20.0);
    k = tan(M_PI * cutoff);
    q = quality;
    
    n = 1 / (1 + k / q + k * k);
    i->a0 = k * k * n;
    i->a1 = 2 * i->a0;
    i->a2 = i->a0;
    i->b1 = 2 * (k * k - 1) * n;
    i->b2 = (1 - k / q + k * k) * n;
}

static void iir_clear(iir * i)
{
    i->z1 = 0.0;
    i->z2 = 0.0;
}

static double iir_process(iir * i, double in)
{
    double out = in * i->a0 + i->z1;
    i->z1 = in * i->a1 + i->z2 - i->b1 * out;
    i->z2 = in * i->a2 - i->b2 * out;
    return out;
}

static double butterworth(unsigned int order, unsigned int phase)
{
    return -0.5 / cos(M_PI / 2.0 * (1.0 + (1.0 + (2.0 * phase + 1.0) / order)));
}

static void lowpass_set_rate( float * out, unsigned int width, double ratio )
{
    double const filter = (ratio < 1.0) ? 0.49 : (1.0 / ratio) * 0.49;
    
    iir iirf[IIR_ORDER / 2];

    int i, j;
    
    double sample;
    
    memset(iirf, 0, sizeof(iirf));
    
    for (i = 0, j = IIR_ORDER / 2; i < j; ++i)
        iir_reset(iirf + i, filter, butterworth(IIR_ORDER, i), 0.0);
    
    sample = 1.0;
    
    while ( width-- )
    {
        for (i = 0, j = IIR_ORDER / 2; i < j; ++i)
            sample = iir_process(iirf + i, sample);
        
        *out++ = sample;
        
        sample = 1e-25;
    }
}
#endif

enum { history_width = FILTER_WIDTH };

typedef struct resampler
{
	float history[history_width];
	double ratio;
	double position;
	int history_write;
	int history_read;
	void * convolver_in;
	void * convolver_out;
} resampler;

void * resampler_create()
{
	resampler * r = (resampler *) calloc(1, sizeof(resampler));

	return (void *) r;
}

void resampler_delete(void * r_)
{
	resampler * r = (resampler *) r_;
	if (r)
	{
		if (r->convolver_in)
			convolver_delete(r->convolver_in);
		if (r->convolver_out)
			convolver_delete(r->convolver_out);
		free(r);
	}
}

void resampler_set_ratio(void * r_, double ratio)
{
	resampler * r = (resampler *) r_;

	if (r->ratio != ratio)
	{
		float impulse[FILTER_WIDTH];
		float *impulse_ptr = impulse;
		r->ratio = ratio;
		lowpass_set_rate(impulse, FILTER_WIDTH, ratio);
		if (!r->convolver_in)
			r->convolver_in = convolver_create( (const float * const *) &impulse_ptr, FILTER_WIDTH, 1, 1, 1, 1 );
		else
			convolver_restage( r->convolver_in, (const float * const *) &impulse_ptr );
		ratio = 1.0 / ratio;
		lowpass_set_rate(impulse, FILTER_WIDTH, ratio);
		if (!r->convolver_out)
			r->convolver_out = convolver_create( (const float * const *) &impulse_ptr, FILTER_WIDTH, 1, 1, 1, 1 );
		else
			convolver_restage( r->convolver_out, (const float * const *) &impulse_ptr );
	}
}

void resampler_dup_inplace(void * r_out_, const void * r_in_)
{
	resampler * r_out = (resampler *) r_out_;
	const resampler * r_in = (const resampler *) r_in_;

	memcpy(r_out->history, r_in->history, sizeof(r_out->history));
	r_out->position = r_in->position;
	r_out->history_write = r_in->history_write;
	r_out->history_read = r_in->history_read;
	resampler_set_ratio(r_out, r_in->ratio);
}

void * resampler_dup(const void * r_in_)
{
	void * r_out_ = malloc(sizeof(resampler));
	if (r_out_)
	{
		resampler_dup_inplace(r_out_, r_in_);
	}
	return r_out_;
}

void resampler_clear(void * r_)
{
    resampler * r = (resampler *) r_;
    r->history_write = 0;
    r->history_read = 0;
    convolver_clear(r->convolver_in);
    convolver_clear(r->convolver_out);
}

void resampler_write(void * r_, float sample)
{
	resampler * r = (resampler *) r_;
	convolver_write(r->convolver_in, &sample);
}

void resampler_write_fixed(void * r_, int sample, int depth)
{
    float scale = 1.0f / (float)(1 << (depth - 1));
    float sf = (float)sample * scale;
    resampler_write(r_, sf);
}

int resampler_process(void * r_)
{
	resampler * r = (resampler *) r_;

	float sample;
	int ready = convolver_ready(r->convolver_in);

	float * history = r->history;
	int history_write = r->history_write;
	int history_read = r->history_read;
	int history_read_minus_1;
#ifdef CUBIC
	int history_read_minus_2;
	int history_read_minus_3;
#endif
	double ratio = r->ratio;
	double position = r->position;
#ifdef CUBIC
	history_read_minus_3 = (history_read + history_width - 3) % history_width;
	history_read_minus_2 = (history_read_minus_3 + 1) % history_width;
	history_read_minus_1 = (history_read_minus_2 + 1) % history_width;
#else
	history_read_minus_1 = (history_read + history_width - 1) % history_width;
#endif

    if ( history_read == history_write )
    {
        if ( !ready ) return 0;

        while ( ready-- )
        {
            convolver_read(r->convolver_in, &sample);
            history[history_write] = sample;
            history_write = (history_write + 1) % history_width;
        }
    }

    if ( convolver_ready(r->convolver_out) ) return 0;

	for (;;)
	{
#ifdef CUBIC
		float s0, s1, s2, s3;
		double A, B, C, D;

		if (position < 1.0)
		{
			s0 = history[history_read_minus_3];
			s1 = history[history_read_minus_2];
			s2 = history[history_read_minus_1];
			s3 = history[history_read];

			A = s3 - s2 - s0 + s1;
			B = s0 - s1 - A;
			C = s2 - s0;
			D = s1;
		}

		while (position < 1.0)
		{
			sample = A * position * position * position + B * position * position + C * position + D;

			convolver_write(r->convolver_out, &sample);
			position += ratio;

			if ((ready = convolver_ready(r->convolver_out))) break;
		}
#else
		float s0, s1, sd;

		if (position < 1.0)
		{
			s0 = history[history_read_minus_1];
			s1 = history[history_read];
			sd = s1 - s0;
		}

		while (position < 1.0)
		{
			sample = s0 + sd * position;

			convolver_write(r->convolver_out, &sample);
			position += ratio;

			if ((ready = convolver_ready(r->convolver_out))) break;
		}
#endif

		if (ready && position < 1.0) break;

#ifdef CUBIC
		history_read_minus_3 = history_read_minus_2;
		history_read_minus_2 = history_read_minus_1;
#endif
		history_read_minus_1 = history_read;
		history_read = (history_read + 1) % history_width;
		position -= 1.0;

		if (ready || history_read == history_write) break;
	}

	r->position = position;
	r->history_write = history_write;
	r->history_read = history_read;

	return 1;
}

int resampler_get_free_count(void * r_)
{
	resampler * r = (resampler *) r_;
	return convolver_ready(r->convolver_in) ? 0 : convolver_get_free_count(r->convolver_in);
}

int resampler_ready(void * r_)
{
	resampler * r = (resampler *) r_;
	return convolver_ready(r->convolver_out);
}

float resampler_read(void * r_)
{
	float sample;
	resampler * r = (resampler *) r_;
	convolver_read(r->convolver_out, &sample);
	return sample;
}

#ifdef MAIN
int main(void)
{
	unsigned char header[58];

	unsigned int input_rate;

	unsigned int output_rate = 44100;

	double ratio;

	FILE * f_in, * f_out;

	size_t length;

	size_t samples_in;
	size_t samples_out;
	size_t samples_to_generate;

	int eof;

	void * resampler = resampler_create();


	f_in = fopen("sweep96.wav", "rb");
	f_out = fopen("sweep441.wav", "wb");

	fread(header, 1, 58, f_in);
	fwrite(header, 1, 58, f_out);

	input_rate = ((unsigned int *)(header + 24))[0];

	ratio = (double)(input_rate) / (double)(output_rate);

	resampler_set_ratio(resampler, ratio);

	samples_in = 0;
	samples_out = 0;
	samples_to_generate = 0;

	eof = 0;

	while (!eof || samples_out < samples_to_generate)
	{
		float sample;

		int samples_ready = resampler_ready( resampler );

		if ( !samples_ready )
		{
			int to_fill = resampler_get_free_count( resampler );

			while ( to_fill-- )
			{
				if ( !eof )
				{
					if (fread(&sample, sizeof(sample), 1, f_in) == 1)
					{
						++samples_in;
					}
					else
					{
						samples_to_generate = (size_t)((double)samples_in / ratio + 0.5);
						eof = 1;
						sample = 0.0;
					}
				}
				else sample = 0.0;

				resampler_write( resampler, sample );
			}

			while ( resampler_process( resampler ) );

			samples_ready = resampler_ready( resampler );
		}

		while ( samples_ready-- )
		{
			sample = resampler_read( resampler );

			fwrite(&sample, sizeof(sample), 1, f_out);

			++samples_out;

			if (eof && samples_out == samples_to_generate) break;
		}
	}

	resampler_delete(resampler);

	fclose(f_in);

	length = ftell(f_out);

	((unsigned int *)(header + 24))[0] = output_rate;
	((unsigned int *)(header + 28))[0] = output_rate * sizeof(float);

	((unsigned int *)(header + 46))[0] = samples_out;

	((unsigned int *)(header + 4))[0] = length - 8;

	((unsigned int *)(header + 54))[0] = length - 58;

	fseek(f_out, 0, SEEK_SET);

	fwrite(header, 1, 58, f_out);

	fclose(f_out);

	return 0;
}
#endif

