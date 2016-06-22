#include <stdlib.h>
#include <string.h>

#include "resampler.h"

#include "fir_resample.h"

enum { resampler_buffer_size = 1024 };

static int initialized = 0;

void resampler_init(void)
{
    if (!initialized)
    {
        initialized = 1;
        fir_init();
    }
}

typedef struct resampler
{
    int write_pos, write_filled;
    int read_pos, read_filled;
    int delay;
    float phase;
    float phase_inc;
    float buffer_in[resampler_buffer_size * 2];
    float buffer_out[resampler_buffer_size];
    
    /* temporary */
    float * out;
} resampler;

void * resampler_create(void)
{
    resampler * r = ( resampler * ) malloc( sizeof(resampler) );
    if ( !r ) return 0;

    r->write_pos = 0;
    r->write_filled = 0;
    r->read_pos = 0;
    r->read_filled = 0;
    r->delay = 0;
    r->phase = 0;
    r->phase_inc = 0;
    memset( r->buffer_in, 0, sizeof(r->buffer_in) );
    memset( r->buffer_out, 0, sizeof(r->buffer_out) );

    return r;
}

void resampler_delete(void * _r)
{
    free( _r );
}

void * resampler_dup(const void * _r)
{
    void * r_out = malloc( sizeof(resampler) );
    if ( !r_out ) return 0;

    resampler_dup_inplace(r_out, _r);

    return r_out;
}

void resampler_dup_inplace(void *_d, const void *_s)
{
    const resampler * r_in = ( const resampler * ) _s;
    resampler * r_out = ( resampler * ) _d;

    r_out->write_pos = r_in->write_pos;
    r_out->write_filled = r_in->write_filled;
    r_out->read_pos = r_in->read_pos;
    r_out->read_filled = r_in->read_filled;
    r_out->delay = r_in->delay;
    r_out->phase = r_in->phase;
    r_out->phase_inc = r_in->phase_inc;
    memcpy( r_out->buffer_in, r_in->buffer_in, sizeof(r_in->buffer_in) );
    memcpy( r_out->buffer_out, r_in->buffer_out, sizeof(r_in->buffer_out) );
}

void resampler_set_quality(void *_r, int quality)
{
    (void)_r;
    (void)quality;
}

int resampler_get_free_count(void *_r)
{
    resampler * r = ( resampler * ) _r;
    return resampler_buffer_size - r->write_filled;
}

int resampler_get_padding_size()
{
    return 256;
}

static int resampler_min_filled(resampler *r)
{
    return fir_latency(r->phase_inc);
}

int resampler_ready(void *_r)
{
    resampler * r = ( resampler * ) _r;
    return r->write_filled > resampler_min_filled(r);
}

void resampler_clear(void *_r)
{
    resampler * r = ( resampler * ) _r;
    r->write_pos = 0;
    r->write_filled = 0;
    r->read_pos = 0;
    r->read_filled = 0;
    r->delay = 0;
    r->phase = 0;
}

void resampler_set_rate(void *_r, double new_factor)
{
    resampler * r = ( resampler * ) _r;
    r->phase_inc = new_factor;
}

void resampler_write_sample(void *_r, short s)
{
    resampler * r = ( resampler * ) _r;
    
    if ( !r->delay )
    {
        int delay = fir_latency( r->phase_inc );
        r->write_pos = delay;
        r->write_filled = delay;
        memset( r->buffer_in, 0, delay * sizeof(*r->buffer_in) );
        r->delay = 1;
    }

    if ( r->write_filled < resampler_buffer_size )
    {
        double s32 = s;
        s32 *= 256.0;

        r->buffer_in[ r->write_pos ] = s32;
        r->buffer_in[ r->write_pos + resampler_buffer_size ] = s32;

        ++r->write_filled;

        r->write_pos = ( r->write_pos + 1 ) % resampler_buffer_size;
    }
}

void resampler_write_sample_fixed(void *_r, int s, unsigned char depth)
{
    resampler * r = ( resampler * ) _r;
    
    if ( !r->delay )
    {
        int delay = fir_latency( r->phase_inc );
        r->write_pos = delay;
        r->write_filled = delay;
        memset( r->buffer_in, 0, delay * sizeof(*r->buffer_in) );
        r->delay = 1;
    }
    
    if ( r->write_filled < resampler_buffer_size )
    {
        double s32 = s;
        s32 /= (double)(1 << (depth - 1));
        s32 += 1e-25;
        
        r->buffer_in[ r->write_pos ] = s32;
        r->buffer_in[ r->write_pos + resampler_buffer_size ] = s32;
        
        ++r->write_filled;
        
        r->write_pos = ( r->write_pos + 1 ) % resampler_buffer_size;
    }
}

static float resampler_getter(void * r_, size_t channel, size_t offset)
{
    resampler * r = (resampler *) r_;
    return r->buffer_in[resampler_buffer_size + r->write_pos - r->write_filled + offset];
}

static void resampler_putter(void * r_, size_t channel, size_t offset, float sample)
{
    resampler * r = (resampler *) r_;
    r->out[offset] = sample;
}

static int resampler_run(resampler * r, float ** out_, float * out_end)
{
    int in_size = r->write_filled;
    float * out = *out_;
    int out_size = out_end - out;
    int required_input = fir_required_input(out_size, r->phase_inc, r->phase);
    int in_used;
    if (required_input > in_size)
    {
        int new_out_size = fir_required_output( in_size, r->phase_inc, r->phase );
        if (new_out_size < out_size)
            out_size = new_out_size;
        if (!out_size) return 0;
    }
    r->out = out;
    in_used = fir_resample(resampler_getter, resampler_putter, r, 1, out_size, r->phase_inc, &r->phase);
    *out_ += out_size;
    return in_used;
}

static void resampler_fill(resampler * r)
{
    int min_filled = resampler_min_filled(r);
    while ( r->write_filled > min_filled &&
            r->read_filled < resampler_buffer_size )
    {
        int in_used;
        int write_pos = ( r->read_pos + r->read_filled ) % resampler_buffer_size;
        int write_size = resampler_buffer_size - write_pos;
        float * out = r->buffer_out + write_pos;
        if ( write_size > ( resampler_buffer_size - r->read_filled ) )
            write_size = resampler_buffer_size - r->read_filled;
        in_used = resampler_run( r, &out, out + write_size );
        r->write_filled -= in_used;
        r->read_filled += out - r->buffer_out - write_pos;
        if (!in_used) break;
    }
}

static void resampler_fill_and_remove_delay(resampler * r)
{
    resampler_fill( r );
}

int resampler_get_sample_count(void *_r)
{
    resampler * r = ( resampler * ) _r;
    if ( r->read_filled < 1 )
        resampler_fill_and_remove_delay( r );
    return r->read_filled;
}

int resampler_get_sample(void *_r)
{
    resampler * r = ( resampler * ) _r;
    if ( r->read_filled < 1 && r->phase_inc)
        resampler_fill_and_remove_delay( r );
    if ( r->read_filled < 1 )
        return 0;
    return (int)r->buffer_out[ r->read_pos ];
}

float resampler_get_sample_float(void *_r)
{
    resampler * r = ( resampler * ) _r;
    if ( r->read_filled < 1 && r->phase_inc)
        resampler_fill_and_remove_delay( r );
    if ( r->read_filled < 1 )
        return 0;
    return r->buffer_out[ r->read_pos ];
}

void resampler_remove_sample(void *_r, int decay)
{
    resampler * r = ( resampler * ) _r;
    if ( r->read_filled > 0 )
    {
        --r->read_filled;
        r->read_pos = ( r->read_pos + 1 ) % resampler_buffer_size;
    }
}
