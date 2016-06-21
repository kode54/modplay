#ifndef _RESAMPLER_H
#define _RESAMPLER_H

#ifdef __cplusplus
extern "C" {
#endif

void * resampler_create();
void resampler_delete(void *);

void * resampler_dup(const void *);
void resampler_dup_inplace(void *, const void *);
    
void resampler_clear(void *);

void resampler_set_ratio(void *, double ratio);

int resampler_get_free_count(void *);
void resampler_write(void *, float sample);
void resampler_write_fixed(void *, int sample, int depth);

int resampler_process(void *);
    
int resampler_ready(void *);
float resampler_read(void *);

#ifdef __cplusplus
}
#endif

#endif

