#ifndef _SIMPLE_CONVOLVER_H_
#define _SIMPLE_CONVOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

void * convolver_create(const float * const* impulses, int impulse_size, int input_channels, int output_channels, int impulse_count, int channels_per_impulse);
void convolver_restage(void *, const float * const* impulses);
void convolver_delete(void *);
void convolver_clear(void *);
int convolver_get_free_count(void *);
void convolver_write(void *, const float *);
int convolver_ready(void *);
void convolver_read(void *, float *);

#ifdef __cplusplus
}
#endif

#endif
