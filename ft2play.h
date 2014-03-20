#ifndef _FT2PLAY_H_
#define _FT2PLAY_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
    
void * ft2play_Alloc(uint32_t _samplingFrequency, int8_t interpolation);
void ft2play_Free(void *);

int8_t ft2play_LoadModule(void *, const int8_t *buffer, size_t size);

void ft2play_PlaySong(void *, int32_t startOrder);
    
void ft2play_RenderFixed(void *, int32_t *buffer, int32_t count);
void ft2play_Render16(void *, int16_t *buffer, int32_t count);
    
uint32_t ft2play_GetLoopCount(void *);

typedef struct
{
    uint16_t order;
    uint16_t pattern;
    uint16_t row;
    uint16_t speed;
    uint16_t tempo;
    uint8_t channels_playing;
} ft2_info;
    
void ft2play_GetInfo(void *, ft2_info *);
    
#ifdef __cplusplus
}
#endif

#endif
