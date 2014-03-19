#ifndef _ST3PLAY_H_
#define _ST3PLAY_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void * st3play_Alloc(uint32_t outputFreq, int8_t interpolation);
void st3play_Free(void *);

int8_t st3play_LoadModule(void *, const uint8_t *module, size_t size);
void st3play_PlaySong(void *, int16_t startOrder);
    
int32_t st3play_GetLoopCount(void *);

void st3play_RenderFixed(void *, int32_t *buffer, int32_t count);
void st3play_Render16(void *, int16_t *buffer, int32_t count);
    
typedef struct
{
    uint16_t order;
    uint16_t pattern;
    uint16_t row;
    uint8_t speed;
    uint8_t tempo;
    uint8_t channels_playing;
} st3_info;
    
void st3play_GetInfo(void *, st3_info *);
    
#ifdef __cplusplus
}
#endif

#endif
