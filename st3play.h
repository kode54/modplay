#ifndef _ST3PLAY_H_
#define _ST3PLAY_H_

#ifdef __cplusplus
extern "C" {
#endif

void * st3play_Alloc(uint32_t outputFreq, int8_t interpolation);
void st3play_Free(void *);

int8_t st3play_LoadModule(void *, const uint8_t *module, size_t size);
void st3play_PlaySong(void *);

void st3play_RenderFixed(void *_p, int32_t *buffer, int32_t count);
void st3play_Render16(void *_p, int16_t *buffer, int32_t count);
    
#ifdef __cplusplus
}
#endif

#endif
