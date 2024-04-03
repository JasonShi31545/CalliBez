#ifndef SETUP_H_
#define SETUP_H_

#include "bezier.h"

#define FPS 2000
#define FTT (1000 / FPS)

SDL_Window *initialize_window(const char *title);

SDL_Renderer *initialize_renderer(SDL_Window *window);

SDL_Surface *initialize_surface(const char *image_file_name);

SDL_Texture *initialize_texture(SDL_Surface *surface, SDL_Renderer *renderer);

int initialize_frame(SDL_Window **w, SDL_Renderer **r, SDL_Surface **s, SDL_Texture **t);

void destroy(SDL_Window *w, SDL_Renderer *r);

void SetupAndLoop(void (*calcAndUpdate)(PixelGrid *,float));

float TimeTransform(float t);

#endif // SETUP_H_
