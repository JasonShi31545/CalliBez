#include "bezier.h"

#define FPS 2000
#define FTT (1000 / FPS)


SDL_Window *initialize_window(const char *title) {
    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        fprintf(stderr, "Error initialising SDL\n");
        return NULL;
    }
    SDL_Window *window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, W_WIDTH, W_HEIGHT, SDL_WINDOW_ALLOW_HIGHDPI);
    return window;
}

SDL_Renderer *initialize_renderer(SDL_Window *window) {
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    return renderer;
}

SDL_Surface *initialize_surface(const char *image_file_name) {
    SDL_Surface *image = SDL_LoadBMP(image_file_name);
    return image;
}

SDL_Texture *initialize_texture(SDL_Surface *surface, SDL_Renderer *renderer) {
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);
    return texture;
}

int initialize_frame(SDL_Window **w, SDL_Renderer **r, SDL_Surface **s, SDL_Texture **t) {
    SDL_Window *window = initialize_window("");
    if (!window) {
        return -1;
    }
    SDL_Renderer *renderer = initialize_renderer(window);

    *w = window;
    *r = renderer;
    *s = initialize_surface("");
    *t = initialize_texture(*s,*r);

    return 0;
}

void destroy(SDL_Window *w, SDL_Renderer *r) {
    SDL_DestroyRenderer(r);
    SDL_DestroyWindow(w);
    SDL_Quit();
}


float t;

int main(int argc, const char *argv[]) {
    using namespace std;

    // Time Setup
    t = 0.0f;
    int last_frame_time = SDL_GetTicks();
    srand(time(NULL));

    // Pixel Canvas Setup
    PixelGrid *grid = new vector<vector<uint32_t>> ();
    grid->resize(W_WIDTH);
    for (size_t i = 0 ; i < W_WIDTH; i++) {
        (*grid)[i].resize(W_HEIGHT);
        for (size_t j = 0; j < W_HEIGHT; j++) {
            (*grid)[i][j] = (uint32_t)0;
        }
    }

    // SDL Setup
    SDL_Window *win;
    SDL_Renderer *ren;
    SDL_Surface *sur;
    SDL_Texture *tex;

    assert(initialize_frame(&win,&ren,&sur,&tex) == 0);

    // Main Event Loop
    bool running = true;
    SDL_Event event;
    while (running) {
        // Poll Events and interrupts
        SDL_PollEvent(&event);
        switch (event.type) {
            case SDL_QUIT:
                running = false;
                break;
            case SDL_KEYDOWN:
                if (event.key.keysym.sym == SDLK_ESCAPE || event.key.keysym.sym == SDLK_q) {
                    running = false;
                }
                break;
            default:
                break;
        }

        // Setup
        SDL_SetRenderDrawColor(ren, 0, 0, 0, 0);
        SDL_RenderClear(ren);


        // Timing

        last_frame_time = SDL_GetTicks();

        t = last_frame_time - 280;
        t = TimeTransform(t);

        // Calculate & Update


        // Render

        for (int i = 0; i < W_WIDTH; i++) {
            for (int j = 0; j < W_HEIGHT; j++) {
                uint32_t pixel_val = (*grid)[i][j];
                if (pixel_val != (uint32_t)0) {
                    // top byte is r, then g, then b, then the bottom is a

                    Uint8 red,g,b,a;
                    red = (Uint8)((pixel_val & 0xFF000000) >> 24U);
                    g = (Uint8)((pixel_val & 0x00FF0000) >> 16U);
                    b = (Uint8)((pixel_val & 0x0000FF00) >> 8U);
                    a = (Uint8)((pixel_val & 0x000000FF) >> 0U);

                    SDL_SetRenderDrawColor(ren, red, g, b, a);
                    SDL_RenderDrawPoint(ren, i, j);
                }
            }
        }
        SDL_SetRenderDrawColor(ren, 0, 0, 0, 0);

        SDL_RenderPresent(ren);

    }

    // Clean up

    destroy(win,ren);
    delete grid;

    return 0;
}
