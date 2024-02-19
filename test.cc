#include "bezier.h"

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


float t_;

int main(int argc, const char *argv[]) {
    using namespace std;

    t_ = 0.0f;
    vector<Point> *points_plotted = new vector<Point>();
    unique_ptr<vector<Point>> pp(points_plotted);
    points_plotted = nullptr;

    int last_frame_time = SDL_GetTicks();

    srand(time(NULL));

    SDL_Window *w;
    SDL_Renderer *r;
    SDL_Surface *s;
    SDL_Texture *t;

    assert(initialize_frame(&w,&r,&s,&t) == 0);
    bool running = true;


    SDL_Event event;
    while (running) {
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

        SDL_SetRenderDrawColor(r, 0, 0, 0, 0);

        SDL_RenderClear(r);

        // UPDATE

        int time_to_wait = FTT - (SDL_GetTicks() - last_frame_time);
        if (time_to_wait > 0 && time_to_wait <= FTT) {
            SDL_Delay(time_to_wait);
        }
        // this is cool
        float delta_time = (SDL_GetTicks() - last_frame_time);
        last_frame_time = SDL_GetTicks();

        t_ = last_frame_time - 280;
        t_ = TimeTransform(t_);
        // std::cerr << "Delta time: " << delta_time << std::endl;
        // std::cerr << "Time: " << SDL_GetTicks() << std::endl;

        // Render

        SDL_SetRenderDrawColor(r, 255, 200, 50, 255);


        // printf("%f\n", t_);

        Point a,b,c,d;
        a = Point(100,100);
        d = Point(500, 100);
        b = Point(180, 300);
        c = Point(400, 20);

        DrawPoint(r, a);
        DrawPoint(r, b);
        DrawPoint(r, c);
        DrawPoint(r, d);


        Point f = BersteinCubicSpline(a,b,c,d,t_);

        for (size_t i = 0; i < (*pp).size(); i++) {
            DrawPoint(r, (*pp)[i]);
        }
        DrawPoint(r, f);
        (*pp).push_back(f);


        SDL_RenderPresent(r);

    }

    destroy(w,r);

    return 0;
}
