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


float t_;

int main(int argc, const char *argv[]) {
    using namespace std;

    t_ = 0.0f;

    PixelGrid *grid = new vector<vector<uint32_t>> ();
    grid->resize(W_WIDTH);
    for (size_t i = 0 ; i < W_WIDTH; i++) {
        (*grid)[i].resize(W_HEIGHT);
        for (size_t j = 0; j < W_HEIGHT; j++) {
            (*grid)[i][j] = (uint32_t)0;
        }
    }


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

        // UPDATE 1

        // int time_to_wait = FTT - (SDL_GetTicks() - last_frame_time);
        // if (time_to_wait > 0 && time_to_wait <= FTT) {
        //     SDL_Delay(time_to_wait);
        // }

        // this is cool
        // float delta_time = (SDL_GetTicks() - last_frame_time);
        last_frame_time = SDL_GetTicks();

        t_ = last_frame_time - 280;
        t_ = TimeTransform(t_);
        // std::cerr << "Delta time: " << delta_time << std::endl;
        // std::cerr << "Time: " << SDL_GetTicks() << std::endl;


        // Update

        // Point a,b,c,d;
        // a = Point(100,100);
        // d = Point(500, 100);
        // b = Point(180, 300);
        // c = Point(400, 20);

        // DrawPoint((*grid), a);
        // DrawPoint((*grid), b);
        // DrawPoint((*grid), c);
        // DrawPoint((*grid), d);


        // Point f = BersteinCubicSpline(a,b,c,d,t_);
        // Vector2 v = BersteinCubicVelocity(a, b, c, d, t_);


        // Point w1, w2, w3, w4;
        // w1 = Point(0,0);
        // w2 = Point(100,100);
        // w3 = Point(200,100);
        // w4 = Point(300,0);

        // Point w = BersteinCubicSpline(w1, w2, w3, w4, t_);

        // // DrawPoint(r, f);
        // DrawWidth((*grid), f, v, (1.0f/20.0f)*w.y);

        // float x = 25 / cosf(30* t_) + 300;
        // float y = 2 * tanf(30*t_) + 300;
        // float x = 100 * cosf(20*t_*2*M_PI);
        // float y = 100 * sinf(20*t_*2*M_PI);
        // Vector2 v(x,y);
        // v = RotateV2(v, (M_PI/180)*35);

        // Point p = v.toPoint();
        // p.x += 300;
        // p.y += 300;
        // DrawPoint((*grid), p);

        Point p0 = Point(100,100);
        Point p1 = Point(150, 100);
        Point p2 = Point(150,50);
        vector<float> weights = {1.0f, sqrtf(2.0)/2.0f, 1.0f};
        vector<Point> points = {p0, p1, p2};

        Point res = BezierCurveRationalWeighted(2, points, weights, 25*t_);
        DrawPoint((*grid), res);



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

                    SDL_SetRenderDrawColor(r, red, g, b, a);
                    SDL_RenderDrawPoint(r, i, j);
                }
            }
        }
        SDL_SetRenderDrawColor(r, 0, 0, 0, 0);

        SDL_RenderPresent(r);

    }

    destroy(w,r);
    delete grid;

    return 0;
}
