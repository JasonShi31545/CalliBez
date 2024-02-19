#ifndef BEZIER_H
#define BEZIER_H

#include "header.h"

#define W_WIDTH 1200
#define W_HEIGHT 700
#define FPS 140
#define FTT (1000 / FPS)

// typedef struct {
//     float x;
//     float y;
// } Point;

// typedef struct {
//     Point p0;
//     Point p1;
//     float thickness;
// } Line;

// float LineGetLength(Line l) {
//     return sqrtf((l.p0.x - l.p1.x) * (l.p0.x - l.p1.x) + (l.p0.y - l.p1.y) * (l.p0.y - l.p1.y));
// }


// Point lerp(Point p0, Point p1, float t) {
//     Point res{};
//     res.x = (1-t)*p0.x + t*p1.x;
//     res.y = (1-t)*p0.y + t*p1.y;
//     return res;
// }


class Point {
    public:
        float x,y;
        bool isjoin;
        Point(): x(0), y(0), isjoin(false) {};
        Point(float a, float b): x(a), y(b), isjoin(false) {};

        Point operator+(Point p) {
            return Point(x+p.x, y+p.y);
        }
        Point operator-(Point p) {
            return Point(x-p.x, y-p.y);
        }
        Point operator*(Point p) {
            return Point(x*p.x, y-p.y);
        }
};

class Vector2 {
    public:
        float x, y;
        Vector2(): x(0), y(0) {};
        Vector2(float a, float b) : x(a), y(b) {};
        Vector2(Point &&p) {
            x = p.x;
            y = p.y;
        }
        Vector2(Point &p) {
            x = p.x;
            y = p.y;
        }
        Vector2(Point p) {
            x = p.x;
            y = p.y;
        }
        Vector2 operator+(Vector2 v) {
            return Vector2(x+v.x, y+v.y);
        }
        Vector2 operator-(Vector2 v) {
            return Vector2(x-v.x, y-v.y);
        }
        float magnitude() {
            return sqrtf(x*x + y*y);
        }
};



float TimeTransform(float t) {
//    return (1 - expf(-0.05f * t));
    return (1 - expf(-0.0005f * t));
}


Point lerp(Point s, Point e, float t) {
    assert((0.0f <= t) && (t <= 1.0f));
    float x, y;
    x = (1-t) * s.x + t * e.x;
    y = (1-t) * s.y + t * e.y;
    return Point(x,y);
}

void DrawPoint(SDL_Renderer *r, Point p) {
    SDL_RenderDrawPointF(r, p.x, W_HEIGHT - p.y);
}

void DrawLine(SDL_Renderer *r, Point s, Point e) {
    SDL_RenderDrawLineF(r, s.x, W_HEIGHT - s.y, e.x, W_HEIGHT - e.y);
}




Point BersteinCubicSpline(Point p0, Point p1, Point p2, Point p3, float t) {
    float x,y;
    x = p0.x * (1 - 3*t + 3*t*t - t*t*t) + p1.x * (3*t - 6*t*t + 3*t*t*t) + p2.x * (3*t*t - 3*t*t*t) + p3.x * (t*t*t);
    y = p0.y * (1 - 3*t + 3*t*t - t*t*t) + p1.y * (3*t - 6*t*t + 3*t*t*t) + p2.y * (3*t*t - 3*t*t*t) + p3.y * (t*t*t);
    return Point(x,y);
}

Vector2 BersteinCubicVelocity(Point p0, Point p1, Point p2, Point p3, float t) {
    float x,y;
    x = p0.x * (3 + 6*t - 3*t*t) + p1.x * (3 - 12*t + 9*t*t) + p2.x * (6*t - 9*t*t) + p3.x * (3*t*t);
    y = p0.y * (3 + 6*t - 3*t*t) + p1.y * (3 - 12*t + 9*t*t) + p2.y * (6*t - 9*t*t) + p3.y * (3*t*t);
    return Vector2(x,y);
}


Vector2 BersteinCubicAcceleration(Point p0, Point p1, Point p2, Point p3, float t) {
    float x,y;
    x = p0.x * (6 - 6*t) + p1.x * (-12 + 18*t) + p2.x * (6-18*t) + p3.x * (6*t);
    y = p0.y * (6 - 6*t) + p1.y * (-12 + 18*t) + p2.y * (6-18*t) + p3.y * (6*t);
    return Vector2(x,y);
}

float Curvature(Point p0, Point p1, Point p2, Point p3, float t) {
    Vector2 v = BersteinCubicVelocity(p0, p1, p2, p3, t);
    Vector2 a = BersteinCubicAcceleration(p0, p1, p2, p3, t);
    float det = (v.x * a.y) - (a.x * v.y); // determinant
    float m = v.magnitude();
    return det / (m*m*m);
}

float Radius(Point p0, Point p1, Point p2, Point p3, float t) {
    return 1.0f / Curvature(p0, p1, p2, p3, t);
}

unsigned int factorial(unsigned int n) {
    for (int i = n-1; i >= 1; i--) {
        n *= i;
    }
    return n;
}

unsigned int choose(unsigned int n, unsigned int k) {
    // n choose k
    // n! / k!(n-k)!
    // => (n) * (n-1) * (n-2) * ... * (n-k+1) / (1 * 2 * 3 ... * k)
    unsigned int res = 1;
    for (unsigned int i = n; i > n-k; i--) {
        res *= i;
    }
    return res / factorial(k);
}

float BersteinCoefficient(unsigned int n, unsigned int i, float t) {
    return choose(n,i) * expf(logf(t) * (float)i) * expf(logf(1-t) * (float)(n-i));
}

Point BezierCurveUnweightedN(unsigned int n, std::vector<Point> points, float t) {
    assert(points.size() == n+1);
    float resx = 0, resy = 0;
    for (int i = 0; i <= n; i++) {
        resx += BersteinCoefficient(n, i, t) * points[i].x;
        resy += BersteinCoefficient(n, i, t) * points[i].y;
    }
    return Point(resx, resy);
}

Point BezierCurveWeightedN(unsigned int n, std::vector<Point> points, std::vector<float> weights, float t) {
    assert(points.size() == n+1);
    assert(weights.size() == n+1);
    float resx = 0, resy = 0;
    for (int i = 0; i <= n; i++) {
        resx += BersteinCoefficient(n, i, t) * points[i].x * weights[i];
        resy += BersteinCoefficient(n, i, t) * points[i].y * weights[i];
    }
    return Point(resx, resy);
}

Point BezierCurveRationalWeighted(unsigned int n, std::vector<Point> points, std::vector<float> weights, float t) {
    assert(points.size() == n+1);
    assert(weights.size() == n+1);
    float resx = 0, resy = 0;
    float den = 0;
    for (int i = 0; i <= n; i++) {
        den += BersteinCoefficient(n, i, t) * weights[i];
        resx += den * points[i].x;
        resy += den * points[i].y;
    }
    resx /= den;
    resy /= den;
    return Point(resx, resy);
}


#endif // BEZIER_H
