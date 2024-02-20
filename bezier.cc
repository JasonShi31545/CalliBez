#include "bezier.h"

float Vector2::magnitude() {
    return sqrtf(x*x + y*y);
}

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


void DrawCircle(SDL_Renderer *r, Point o, float radius) {
    float x,y;
    for (float t = 0.0f; t <= 2*M_PI; t += 0.01f) {
        x = radius*cosf(t) + o.x;
        y = radius*sinf(t) + o.y;
        DrawPoint(r, Point(x,y));
    }
}

void DrawWidth(SDL_Renderer *r, Point o, Vector2 v, float w) {
    // v for velocity vector
    // w for width
    Vector2 nn = NormalizeV2(NormalV2(v)); // normalised normal
    const float epsilon = 0.0007f; // width resolution
    for (float ae = 0.0f; ae <= w/2; ae += epsilon) {
        // ae: accumulated epsilon
        // w/2 as symmetrical
        float m = nn.y / nn.x;
        auto func = [&m, &o](float x)  {
            return m*x - o.x * m + o.y;
        };
        float x1 = o.x - ae;
        float x2 = o.x + ae;
        float y1 = func(x1);
        float y2 = func(x2);
        DrawPoint(r, Point(x1, y1));
        DrawPoint(r, Point(x2, y2));
    }
}

Point BersteinCubicSpline(Point p0, Point p1, Point p2, Point p3, float t) {
    float x,y;
    float s,c;
    s = t*t;
    c = t*t*t;
    x = p0.x * (1 - 3*t + 3*s - c) + p1.x * (3*t - 6*s + 3*c) + p2.x * (3*s - 3*c) + p3.x * (c);
    y = p0.y * (1 - 3*t + 3*s - c) + p1.y * (3*t - 6*s + 3*c) + p2.y * (3*s - 3*c) + p3.y * (c);
    return Point(x,y);
}

Vector2 BersteinCubicVelocity(Point p0, Point p1, Point p2, Point p3, float t) {
    float x,y;
    float s;
    s = t*t;
    x = p0.x * (3 + 6*t - 3*s) + p1.x * (3 - 12*t + 9*s) + p2.x * (6*t - 9*s) + p3.x * (3*s);
    y = p0.y * (3 + 6*t - 3*s) + p1.y * (3 - 12*t + 9*s) + p2.y * (6*t - 9*s) + p3.y * (3*s);
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
    if (n == 0) return 1;
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
        float bc = BersteinCoefficient(n, i, t) * weights[i];
        den += bc;
        resx += den * points[i].x;
        resy += den * points[i].y;
    }
    resx /= den;
    resy /= den;
    return Point(resx, resy);
}

// Vectors

Vector2 NormalizeV2(Vector2 v) {
    float mag = v.magnitude();
    return Vector2(v.x / mag, v.y / mag);
}


float AngleV2(Vector2 v) {
    return atanf(v.y / v.x);
}

Vector2 RotateV2(Vector2 v, float angle) {
    // Angle in rad
    float x = v.x * cosf(angle) - v.y * sinf(angle);
    float y = v.x * sinf(angle) - v.y * cosf(angle);
    return Vector2(x,y);
}

Vector2 NormalV2(Vector2 v) {
    return Vector2(-v.y, v.x);
}
