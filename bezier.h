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
        float magnitude();
};


float TimeTransform(float t);
Point lerp(Point s, Point e, float t);
void DrawPoint(SDL_Renderer *r, Point p);
void DrawLine(SDL_Renderer *r, Point s, Point e);
void DrawCircle(SDL_Renderer *r, Point o, float r);
void DrawWidth(SDL_Renderer *r, Point o, Vector2 v);
Point BersteinCubicSpline(Point p0, Point p1, Point p2, Point p3, float t);
Vector2 BersteinCubicVelocity(Point p0, Point p1, Point p2, Point p3, float t);
Vector2 BersteinCubicAcceleration(Point p0, Point p1, Point p2, Point p3, float t);
float Curvature(Point p0, Point p1, Point p2, Point p3, float t);
float Radius(Point p0, Point p1, Point p2, Point p3, float t);
unsigned int choose(unsigned int n, unsigned int k);
float BersteinCoefficient(unsigned int n, unsigned int i, float t);
Point BezierCurveUnweightedN(unsigned int n, std::vector<Point> points, float t);
Point BezierCurveWeightedN(unsigned int n, std::vector<Point> points, std::vector<float> weights, float t);
Point BezierCurveRationalWeighted(unsigned int n, std::vector<Point> points, std::vector<float> weights, float t);

Vector2 NormalizeV2(Vector2 v);
float AngleV2(Vector2 v);
Vector2 RotateV2(Vector2 v, float angle);
Vector2 NormalV2(Vector2 v);

unsigned int factorial(unsigned int n);

#endif // BEZIER_H
