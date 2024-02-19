#ifndef BEZIER_H
#define BEZIER_H

#include "header.h"

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
    private:
        float x,y;
        bool isjoin;
    public:
        Point(): x(0), y(0), isjoin(false) {};
        Point(float a, float b): x(a), y(b), isjoin(false) {};
        Point(Point &&p) { x = p.x; y = p.y; isjoin = p.isjoin; };

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




#endif // BEZIER_H
