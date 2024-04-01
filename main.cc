#include "setup.h"


void Calculate(PixelGrid *grid, float t) {

    Point _p0 = Point{100, 100};
    float alpha = 1.05f;
    std::tuple<HomogeneousPoint, Point, float> circ = ConstructArc(_p0, alpha);
    HomogeneousPoint p1 = std::get<0>(circ);
    Point _p1 = Point{p1.x, p1.y};
    float w = p1.w;

    Point _p2 = std::get<1>(circ);

    std::vector<Point> points = {_p0, _p1, _p2};
    std::vector<float> weights1 = {1.0f, w, 1.0f};
    std::vector<float> weights2 = {1.0f, -w, 1.0f};

    Point p = BezierCurveRationalWeighted(2, points, weights1, t);
    p = LinearScale(p, 2.0f, 1.0f);
    p = ShiftCoordinate(p, 300, 300);

    Point q = BezierCurveRationalWeighted(2, points, weights2, t);
    q = LinearScale(q, 2.0f, 1.0f);
    q = ShiftCoordinate(q, 300, 300);

    DrawPoint(*grid, p);
    DrawPoint(*grid, q);
}



int main(int argc, const char *argv[]) {
    SetupAndLoop(Calculate);
    return 0;
}
