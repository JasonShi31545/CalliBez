#include "setup.h"

Point pp, pq;
Point p,q;
Point _p0, _p1, _p2;
float w;

Point a,b,c;
Point t11,t12,t13;

const float SHIFT_X = W_WIDTH/2.0f;
const float SHIFT_Y = W_HEIGHT/2.0f;

const float SCALE_X = 2.0f;
const float SCALE_Y = 1.0f;

const float SHEAR_X = 2.0f;

const float ROTATE_ANGLE = M_PI_4f;



std::vector<Point> interpolatedPoints;
std::vector<Point> ipOutputs;



void Update(PixelGrid *grid, float t) {

    p = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, w, t), SCALE_X, SCALE_Y), ROTATE_ANGLE), SHIFT_X, SHIFT_Y);
    q = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, -w, t), SCALE_X, SCALE_Y), ROTATE_ANGLE), SHIFT_X, SHIFT_Y);


    DrawLine(*grid, p, pp);
    DrawLine(*grid, q, pq);

    pp = p;
    pq = q;


    for (size_t i = 0; i < interpolatedPoints.size(); i++) {
        DrawCircle(*grid, ShiftCoordinate(interpolatedPoints[i], 50.0f, 50.0f), 30.0f);
    }



    for (size_t i = 0; i < interpolatedPoints.size() - 1; i++) {
        Point p0, a, b, p2;
        p0 = interpolatedPoints[i];
        p2 = interpolatedPoints[i+1];
        a = ipOutputs[i * 2 + 0];
        b = ipOutputs[i * 2 + 1];
        Point res = BernsteinCubicSpline(p0, a, b, p2, t);
        res = ShiftCoordinate(res, 50.0f, 50.0f);
        DrawPoint(*grid, res);

    }

}



int main(int argc, const char *argv[]) {

    _p0 = Point{100, 100};
    float alpha = 1.05f;
    std::tuple<HomogeneousPoint, Point, float> circ = ConstructArc(_p0, alpha);
    HomogeneousPoint p1 = std::get<0>(circ);
    _p1 = Point{p1.x, p1.y};
    w = p1.w;

    _p2 = std::get<1>(circ);

    pp = ShiftCoordinate(LinearRotate(LinearScale(_p0, SCALE_X, SCALE_Y), ROTATE_ANGLE), SHIFT_X, SHIFT_Y);
    pq = pp;

    interpolatedPoints = {
        Point{200.0f , 300.0f},
        Point{300.0f, 150.0f},
        Point{700.0f, 400.0f},
        Point{800.0f, 75.0f},
        Point{900.0f, 400.0f}
    };


    std::vector<std::pair<Point, Point>> interpolates = CubicSplineInterpolation(interpolatedPoints.size(), interpolatedPoints, 0.0f, 0.0f);

    ipOutputs.resize(2*(interpolatedPoints.size() - 1));
    for (size_t i = 0; i < interpolatedPoints.size() - 1; i++) {
        ipOutputs[i*2 + 0] = interpolates[i].first;
        ipOutputs[i*2 + 1] = interpolates[i].second;
    }

    // for (size_t i = 0; i < interpolatedPoints.size() - 1; i++) {
    //     std::cerr << "a (" << i << "): " << "x: " << ipOutputs[i*2 + 0].x << " y: " << ipOutputs[i*2 + 0].y << std::endl;
    //     std::cerr << "b (" << i << "): " << "x: " << ipOutputs[i*2 + 1].x << " y: " << ipOutputs[i*2 + 1].y << std::endl;
    // }


    a = Point(5.244f, 1.943);
    b = Point(4.931f, 1.499f);
    c = Point(4.6f, 1.05f);


    t11 = Point(0.0f, 0.0f);
    t12 = Point(0.5f, (6.0f/16.0f));
    t13 = Point(1.0f, 0.0f);

    a = LinearScale(a, 10, 10);
    b = LinearScale(b, 10, 10);
    c = LinearScale(c, 10, 10);

    t11 = LinearScale(t11, 10, 10);
    t12 = LinearScale(t12, 10, 10);
    t13 = LinearScale(t13, 10, 10);

    auto rline = CubicSplineInterpolation(3, {a,b,c}, 3, 3);





    SetupAndLoop(Update);
    return 0;
}
