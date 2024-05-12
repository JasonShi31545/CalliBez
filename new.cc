#include "setup.h"

using namespace std;

Point a,b,c;
vector<Point> points;

Point pp, pq;
Point p,q;
Point _p0, _p1, _p2;
float w;

Point t11,t12,t13;

const float SHIFT_X = W_WIDTH/2.0f;
const float SHIFT_Y = W_HEIGHT/2.0f;

const float SCALE_X = 2.0f;
const float SCALE_Y = 1.0f;

const float SHEAR_X = 2.0f;

const float ROTATE_ANGLE = 55.0f * (M_PI / 180.0f);



std::vector<Point> interpolatedPoints;
std::vector<Point> ipOutputs;




void Update(PixelGrid *grid, float t) {

    p = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, w, t), SCALE_X, SCALE_Y), ROTATE_ANGLE), SHIFT_X, SHIFT_Y);
    q = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, -w, t), SCALE_X, SCALE_Y), ROTATE_ANGLE), SHIFT_X, SHIFT_Y);


    float ti = BSRQC(_p0, _p1, _p2, w, t), tj = BSRQC(_p0, _p1, _p2, -w, t);

    printf("ti: %f, tj: %f\n", ti, tj);


    DrawLine(*grid, p, pp);
    DrawLine(*grid, q, pq);

    pp = p;
    pq = q;


    float tmp;
    Point curve, curve2;
    for (tmp = 0.0f; tmp <= 1.0f; tmp+=0.001f) {
        curve = BernsteinCubicSpline(a, points[1], points[2], b, tmp);
        curve2 = BernsteinCubicSpline(b, points[4], points[5], c, tmp);
        DrawPoint(*grid, curve);
        DrawPoint(*grid, curve2);
    }

    DrawCircle(*grid, a, 2.0f);
    DrawCircle(*grid, b, 2.0f);
    DrawCircle(*grid, c, 2.0f);

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

    float i,j;
    i = 0, j = 0;

    float ax,ay,bx,by,cx,cy;
    ax = 5.244;
    ay = 1.943;
    bx = 4.931;
    by = 1.499;
    cx = 4.600;
    cy = 1.050;


    // printf("Enter a: ");
    // scanf("%f%f", &ax, &ay);
    // printf("Enter b: ");
    // scanf("%f%f", &bx, &by);
    // printf("Enter c: ");
    // scanf("%f%f", &cx, &cy);

    // printf("Initial curvature: ");
    // scanf("%f", &i);
    // printf("Final curvature: ");
    // scanf("%f", &j);

    a = Point{ax,ay};
    b = Point{bx,by};
    c = Point{cx,cy};

    a = LinearScale(a, 100, 100);
    b = LinearScale(b, 100, 100);
    c = LinearScale(c, 100, 100);

    a = ShiftCoordinate(a, 100, 200);
    b = ShiftCoordinate(b, 100, 200);
    c = ShiftCoordinate(c, 100, 200);

    auto res = CubicSplineInterpolation(3, vector<Point>({a,b,c}), i, j);

    points = {a, res[0].first, res[0].second, b, res[1].first, res[1].second, c};
    SetupAndLoop(Update);
    return 0;
}
