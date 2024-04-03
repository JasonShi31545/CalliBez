#include "setup.h"

Point pp, pq, pr, ps, pu, pv, py, pz;
Point p, q, r, s, u, v, y, z;
Point _p0, _p1, _p2;
float w;

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

    r = ShiftCoordinate(LinearScale(BSRQS(_p0, _p1, _p2, w, t), SCALE_X, SCALE_Y), SHIFT_X, SHIFT_Y);
    s = ShiftCoordinate(LinearScale(BSRQS(_p0, _p1, _p2, -w, t), SCALE_X, SCALE_Y), SHIFT_X, SHIFT_Y);

    u = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, w, t), SCALE_X, SCALE_Y), -ROTATE_ANGLE), SHIFT_X, SHIFT_Y);
    v = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, -w, t), SCALE_X, SCALE_Y), -ROTATE_ANGLE), SHIFT_X, SHIFT_Y);

    y = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, w, t), SCALE_X, SCALE_Y), M_PI_2f), SHIFT_X, SHIFT_Y);
    z = ShiftCoordinate(LinearRotate(LinearScale(BSRQS(_p0, _p1, _p2, -w, t), SCALE_X, SCALE_Y), M_PI_2f), SHIFT_X, SHIFT_Y);

    DrawLine(*grid, p, pp);
    DrawLine(*grid, q, pq);
    DrawLine(*grid, r, pr);
    DrawLine(*grid, s, ps);
    DrawLine(*grid, u, pu);
    DrawLine(*grid, v, pv);
    DrawLine(*grid, y, py);
    DrawLine(*grid, z, pz);

    pp = p;
    pq = q;
    pr = r;
    ps = s;
    pu = u;
    pv = v;
    py = y;
    pz = z;
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
    pr = ShiftCoordinate(LinearScale(_p0, SCALE_X, SCALE_Y), SHIFT_X, SHIFT_Y);
    ps = pr;
    pu = ShiftCoordinate(LinearRotate(LinearScale(_p0, SCALE_X, SCALE_Y), -ROTATE_ANGLE), SHIFT_X, SHIFT_Y);
    pv = pu;
    py = ShiftCoordinate(LinearRotate(LinearScale(_p0, SCALE_X, SCALE_Y), M_PI_2f), SHIFT_X, SHIFT_Y);
    pz = py;

    interpolatedPoints = {
        Point{0.0f , 40.0f},
        Point{5.0f , 50.0f},
        Point{6.0f , 20.0f},
        Point{8.0f , 00.0f},
        Point{12.0f, 80.0f},
        Point{16.0f, 60.0f},
        Point{20.0f, 90.0f},
        Point{22.0f, 70.0f},
        Point{23.0f, 30.0f},
        Point{30.0f, 70.0f},
    };


    std::vector<std::pair<Point, Point>> interpolates = CubicSplineInterpolation(interpolatedPoints.size(), interpolatedPoints, 0.0f, 0.0f);

    ipOutputs.resize(2*interpolatedPoints.size());
    for (size_t i = 0; i < interpolatedPoints.size(); i++) {
        ipOutputs[i*2 + 0] = interpolates[i].first;
        ipOutputs[i*2 + 1] = interpolates[i].second;
    }

    for (size_t i = 0; i < interpolatedPoints.size(); i++) {
        std::cerr << "Point (" << i << "): " << "x: " << interpolatedPoints[i].x << " y: " << interpolatedPoints[i].y << std::endl;
        std::cerr << "a (" << i << "): " << "x: " << ipOutputs[i*2 + 0].x << "y: " << ipOutputs[i*2 + 0].y << std::endl;
        std::cerr << "b (" << i << "): " << "x: " << ipOutputs[i*2 + 1].x << "y: " << ipOutputs[i*2 + 1].y << std::endl;
    }



    SetupAndLoop(Update);
    return 0;
}
