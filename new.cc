#include "setup.h"
using namespace std;

vector<Curve<Point>> todraw;

void Update(PixelGrid *grid, float t) {
    for (float t = 0.0f; t <= 1.0f; t += 0.005f) {
        Point f = TruncateCurve(todraw, t);
        DrawPoint(*grid, f);
    }
}



int main(void) {

    Point A = Point{7.339,2.928};
    Point B = Point{7.076,2.455};
    Point C = Point{6.466,2.046};
    Point D = Point{6.69,2.087};
    Point E = Point{7.02,2.088};

    vector<Point> interpoints = {A,B,C,D,E};

    auto interpolate = CubicSplineInterpolation(5, interpoints, 0,0);

    vector<Point> finalpoints;
    for (int i = 0; i < 4; i++) {
        finalpoints.push_back(interpoints[i]);
        finalpoints.push_back(interpolate[i].first);
        finalpoints.push_back(interpolate[i].second);
        finalpoints.push_back(interpoints[i+1]);
    }




    for (size_t i = 0; i < finalpoints.size(); i += 4) {
        if (i + 3 >= finalpoints.size()) break;
        Point p,q,r,s;
        p = finalpoints[i];
        q = finalpoints[i+1];
        r = finalpoints[i+2];
        s = finalpoints[i+3];
        todraw.push_back(Curve<Point>([p,q,r,s](float t) { return LinearScale(BernsteinCubicSpline(p, q, r, s, t), 80.0, 80.0);}));
    }


    SetupAndLoop(Update);
    return 0;

}
