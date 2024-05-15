#include "setup.h"

#ifdef MAIN

using namespace std;

vector<vector<Point>> aps, fps, rps, eps;
vector<vector<Point>> aws, fws, rws, ews;

vector<vector<Point>> Fpoints() {

    Point A = Point{7.339,2.928};
    Point A1 = Point{7.292,2.913};
    Point A2 = Point{7.14,2.6};
    Point B = Point{7.076,2.455};
    Point B1 = Opposing(B,A2);
    Point B2 = Point{6.925,2.037};
    Point C = Point{6.466,2.046};
    Point C1 = Opposing(C,B2);
    Point C2 = Point{6.439,2.09};
    Point D = Point{6.69,2.087};
    Point D1 = Opposing(D, C2);
    Point D2 = Point{6.828,2.091};
    Point E = Point{7.02,2.088};
    Point F = Point{7.087,2.742};
    Point F1 = Point{7.014,2.579};
    Point F2 = Point{6.873,2.369};
    Point G = Point{6.775,2.355};
    Point G1 = Opposing(G,F2);
    Point G2 = Point{6.672,2.499};
    Point H = Point{6.685,2.556};
    Point H1 = Opposing(H,G2);
    Point H2 = Point{6.76,2.86};
    Point I = Point{7.083,2.82};
    Point I1 = Opposing(I, H2);
    Point I2 = Point{7.486,2.812};
    Point J = Point{7.528,2.855};
    //Point J1 = Opposing(J, I2);
    Point K = Point{7.572,2.912};
    Point K1 = Point{7.505,2.865};
    Point K2 = Point{7.531,2.752};
    Point L = Point{7.459,2.724};
    //Point L1 = Opposing(L, K2);
    Point M = Point{6.984,2.465};
    Point M1 = Point{7.032,2.501};
    Point M2 = Point{7.12,2.44};
    Point N = Point{7.16,2.47};
    //Point N1 = Opposing(N,M2);
    Point O = Point{7.174,2.485};
    Point O1 = Point{7.146,2.4547};
    Point O2 = Point{7.137,2.4304};
    Point P = Point{7.127,2.408};
    //Point P1 = Opposing(P, O2);



     return vector<vector<Point>>({
         vector<Point>({A,A1,A2,B,B1,B2,C,C1,C2,D,D1,D2,E}),
         vector<Point>({F,F1,F2,G,G1,G2,H,H1,H2,I,I1,I2,J}),
         vector<Point>({K,K1,K2,L}),
         vector<Point>({M,M1,M2,N}),
         vector<Point>({O,O1,O2,P})
     });
}

vector<vector<Point>> Epoints() {


    Point A = Point{8.06,2.459};
    Point A1 = Point{7.998,2.664};
    Point A2 = Point{7.81,2.68};
    Point B = Point{7.774,2.68};
    Point B1 = Opposing(B,A2);
    Point B2 = Point{7.5954,2.655};
    Point C = Point{7.59,2.54};
    Point C1 = Opposing(C,B2);
    Point C2 = Point{7.6855,2.3764};
    Point D = Point{7.79,2.41};
    Point D1 = Opposing(D,C2);
    Point D2 = Point{8.146,2.558};
    Point E = Point{8.232,2.704};
    Point E1 = Opposing(E,D2);
    Point E2 = Point{8.2044,2.859};
    Point F = Point{8.144,2.843};
    Point F1 = Opposing(F,E2);
    Point F2 = Point{7.9355,2.754};
    Point G = Point{7.882,2.57};
    Point G1 = Opposing(G,F2);
    Point G2 = Point{7.833,2.2184};
    Point H = Point{7.879,2.206};
    Point H1 = Opposing(H,G2);
    Point H2 = Point{7.957,2.217};
    Point I = Point{7.96,2.237};
    Point I1 = Opposing(I,H2);
    Point I2 = Point{7.954,2.291};
    Point J = Point{7.903,2.276};
    Point J1 = Opposing(J,I2);
    Point J2 = Point{7.633,2.142};
    Point K = Point{7.636,1.927};
    Point K1 = Opposing(K,J2);
    Point K2 = Point{7.74,1.7};
    Point L = Point{7.785,1.701};
    Point L1 = Opposing(L,K2);
    Point L2 = Point{7.962,1.759};
    Point M = Point{8.028,1.91};
    Point M1 = Opposing(M,L2);
    Point M2 = Point{8.0077,2.0995};
    Point N = Point{8,2.1};
    Point N1 = Opposing(N,M2);
    Point N2 = Point{7.912,2.0976};
    Point O = Point{7.867,1.9};
    //Point O1 = Opposing(O,N2);
    //Point O2 = Point{-1,-1};
    return vector<vector<Point>>({vector<Point>({A,A1,A2,B,B1,B2,C,C1,C2,D,D1,D2,E,E1,E2,F,F1,F2,G,G1,G2,H,H1,H2,I,I1,I2,J,J1,J2,K,K1,K2,L,L1,L2,M,M1,M2,N,N1,N2,O})});

}

vector<vector<Point>> Apoints() {


    Point A = Point{ 4.444,1.92};
    Point A1 = Point{ 4.346,1.93};
    Point A2 = Point{ 4.1662,1.6017};
    Point B = Point{ 4.094,1.468};
    Point B1 = Opposing(B,A2);
    Point B2 = Point{ 3.872,1.062};
    Point C = Point{ 3.718,1.046};
    Point C1 = Point{ 3.589,1.048};
    Point C2 = Point{ 3.629,1.306};
    Point D = Point{ 3.668,1.396};
    Point D1 = Opposing(D, C2);
    Point D2 = Point{ 3.834,1.67};
    Point E = Point{ 3.9465,1.642};
    Point E1 = Opposing(E, D2);
    Point E2 = Point{ 3.932,1.256};
    Point F = Point{ 3.746,1.219};
    Point G = Point{ 4.048,1.049};
    Point H = Point{ 4.212,1.43};
    Point H1 = Point{ 4.214,1.4816};
    Point H2 = Point{ 4.1835,1.54};
    Point I = Point{ 4.152,1.564};
    Point I1 = Opposing(I, H2);
    Point I2 = Point{ 4.051,1.596};
    Point J = Point{ 4.031,1.531};
    Point J1 = Opposing(J, I2);
    Point J2 = Point{ 4.036,1.425};
    Point K = Point{ 4.053,1.399};
    Point K1 = Opposing(K, J2);
    Point K2 = Point{ 4.1144,1.313};
    Point L = Point{ 4.173,1.303};
    Point L1 = Opposing(L, K2);
    Point L2 = Point{ 4.346,1.3244};
    Point M = Point{ 4.414,1.462};

      return vector<vector<Point>>({
          vector<Point>({A,A1,A2,B,B1,B2,C,C1,C2,D,D1,D2,E,E1,E2,F}),
          vector<Point>({A,A,G,G}),
          vector<Point>({H,H1,H2,I,I1,I2,J,J1,J2,K,K1,K2,L,L1,L2,M})
      });
}

vector<vector<Point>> Rpoints() {

    Point A = Point{1.48,5.756};
    Point A1 = Point{1.4316,5.756};
    Point A2 = Point{1.3556,5.585};
    Point B = Point{1.317,5.515};
    Point B1 = Opposing(B,A2);
    Point B2 = Point{1.1825,5.245};
    Point C = Point{1.108,5.296};
    //Point C1 = Opposing(C, B2);
    Point D = Point{1.312,5.658};
    Point D1 = Point{1.28,5.6};
    Point D2 = Point{1.227,5.498};
    Point E = Point{1.19,5.505};
    Point E1 = Opposing(E,D2);
    Point E2 = Point{1.1725,5.565};
    Point F = Point{1.174,5.574};
    Point F1 = Opposing(F,E2);
    Point F2 = Point{1.229,5.74};
    Point G = Point{1.36,5.735};
    Point G1 = Opposing(G,F2);
    Point G2 = Point{1.4727,5.6124};
    Point H = Point{1.4712,5.6057};
    Point H1 = Opposing(H,G2);
    Point H2 = Point{1.4276,5.4806};
    Point I = Point{1.367,5.497};
    Point I1 = Opposing(I, H2);
    Point I2 = Point{1.327,5.5269};
    Point J = Point{1.344,5.513};
    Point J1 = Opposing(J, I2);
    Point J2 = Point{1.3656,5.4785};
    Point K = Point{1.347,5.403};
    Point K1 = Opposing(K, J2);
    Point K2 = Point{1.3424,5.3195};
    Point L = Point{1.354,5.317};
    Point L12 = Opposing(L, K2);
    Point L2 = Point{1.412,5.331};
    Point M = Point{1.444,5.424};
    //Point M1 = Opposing(M,L2);
    //Point M2 = Point{2.5,5.6};

    return vector<vector<Point>>({vector<Point>({A,A1,A2,B,B1,B2,C}),vector<Point>({D,D1,D2,E,E1,E2,F,F1,F2,G,G1,G2,H,H1,H2,I,I1,I2,J,J1,J2,K,K1,K2,L,L12,L2,M})});
}

vector<vector<Point>> Fweights() {
    return
        {{
        Point{0.115,1.417},
        Point{0.615,1.439},
            Point{0.404,0.024}},{
        Point{0.007,2.34},
        Point{0.287,0.027},
            Point{0.029,0.042}},{
        Point{0.018,1.22},
        Point{-0.013,1.643},
            Point{0.192,0.044}},{
        Point{0,0.3},
        Point{0.5,0},
            Point{1,0.3}},{
        Point{0.017,0.992},
        Point{0.01,1.86},
            Point{0.211,0.019}}
};

}

vector<vector<Point>> Eweights() {
    return {{Point{0.226,0.323},
        Point{0.645,1.879},
            Point{0.944,0.78}},{
        Point{0,0.1},
        Point{0.449,0.074},
            Point{0.996,0.119}},{
        Point{0.67,-0.008},
        Point{0.49,2.66},
            Point{0.3,0}},{
        Point{0.55,-0.1},
        Point{0.21,2.708},
            Point{0.4,0}},{
        Point{0.084,0.86},
        Point{1.138,1.027},
        Point{0.998,1.488}}};
}
vector<vector<Point>> Aweights() {
    return {{Point{0.55,0},
        Point{0.5,2.67},
                Point{0.448,0}},{

        Point{0.647,0.02},
        Point{0.737,0.02},
            Point{0.7,2.35}},{

        Point{-0.018,1.92},
        Point{0.476,0.04},
        Point{1.015,1.91}},{

        Point{0.6,-0.08},
        Point{0.335,2.704},
            Point{0.208,-0.037}}};
}
vector<vector<Point>> Rweights() {
    return {{Point{0.43,0.018},
        Point{0.503,2.62},
            Point{0.57,0.017}}
            ,{
        Point{0.064,2.27},
        Point{0.2,0.1},
            Point{0.06,0.04}}
            ,{
        Point{0.6,0},
        Point{0.518,2.636},
            Point{0.4,0}}
            ,{
        Point{0.056,2.868},
        Point{1.06,-1.245},
        Point{0.067,0.36}}};
}


void Draw(PixelGrid *grid, float t, vector<vector<Point>> &ps, vector<vector<Point>> &ws, float wf, float scalex, float scaley, float shiftx, float shifty, char ch) {

    for (float tmp = 0.0f; tmp < 1.0f; tmp += 0.005f) {
        Point o,e,p1,p2;
        float width;


        // this draws every single stroke
        for (size_t stroke = 0; stroke < ps.size(); stroke++) {
            vector<Point> st = ps[stroke];
            assert(ws[stroke].size() == 3);

            const Point i = Point(0,0), f = Point(1,0);
            Point j,k,l;

            // figure out the stroke widths given the time
            switch (ch) {
                case 'F': {
                    assert(ps.size() == ws.size());
                    j = ws[stroke][0], k = ws[stroke][1], l = ws[stroke][2];
                    width = wf * BernsteinQuarticSpline(i, j, k, l, f, tmp).y;
                    break;
                }
                case 'E': {
                    assert(ps.size() == 1);
                    vector<Curve<Point>> curves;
                    for (size_t estrokes = 0; estrokes < 5; estrokes++) {
                        j = ws[estrokes][0], k = ws[estrokes][1], l = ws[estrokes][2];
                        curves.push_back(Curve<Point>{[i,j,k,l,f](float time) { return BernsteinQuarticSpline(i, j, k, l, f, time);}});
                    }
                    width = wf * ((Point)TruncateCurve(curves, tmp)).y;
                    break;
                }
                case 'A': {
                    assert(ps.size() == 3);
                    if (stroke == 0) {
                        vector<Curve<Point>> curves;
                        for (size_t astrokes = 0; astrokes <= 2; astrokes ++) {
                            j = ws[astrokes][0], k = ws[astrokes][1], l = ws[astrokes][2];
                            curves.push_back(Curve<Point>{[i,j,k,l,f](float time) { return BernsteinQuarticSpline(i, j, k, l, f, time);}});
                        }
                        width = wf * ((Point)TruncateCurve(curves, tmp)).y;
                    } else {
                        j = ws[stroke+1][0], k = ws[stroke+1][1], l = ws[stroke+1][2];
                        width = wf * BernsteinQuarticSpline(i, j, k, l, f, tmp).y;
                    }
                    break;
                }
                case 'R': {
                    assert(ps.size() == 2);
                    if (stroke != 0) {
                        vector<Curve<Point>> curves;
                        for (size_t rstrokes = 1; rstrokes < 4; rstrokes++) {
                            j = ws[rstrokes][0], k = ws[rstrokes][1], l = ws[rstrokes][2];
                            curves.push_back(Curve<Point>{[i,j,k,l,f](float time) { return BernsteinQuarticSpline(i, j, k, l, f, time);}});
                        }
                        width = wf * ((Point)TruncateCurve(curves, tmp)).y;

                    } else { // stroke == 0
                        j = ws[stroke][0], k = ws[stroke][1], l = ws[stroke][2];
                        width = wf * BernsteinQuarticSpline(i, j, k, l, f, tmp).y;
                    }
                    break;
                }
                default:
                    return;
            }

            // figure out the strokes, given the widths.
            assert((st.size() - 1)%3==0);

            Point curve;
            Vector2 velocity;

            vector<Curve<Point>> pointCurves;
            vector<Curve<Vector2>> vectorCurves;

            for (size_t i = 0; i < st.size(); i+=3) {
                if (i+3 >= st.size()) break;
                o = st[i];
                p1 = st[i+1];
                p2 = st[i+2];
                e = st[i+3];
                Curve<Point> c{[o,p1,p2,e, scalex, scaley, shiftx, shifty](float time) {
                    Point cv = BernsteinCubicSpline(o, p1, p2, e, time);
                    cv = LinearScale(cv, scalex, scaley);
                    cv = ShiftCoordinate(cv, shiftx, shifty);
                    return cv;
                }};
                Curve<Vector2> v{[o,p1,p2,e](float time) {
                    return BernsteinCubicVelocity(o, p1, p2, e, time);
                }};
                pointCurves.push_back(c);
                vectorCurves.push_back(v);
            }


            curve = TruncateCurve(pointCurves, tmp);
            velocity = TruncateCurve(vectorCurves, tmp);

            cerr << "X: " << curve.x << " Y: " << curve.y << endl;
            cerr << "VX: " << velocity.x << " VY: " << velocity.y << endl;

            DrawWidth(*grid, curve, velocity, width);

            // for (size_t i = 0; i < st.size(); i+=3) {
            //     if (i+3 >= st.size()) break;
            //     o = st[i];
            //     p1 = st[i+1];
            //     p2 = st[i+2];
            //     e = st[i+3];
            //     Point curve = BernsteinCubicSpline(o, p1,p2 ,e ,tmp);
            //     Vector2 vel = BernsteinCubicVelocity(o, p1, p2, e, tmp);
            //     curve = LinearScale(curve, scalex, scaley);
            //     curve = ShiftCoordinate(curve, shiftx, shifty);
            //     // DrawPoint(*grid, curve);
            //     DrawWidth(*grid, curve, vel, width);
            // }
        }

    }
}

bool drawn;
void Update(PixelGrid *grid, float t) {

    if (drawn) return;
    fps = Fpoints();
    eps = Epoints();
    aps = Apoints();
    rps = Rpoints();
    fws = Fweights();
    ews = Eweights();
    aws = Aweights();
    rws = Rweights();

    Draw(grid, t, fps, fws, 2.0f, 100.0f, 100.0f, -30, 30, 'F');
    Draw(grid, t, eps, ews, 2.0f, 100.0f, 100.0f, -30, 30, 'E');
    Draw(grid, t, aps, aws, 2.0f, 100.0f, 100.0f, -30, 30, 'A');
    Draw(grid, t, rps, rws, 2.0f, 100.0f, 100.0f, -30, 30, 'R');




    drawn = true;

}
int main(int argc, const char *argv[]) {
    

    drawn = false;

    SetupAndLoop(Update);
    return 0;
}

#endif

#ifdef TESTING

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

    Point A = Point{1.48,5.756};
    Point B = Point{1.317,5.515};
    Point C = Point{1.108,5.296};

    vector<Point> interpoints = {A,B,C};

    auto interpolate = CubicSplineInterpolation(interpoints.size(), interpoints, 0,0);

    vector<Point> finalpoints;
    for (size_t i = 0; i < interpoints.size()-1; i++) {
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
        todraw.push_back(Curve<Point>([p,q,r,s](float t) { return LinearScale(BernsteinCubicSpline(p, q, r, s, t), 50.0, 50.0);}));
    }


    SetupAndLoop(Update);
    return 0;

}

#endif
