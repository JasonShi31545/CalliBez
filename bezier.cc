#include "bezier.h"

// Vectors and Matrices

float Vector2::magnitude() {
    return sqrtf(x*x + y*y);
}

float Vector3::magnitude() {
    return sqrtf(x*x + y*y + z*z);
}
std::vector<std::vector<float>> Matrix2::getVectorForm() {
    return iv;
}
std::vector<std::vector<float>> Matrix3::getVectorForm() {
    return iv;
}

float dot(Vector2 a, Vector2 b) {
    return a.x * b.x + a.y * b.y;
}

float dot(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector2 dot(Matrix2 m, Vector2 v) {
    float a,b,c,d,x,y;
    a = m.get(0, 0);
    b = m.get(0, 1);
    c = m.get(1, 0);
    d = m.get(1, 1);
    x = v.x;
    y = v.y;
    return Vector2(a*x+b*y, c*x+d*y);
}

Vector3 dot(Matrix3 m, Vector3 v) {
    std::vector<float> res(3,0);
    for (int i = 0; i < 3; i++) { // row
        for (int j = 0; j < 3; j++) { // column
            float n = m.get(i,j);
            res[i] += n * v[i];
        }
    }
    return Vector3(res);
}

Matrix2 dot(Matrix2 m1, Matrix2 m2) {
    Matrix2 res{};
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res.set(i,j, dot(m1.getRow(i), m2.getColumn(j)));
        }
    }
    return res;
}

Matrix3 dot(Matrix3 m1, Matrix3 m2) {
    Matrix3 res{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            res.set(i,j, dot(m1.getRow(i), m2.getColumn(j)));
        }
    }
    return res;
}

float TimeTransform(float t) {
//    return (1 - expf(-0.05f * t));
//    return (1 - expf(-0.00005f * t));

    float val = t/(10000.0f);
    if (val <= 0.0f) {
        return 0.0f;
    } else if (val >= 1.0f) {
        return 1.0f;
    } else {
        return val;
    }
}


Point lerp(Point s, Point e, float t) {
    assert((0.0f <= t) && (t <= 1.0f));
    float x, y;
    x = (1-t) * s.x + t * e.x;
    y = (1-t) * s.y + t * e.y;
    return Point(x,y);
}

void DrawPoint(PixelGrid &g, Point p) {
    // SDL_RenderDrawPointF(r, p.x, W_HEIGHT - p.y)
    int rx = roundf(p.x);
    int ry = W_HEIGHT - round(p.y);
    if (rx < 0 || rx > W_WIDTH || ry < 0 || ry > W_HEIGHT) return;
    g[rx][ry] = (uint32_t)0xFFFFFFFF;
}

void DrawLine(PixelGrid &g, Point s, Point e) {
    // SDL_RenderDrawLineF(r, s.x, W_HEIGHT - s.y, e.x, W_HEIGHT - e.y);
    const float epsilon = 0.01f;
    for (float t = 0; t <= 1.0f; t += epsilon) {
        Point p = lerp(s,e,t);
        DrawPoint(g, p);
    }
}


void DrawCircle(PixelGrid &g, Point o, float radius) {
    float x,y;
    for (float t = 0.0f; t <= 2*M_PI; t += 0.01f) {
        x = radius*cosf(t) + o.x;
        y = radius*sinf(t) + o.y;
        // DrawPoint(r, Point(x,y));
        DrawPoint(g, Point(x,y));
    }
}

void DrawWidth(PixelGrid &g, Point o, Vector2 v, float w) {

    // v for velocity vector
    // w for width
    if (w <= 0.0f) {
        w = 0.5f;
    }
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
        // DrawPoint(r, Point(x1, y1));
        DrawPoint(g, Point(x1,y1));
        // DrawPoint(r, Point(x2, y2));
        DrawPoint(g, Point(x2,y2));
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
    for (unsigned int i = 0; i <= n; i++) {
        resx += BersteinCoefficient(n, i, t) * points[i].x;
        resy += BersteinCoefficient(n, i, t) * points[i].y;
    }
    return Point(resx, resy);
}

Point BezierCurveWeightedN(unsigned int n, std::vector<Point> points, std::vector<float> weights, float t) {
    assert(points.size() == n+1);
    assert(weights.size() == n+1);
    float resx = 0, resy = 0;
    for (unsigned int i = 0; i <= n; i++) {
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
    for (unsigned int i = 0; i <= n; i++) {
        float bc = BersteinCoefficient(n, i, t) * weights[i];
        den += bc;
        resx += bc * points[i].x;
        resy += bc * points[i].y;
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


float StandardizeWeight(float w0, float w1, float w2) {
    return w1 * sqrtf(1/ (w0 * w2));
}

float WeightFromShapeCoefficient(float k) {
    return k / (1 - k);
}

// Matrix



// 3D
