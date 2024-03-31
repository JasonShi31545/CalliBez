#ifndef BEZIER_H
#define BEZIER_H

#include "header.h"

#define W_WIDTH 1800
#define W_HEIGHT 1050

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

class HomogeneousPoint: public Point {
    public:
        float w;
        HomogeneousPoint (float x, float y, float z): Point() {
            this->x = x;
            this->y = y;
            w = z;
        }
        HomogeneousPoint(float z): Point() {
            w = z;
        }
        HomogeneousPoint(Point &p, float z): Point() {
            x = p.x;
            y = p.y;
            w = z;
        }
        HomogeneousPoint(Point &p): Point() {
            w = 1.0f;
        }
        HomogeneousPoint() : Point() {
            w = 0.0f;
        }
        Point projected() {
            return Point(x/w, y/w);
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
        Vector2(std::vector<float> sv) {
            assert(sv.size() == 2);
            x = sv[0]; y = sv[1];
        }
        Vector2 operator+(Vector2 v) {
            return Vector2(x+v.x, y+v.y);
        }
        Vector2 operator-(Vector2 v) {
            return Vector2(x-v.x, y-v.y);
        }
        float &operator[](size_t pos) {
            assert(pos >= 0 && pos < 2);
            if (pos == 0) return x;
            if (pos == 1) return y;
            throw "Vector2 Index Error";
        }
        const float &operator[](size_t pos) const {
            assert(pos >= 0 && pos < 2);
            if (pos == 0) return x;
            if (pos == 1) return y;
            throw "Vector2 Index Error";
        }
        float magnitude();
        Point toPoint() {
            return Point(x,y);
        }
};

class Vector3 {
    public:
        float x, y, z;
        Vector3(): x(0), y(0), z(0) {};
        Vector3(float a, float b, float c) : x(a), y(b), z(c) {};
        Vector3(std::vector<float> sv) {
            assert(sv.size() == 3);
            x = sv[0]; y = sv[1]; z = sv[2];
        }
        Vector3(HomogeneousPoint &p): x(p.x), y(p.y), z(p.w) {};
        Vector3 operator+(Vector3 v) {
            return Vector3(x+v.x, y+v.y, z+v.z);
        }
        Vector3 operator-(Vector3 v) {
            return Vector3(x-v.x, y-v.y, z+v.z);
        }
        float &operator[](size_t pos) {
            assert(pos >= 0 && pos < 3);
            if (pos == 0) return x;
            if (pos == 1) return y;
            if (pos == 2) return z;
            throw "Vector2 Index Error";
        }
        const float &operator[](size_t pos) const {
            assert(pos >= 0 && pos < 3);
            if (pos == 0) return x;
            if (pos == 1) return y;
            if (pos == 2) return z;
            throw "Vector2 Index Error";
        }
        float magnitude();
        HomogeneousPoint toPoint() {
            return HomogeneousPoint(x,y,z);
        }
};

class Matrix2 {
    // dimension 2x2
    // Matrix2[i][j] is ith row, jth column
    private:
        std::vector<std::vector<float>> iv;
    public:
        Matrix2() {
            iv.resize(2);
            for (int i = 0; i < 2; i++) {
                iv[i].resize(2);
                iv[i][0] = 0;
                iv[i][1] = 0;
            }
        }
        Matrix2(std::vector<std::vector<float>> v) {
            assert(v.size() == 2);
            for (int i = 0; i < 2; i++) {
                assert(v[i].size() == 2);
            }
            iv = v;
        }
        float get(size_t px, size_t py) {
            assert(px < 2 && py < 2);
            return iv[px][py];
        }
        void set(size_t px, size_t py, float n) {
            iv[px][py] = n;
        }
        void setRow(size_t r, Vector2 v) {
            iv[r][0] = v.x;
            iv[r][1] = v.y;
        }
        void setColumn(size_t c, Vector2 v) {
            iv[0][c] = v.x;
            iv[1][c] = v.y;
        }
        Vector2 getRow(size_t r) {
            return Vector2(iv[r][0], iv[r][1]);
        }
        Vector2 getColumn(size_t c) {
            return Vector2(iv[0][c], iv[1][c]);
        }
        std::vector<std::vector<float>> getVectorForm();
};

class Matrix3 {
    // dimension 3x3
    // Matrix3[i][j] is ith row, jth column
    private:
        std::vector<std::vector<float>> iv;
    public:
        Matrix3() {
            iv.resize(3);
            for (int i = 0; i < 3; i++) {
                iv[i].resize(3);
                iv[i][0] = 0;
                iv[i][1] = 0;
                iv[i][2] = 0;
            }
        }
        Matrix3(std::vector<std::vector<float>> v) {
            assert(v.size() == 3);
            for (int i = 0; i < 3; i++) {
                assert(v[i].size() == 3);
            }
            iv = v;
        }
        float get(size_t px, size_t py) {
            assert(px < 3 && py < 3);
            return iv[px][py];
        }
        void set(size_t px, size_t py, float n) {
            iv[px][py] = n;
        }
        void setRow(size_t r, Vector3 v) {
            iv[r][0] = v.x;
            iv[r][1] = v.y;
            iv[r][2] = v.z;
        }
        void setColumn(size_t c, Vector3 v) {
            iv[0][c] = v.x;
            iv[1][c] = v.y;
            iv[1][c] = v.z;
        }
        Vector3 getRow(size_t r) {
            return Vector3(iv[r][0], iv[r][1], iv[r][2]);
        }
        Vector3 getColumn(size_t c) {
            return Vector3(iv[0][c], iv[1][c], iv[2][c]);
        }
        std::vector<std::vector<float>> getVectorForm();
};


typedef std::vector<std::vector<uint32_t>> PixelGrid;

float TimeTransform(float t);
Point lerp(Point s, Point e, float t);
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

// Matrices and Vectors

Vector2 dot(Matrix2 m, Vector2 v);
Vector3 dot(Matrix3 m, Vector3 v);
Matrix2 dot(Matrix2 m1, Matrix2 m2);
Matrix3 dot(Matrix3 m1, Matrix3 m2);

// Rational Bezier Curves and Conic Sections

Point MidPoint(Point p, Point q);
float Distance(Point p, Point q);

const float FINITE_INF = FLT_MAX;

// Circle -- aka. Pain in the arse
Point UnitCircle(float t); // Parametric Unit Circle
float StandardizeWeight(float w0, float w1, float w2);

// float CircularWeightFromPoints(Point p0, Point p1, Point p2);
// float WeightFromShapeCoefficient(float k);
// float ShapeCoefficient(Point m, Point i, HomogeneousPoint c); // m: midpoint, i: intersection, c: P1 control point

std::tuple<HomogeneousPoint, Point, float> ConstructArc(Point p0, float angle);


Point BSRQS(Point p0, Point p1, Point p2, float w, float t); // Berstein Standardized Rational Quadratic Spline
Point BSRQV(Point p0, Point p1, Point p2, float w, float t); // BSRQ Velocity
Point BSRQA(Point p0, Point p1, Point p2, float w, float t); // BSRQ Acceleration
Point BSRQC(Point p0, Point p1, Point p2, float w, float t); // BSRQ Curvature

// Point CircleEllipsePointTangentForm(Point p0, Point p2, Vector2 t0, Vector2 t2); // Use vectors and points to calculate intersections and thus points for projection

// Affine Transformations

HomogeneousPoint AffineTransformation(HomogeneousPoint input, Matrix3 transm);
HomogeneousPoint AffineTranslate(HomogeneousPoint input, float x, float y);
HomogeneousPoint AffineRotate(HomogeneousPoint input, float angle); // rad
HomogeneousPoint AffineScale(HomogeneousPoint input, float x, float y);
HomogeneousPoint AffineReflectionO(HomogeneousPoint input);
HomogeneousPoint AffineReflectionX(HomogeneousPoint input);
HomogeneousPoint AffineReflectionY(HomogeneousPoint input);
HomogeneousPoint AffineShearX(HomogeneousPoint input, float angle);
HomogeneousPoint AffineShearY(HomogeneousPoint input, float angle);


// Linear Transformations

Point LinearTransformation(Point input, Matrix2 tranm);
Point LinearReflectionX(Point input);
Point LinearReflectionY(Point input);
Point LinearReflectionLineAngle(Point input, float angle /* rad */);
Point LinearScale(Point input, float x, float y);
Point LinearShearX(Point input, float k);
Point LinearShearY(Point input, float k);
Point LinearRotate(Point input, float angle /* rad */);

// Basic Translate Coz Im lazy

Point ShiftCoordinate(Point input, float x, float y);

// Drawing functions
void DrawPoint(PixelGrid &g, Point p);
void DrawLine(PixelGrid &g, Point s, Point e);
void DrawCircle(PixelGrid &g, Point o, float radius);
void DrawWidth(PixelGrid &g, Point o, Vector2 v, float w);



#endif // BEZIER_H
