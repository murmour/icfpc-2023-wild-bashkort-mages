#pragma once
#include <cmath>

typedef double gbase;
const gbase geps = 1e-9;


struct Point
{
    gbase x, y;
    Point(gbase x = 0, gbase y = 0) : x(x), y(y) {}
    Point to(const Point &p) const { return Point(p.x - x, p.y - y); }
    gbase dot(const Point &p) const { return p.x * x + p.y * y; }
    gbase cross(const Point &p) const { return x * p.y - y * p.x; }
    gbase len2() const { return x*x + y*y; }
    double len() const { return sqrt(len2()); }
    Point normalized() const { double t = len(); return Point(x/t, y/t); }
    Point rot90() const { return {-y, x}; }
    friend Point operator + (const Point &a, const Point &b) { return Point(a.x + b.x, a.y + b.y); }
    friend Point operator - (const Point &a, const Point &b) { return Point(a.x - b.x, a.y - b.y); }
    friend Point operator * (gbase k, const Point &p) { return Point(p.x * k, p.y * k); }
    friend Point operator * (const Point &p, gbase k) { return Point(p.x * k, p.y * k); }
    friend bool operator < (const Point &u, const Point &v) { return false; }
};


struct Line
{
    gbase a, b, c;
    Line(gbase a, gbase b, gbase c) : a(a), b(b), c(c) {}
    static Line PP(const Point &u, const Point &v) { double a = u.y - v.y, b = v.x - u.x; return Line(a, b, -(a * u.x + b * u.y)); }
    static Line PN(const Point &p, const Point &n) { return Line(n.x, n.y, -(n.x * p.x + n.y * p.y)); }
    gbase at(const Point &p) const { return a * p.x + b * p.y + c; }
    gbase sigdist(const Point &p) const { return at(p) / sqrt(a*a+b*b); }
    Point isect(const Line &other) { double d = a*other.b - b*other.a; return Point((b*other.c - c*other.b) / d, (c*other.a - a*other.c) / d); }
};
