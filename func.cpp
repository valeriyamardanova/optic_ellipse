#include <iostream>
#include <cmath>
#include <iomanip>

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {};
};

double epsilon = 10.0 * std::numeric_limits<double>::epsilon();

class Ellipce {
private:
    double x0, y0, a, b, c;
    Point f1, f2;
public:
    Ellipce(double x0, double y0, double a, double b) : x0(x0), y0(y0), a(a), b(b) {
        c = std::sqrt(std::abs(a * a - b * b));
        if (a > b) {
            f1 = Point(-c, 0);
            f2 = Point(c, 0);
        }
        else {
            f1 = Point(0, -c);
            f2 = Point(0, c);
        }
    }
    bool onEllipce(Point& m) const {
        if (fabs(m.x * m.x / (a * a) + m.y * m.y / (b * b) - 1) < epsilon)return true;
        return false;
    }
    Point fromNormalToShtrih(Point& a) {
        return Point(a.x - x0, a.y - y0);
    }
    Point fromShtrihToNormal(Point& a) {
        return Point(a.x + x0, a.y + y0);
    }
    Point pointPer(Point& m, double& angle) {
        double dx = std::cos(angle);
        double dy = std::sin(angle);
        double A = dx * dx / (a * a) + dy * dy / (b * b);
        double B = 2 * ((m.x * dx / (a * a)) + m.y * dy / (b * b));
        double C = m.x * m.x / (a * a) + m.y * m.y / (b * b) - 1;
        double t1, t2, t = 0;
        double max = std::max(fabs(B * B), fabs(4 * A * C));
        if (max < epsilon) {
            t = 0;
            return m;
        }
        double norm_D = (B * B - 4 * A * C) / max;
        if (norm_D > epsilon) {
            double D = norm_D * max;
            t1 = ((-B + sqrt(D)) / (2 * A));
            t2 = ((-B - sqrt(D)) / (2 * A));
            if (fabs(t1 - t2) > epsilon) {
                if (t1 > epsilon && t2 < epsilon) t = t1;
                else if (t1 < epsilon && t2 > epsilon) t = t2;
                else if (t1 > epsilon && t2 > epsilon) t = std::min(t1, t2);
                else if ((t1 < epsilon && t2 > -epsilon) || (t1 > -epsilon && t2 < epsilon)) t = 0;
                else t = 1;
            }
            else if (t1 < -epsilon) t = 1;
            else t = std::min(t1, t2);
        }
        if (norm_D < -(epsilon)) {
            t = 1;
        }
        if (fabs(norm_D) < epsilon) {
            t = -B / (2 * A);
            if (fabs(t) < epsilon) t = 0;
            else if (t < -epsilon) t = 1;
        }
        return Point(m.x + dx * t, m.y + dy * t);
    }
    Point reflect(Point& m, double& angle) {
        Point p = pointPer(m, angle);
        if (fabs(m.x - p.x) < epsilon && fabs(m.y - p.y) < epsilon) return m;
        double xn = -b * b / p.x;
        double yn = a * a / p.y;
        double xk = p.x / (a * a);
        double yk = p.y / (b * b);
        double t = (p.x - m.x) / (xk - xn);
        Point l;
        l.x = m.x + 2 * xk * t;
        l.y = m.y + 2 * yk * t;
        return l;
    }
    double distFromPointToLine(Point& f, Point& p, Point& l) {
        double A = p.y - l.y;
        double B = l.x - p.x;
        double C = l.y * p.x - l.x * p.y;
        double t = ((l.y - p.y) * (f.y - p.y) + (l.x - p.x) * (f.x - p.x));
        if (t > 0) return (std::abs(A * f.x + B * f.y + C) / std::sqrt(A * A + B * B));
        return std::sqrt((p.x - f.x) * (p.x - f.x) + (p.y - f.y) * (p.y - f.y));
    }
    void resh(Point M, double& angle, double& d1, double& d2, double& x, double& y) {
        M = fromNormalToShtrih(M);
        Point p = pointPer(M, angle);
        Point l = reflect(M, angle);
        if (onEllipce(p) || onEllipce(M)) {
            d1 = distFromPointToLine(f1, p, l);
            d2 = distFromPointToLine(f2, p, l);
        }
        else {
            d1 = distFromPointToLine(f1, M, l);
            d2 = distFromPointToLine(f2, M, l);
        }
        p = fromShtrihToNormal(p);
        x = p.x;
        y = p.y;
    }
};

int main()
{
    double x0, y0, a, b, xi, eta;
    int m;
    std::cin >> x0 >> y0 >> a >> b;
    std::cin >> xi >> eta;
    std::cin >> m;
    Point M(xi, eta);
    Ellipce ellipce(x0, y0, a, b);
    for (int j = 0; j < m; ++j) {
        double d1 = 0, d2 = 0, x, y;
        double angle = 2 * j * acos(-1) / m;
        ellipce.resh(M, angle, d1, d2, x, y);
        std::cout << x << " " << y << " " << d1 << " " << d2 << std::endl;
    }
}
