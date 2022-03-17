#include <cmath>
#include <utility>
#include <vector>
#include <iostream>
#include <initializer_list>
#include <limits>

using floating = double;

const floating PI = std::atan2(0, -1);
const floating INF = 1e9;
const floating EPS = 1e-6;

struct Vector;
class Line;

bool equal(floating x, floating y) {
    return std::abs(x - y) < EPS;
}

struct Point {
    floating x;
    floating y;

    explicit Point(floating x = 0, floating y = 0) : x(x), y(y) {}
    explicit Point(const Vector& a);

    [[nodiscard]] floating distance(const Point& a) const {
        return std::sqrt((x - a.x) * (x - a.x) + (y - a.y) * (y - a.y));
    }
    Point& rotate(Point centre, floating angle);
    Point& reflect(const Point& center);
    Point& reflect(const Line& axis);
};

template<typename V>
const bool has_coordinates = std::is_same_v < V, Vector
> or
std::is_same_v<V, Point>;

template<typename V, typename = std::enable_if_t<has_coordinates<V>>>
bool operator==(const V& a, const V& b) {
    return equal(a.x, b.x) && equal(a.y, b.y);
}

template<typename V, typename = std::enable_if<has_coordinates<V>>>
bool operator!=(const V& a, const V& b) {
    return !(a == b);
}

struct Vector {
    floating x;
    floating y;

    explicit Vector(floating x = 0, floating y = 0) : x(x), y(y) {}
    explicit Vector(const Point& a) : x(a.x), y(a.y) {}
    Vector(const Point& from, const Point& to) : x(to.x - from.x), y(to.y - from.y) {}

    Vector& operator*=(floating b);
    Vector operator-() const;

    Vector& normalize();
    Vector& rotate(floating angle);
    [[nodiscard]] floating length() const;
};

template<typename V, typename = std::enable_if_t<has_coordinates<V>>>
V& operator+=(V& x, const Vector& y) {
    x.x += y.x;
    x.y += y.y;
    return x;
}

template<typename V, typename = std::enable_if_t<has_coordinates<V>>>
V& operator-=(V& x, const Vector& y) {
    x.x -= y.x;
    x.y -= y.y;
    return x;
}

template<typename V, typename = std::enable_if_t<has_coordinates<V>>>
V operator+(const V& x, const Vector& y) {
    V res(x);
    res += y;
    return res;
}

template<typename V, typename = std::enable_if_t<has_coordinates<V>>>
V operator-(const V& x, const Vector& y) {
    V res(x);
    res -= y;
    return res;
}

Vector operator-(const Point& to, const Point& from) { // x - y = vector from y to x
    return {from, to};
}

Point::Point(const Vector& a) : x(a.x), y(a.y) {}

floating Vector::length() const {
    return std::sqrt(x * x + y * y);
}

Vector& Vector::operator*=(floating b) {
    x *= b;
    y *= b;
    return *this;
}

Vector operator*(floating a, const Vector& b) {
    Vector res = b;
    return (res *= a);
}

Vector operator*(const Vector& a, floating b) {
    Vector res = a;
    return (res *= b);
}

Point& Point::reflect(const Point& center) {
    Vector delta = Vector(*this, center) * 2;
    *this += delta;
    return *this;
}

Vector& Vector::normalize() {
    floating len = length();
    x /= len;
    y /= len;
    return *this;
}

Vector& Vector::rotate(floating angle) {
    floating nx = x * cos(angle) - y * sin(angle);
    floating ny = x * sin(angle) + y * cos(angle);
    x = nx;
    y = ny;
    return *this;
}

Vector Vector::operator-() const {
    return Vector(-x, -y);
}

floating operator*(const Vector& a, const Vector& b) {
    return a.x * b.x + a.y * b.y;
}

floating operator^(const Vector& a, const Vector& b) {
    return a.x * b.y - a.y * b.x;
}

Point& Point::rotate(Point centre, floating angle) {
    Vector v = Vector(centre, *this);
    v.rotate(angle);
    centre += v;
    return *this = centre;
}

// Биссектриса острого угла между a и b
Vector bisector(Vector a, Vector b) {
    if ((a * b) < 0)
        b = -b;
    a.normalize();
    b.normalize();
    return a + b;
}

Point midpoint(Point a, Point b) {
    return ((a) + (b - a) * 0.5);
}

// OK


class Line {
private:
    friend Point intersection(const Line& a, const Line& b);
    floating A = 1;
    floating B = 0;
    floating C = 0;
    friend bool operator==(const Line& a, const Line& b);
public:
    void normalize() {
        if (!equal(A, 0)) {
            B /= A;
            C /= A;
            A = 1.0;
        } else {
            C /= B;
            B = 1.0;
        }
    }
    Line() = default;
    Line(const Point& a, const Point& b) {
        A = a.y - b.y;
        B = b.x - a.x;
        C = a.x * b.y - b.x * a.y;
        normalize();
    }
    Line(floating k, floating b) : Line(Point(0, b), Point(1, k + b)) {}
    Line(const Point& a, floating k) : Line(a, Point(0, a.y - a.x * k)) {}

    //xd
    [[nodiscard]] Vector normal() const {
        return Vector(A, B);
    }

    [[nodiscard]] floating at(Point x) const {
        return x.x * A + x.y * B + C;
    }

    [[nodiscard]] floating distance(Point x) const {
        return std::abs(this->at(x)) / std::sqrt(A * A + B * B);
    }
};

Point& Point::reflect(const Line& axis) {
    Vector normal = axis.normal();
    if (axis.at(*this) > 0)
        normal *= -1;
    floating dist = 2 * axis.distance(*this);
    normal *= (dist / normal.length());
    *this += normal;
    return *this;
}

Point intersection(const Line& a, const Line& b) {
    return Point((b.C * a.B - a.C * b.B) / (a.A * b.B - b.A * a.B), (b.A * a.C - a.A * b.C) / (a.A * b.B - b.A * a.B));
}

bool operator==(const Line& a, const Line& b) {
    return equal(a.A, b.A) && equal(a.B, b.B) && equal(a.C, b.C);
}

bool operator!=(const Line& a, const Line& b) {
    return !(a == b);
}

class Shape {
public:
    [[nodiscard]] virtual floating perimeter() const = 0;
    [[nodiscard]] virtual floating area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    [[nodiscard]] virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual void rotate(Point center, floating angle) = 0;
    virtual void reflect(Point center) = 0;
    virtual void reflect(Line axis) = 0;
    virtual void scale(Point center, floating coefficient) = 0;
    virtual ~Shape() = default;
};

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    [[nodiscard]] Point centroid() const {
        Vector center;
        for (auto vertex: vertices) {
            center += Vector(vertex);
        }
        center.x /= static_cast<floating>(verticesCount());
        center.y /= static_cast<floating>(verticesCount());
        return Point(center);
    }
    Polygon() = default;
    explicit Polygon(std::vector<Point> vertices) : vertices(std::move(vertices)) {}
    Polygon(const std::initializer_list<Point>& vertices) : vertices(vertices) {}

    explicit Polygon(const Point& point) {
        vertices.push_back(point);
    }

    template<typename... Args>
    Polygon(const Point& point, Args... args): Polygon(args...) {
        vertices.insert(vertices.begin(), point);
    }

    [[nodiscard]] size_t verticesCount() const {
        return vertices.size();
    }

    [[nodiscard]] const std::vector<Point>& getVertices() const {
        return vertices;
    }

    [[nodiscard]] floating perimeter() const override {
        floating P = (vertices.back() - vertices[0]).length();
        for (size_t i = 1; i < vertices.size(); ++i) {
            P += (vertices[i] - vertices[i - 1]).length();
        }
        return P;
    }

    [[nodiscard]] floating area() const override {
        floating S = Vector(vertices[verticesCount() - 1]) ^ Vector(vertices[0]);
        for (size_t i = 0; i < vertices.size(); ++i) {
            S += (Vector(vertices[i]) ^ Vector(vertices[i + 1]));
        }
        return std::abs(S) / 2.0;
    }

    [[nodiscard]] bool isConvex() const {
        bool positive = (Vector(vertices[0], vertices[1]) ^ Vector(vertices[1], vertices[2])) > 0;
        for (size_t i = 0; i < vertices.size(); ++i) {
            if (((Vector(vertices[i], vertices[(i + 1) % verticesCount()])
                ^ Vector(vertices[(i + 1) % verticesCount()], vertices[(i + 2) % verticesCount()])) > 0) != positive)
                return false;
        }
        return true;
    }

    //Compares Polygons with respect to their position \ rotation in space
    bool operator==(const Shape& otherShape) const override {
        auto ptr = dynamic_cast<const Polygon*>(&otherShape);
        if (ptr == nullptr)
            return false;
        // Checks if other Shape is Polygon too

        const auto& otherPoly = dynamic_cast<const Polygon&>(otherShape);
        if (otherPoly.verticesCount() != verticesCount())
            return false;

        const size_t n = verticesCount();
        size_t ind = 0;
        for (; ind < n; ++ind) {
            if (otherPoly.vertices[ind] == vertices[0])
                break;
        }
        if (ind == n)
            return false;

        int traversalOrder;
        if (otherPoly.vertices[(ind + 1) % n] == vertices[1])
            traversalOrder = 1;
        else
            traversalOrder = -1;
        for (size_t i = 0; i < n; ++i) {
            if (otherPoly.vertices[ind] != vertices[i])
                return false;
            ind = (ind + n + traversalOrder) % n;
        }
        return true;
    }

    bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }

    //Compares Polygons ignoring their position \ rotation in space
    [[nodiscard]] bool isCongruentTo(const Shape& otherShape) const override {
        auto ptr = dynamic_cast<const Polygon*>(&otherShape);
        if (ptr == nullptr)
            return false;
        // Checks if other Shape is Polygon too
        const auto& otherPoly = dynamic_cast<const Polygon&>(otherShape);
        if (otherPoly.verticesCount() != verticesCount())
            return false;
        const size_t n = verticesCount();
        Polygon mySetPoints;
        Polygon otherSetPoints;
        for (size_t i = 0; i < n; ++i) {
            mySetPoints.vertices.emplace_back(std::abs(
                                                  (vertices[(i + 1) % n] - vertices[i]).normalize() ^
                                                      (vertices[(i - 1 + n) % n] - vertices[i]).normalize()
                                              ),
                                              (vertices[i] - vertices[(i + 1) % n]).length()
                                                  * (vertices[i] - vertices[(i - 1 + n) % n]).length());
            otherSetPoints.vertices.emplace_back(std::abs(
                                                     (otherPoly.vertices[(i + 1) % n] - otherPoly.vertices[i]).normalize() ^
                                                         (otherPoly.vertices[(i - 1 + n) % n] - otherPoly.vertices[i]).normalize()
                                                 ),
                                                 (otherPoly.vertices[i] - otherPoly.vertices[(i + 1) % n]).length()
                                                     * (otherPoly.vertices[i]
                                                         - otherPoly.vertices[(i - 1 + n) % n]).length());
        }
        return mySetPoints == otherSetPoints;
    }

    [[nodiscard]] bool containsPoint(Point point) const override {
        floating sum = 0;
        for (size_t i = 0; i < verticesCount(); ++i) {
            Vector a = Vector(point, vertices[i]);
            Vector b = Vector(point, vertices[(i + 1) % verticesCount()]);
            sum += atan2(a ^ b, a * b);
        }
        return !equal(sum, 0);
    }

    [[nodiscard]] bool isSimilarTo(const Shape& otherShape) const override {
        auto ptr = dynamic_cast<const Polygon*>(&otherShape);
        if (ptr == nullptr)
            return false;
        // Checks if other Shape is Polygon too

        const auto& otherPoly = dynamic_cast<const Polygon&>(otherShape);
        if (otherPoly.verticesCount() != verticesCount())
            return false;

        const size_t n = verticesCount();
        floating smallestSideWe = std::numeric_limits<floating>::max();
        floating smallestSideOther = smallestSideWe;
        for (size_t i = 0; i < n; ++i) {
            smallestSideWe = std::min((vertices[(i + 1) % n] - vertices[i]).length(), smallestSideWe);
            smallestSideOther =
                std::min((otherPoly.vertices[(i + 1) % n] - otherPoly.vertices[i]).length(), smallestSideOther);
        }
        Polygon nw = *this;
        nw.scale(Point(0, 0), smallestSideOther / smallestSideWe);
        return nw.isCongruentTo(otherPoly);
    }

    void rotate(Point center, floating angle) override {
        angle = PI * angle / 180.0;
        for (auto& d: vertices) {
            d = Point(center + (d - center).rotate(angle));
        }
    }

    void reflect(Point center) override {
        for (auto& d: vertices) {
            d.reflect(center);
        }
    }

    void reflect(Line axis) override {
        for (auto& d: vertices) {
            d.reflect(axis);
        }
    }

    void scale(Point center, floating coefficient) override {
        for (auto& d: vertices) {
            d = (center + Vector(center, d) * coefficient);
        }
    }
};

class Ellipse : public Shape {
protected:
    std::pair<Point, Point> focus;
    floating sumDistances = 0;
    [[nodiscard]] floating c() const {
        return (focus.first - focus.second).length() / 2;
    }

    [[nodiscard]] floating a() const {
        return sumDistances / 2;
    }

    [[nodiscard]] floating e() const {
        return c() / a();
    }

    [[nodiscard]] floating b() const {
        return std::sqrt(a() * a() - c() * c());
    }
public:
    Ellipse() = default;
    Ellipse(const Point& F1, const Point& F2, floating sumDistances) : focus(std::make_pair(F1, F2)),
                                                                       sumDistances(sumDistances) {}

    [[nodiscard]] std::pair<Point, Point> focuses() const {
        return focus;
    }

    [[nodiscard]] Point center() const {
        return midpoint(focus.first, focus.second);
    }

    [[nodiscard]] floating eccentricity() const {
        return e();
    }

    [[nodiscard]] std::pair<Line, Line> directrices() const {
        Vector resultDirection = Line(focus.first, focus.second).normal();
        Vector myDirection = Vector(focus.first, focus.second).normalize();
        floating k = a() * a() / c();
        return {Line(Point(center() + myDirection * k), Point(center() + myDirection * k + resultDirection)),
                Line(Point(center() - myDirection * k), Point(center() - myDirection * k + resultDirection))};
    }

    //approximation
    [[nodiscard]] floating perimeter() const override {
        return PI * (3 * (a() + b()) - std::sqrt((3 * a() + b()) * (a() + 3 * b())));
    }

    [[nodiscard]] floating area() const override {
        return PI * a() * b();
    }

    bool operator==(const Shape& otherShape) const override {
        if (dynamic_cast<const Ellipse*>(&otherShape) == nullptr)
            return false;

        const auto& otherEllipse = dynamic_cast<const Ellipse&>(otherShape);
        return (equal(sumDistances, otherEllipse.sumDistances))
            and ((otherEllipse.focus == focus)
                or (otherEllipse.focus.second == focus.first and otherEllipse.focus.first == focus.second));
    }

    bool operator!=(const Shape& another) const override {
        return !(*this == another);
    }

    [[nodiscard]] bool isCongruentTo(const Shape& otherShape) const override {
        if (dynamic_cast<const Ellipse*>(&otherShape) == nullptr)
            return false;

        const auto& otherEllipse = dynamic_cast<const Ellipse&>(otherShape);
        return (equal(sumDistances, otherEllipse.sumDistances) and equal(eccentricity(), otherEllipse.eccentricity()));
    }

    [[nodiscard]] bool isSimilarTo(const Shape& otherShape) const override {
        if (dynamic_cast<const Ellipse*>(&otherShape) == nullptr)
            return false;

        const auto& otherEllipse = dynamic_cast<const Ellipse&>(otherShape);
        return equal(eccentricity(), otherEllipse.eccentricity());
    }

    [[nodiscard]] bool containsPoint(Point point) const override {
        return sumDistances > point.distance(focus.first) + point.distance(focus.second) - EPS;
    }

    void rotate(Point center, floating angle) override {
        angle = PI * angle / 180;
        focus.first.rotate(center, angle);
        focus.second.rotate(center, angle);
    }

    void reflect(Point center) override {
        focus.first.reflect(center);
        focus.second.reflect(center);
    }

    void reflect(Line axis) override {
        focus.first.reflect(axis);
        focus.second.reflect(axis);
    }

    void scale(Point center, floating coefficient) override {
        focus.first = Point(center + (focus.first - center) * coefficient);
        focus.second = Point(center + (focus.second - center) * coefficient);
        sumDistances *= std::abs(coefficient);
    }
};

class Circle : public Ellipse {
public:
    Circle() = default;
    Circle(const Point& a, floating r) : Ellipse(a, a, 2 * r) {}

    [[nodiscard]] floating radius() const {
        return a();
    }
};

class Triangle : public Polygon {
public:
    Triangle() = default;
    Triangle(const Point& a, const Point& b, const Point& c) : Polygon({a, b, c}) {}
    [[nodiscard]] Point orthocenter() const {
        Line a = Line(vertices[0], vertices[1]);
        Line b = Line(vertices[0], vertices[2]);
        Line h1 = Line(vertices[2], Point(vertices[2] + a.normal()));
        Line h2 = Line(vertices[1], Point(vertices[1] + b.normal()));
        return intersection(h1, h2);
    }

    [[nodiscard]] Circle circumscribedCircle() const {
        Point centre = orthocenter() + (centroid() - orthocenter()) * 1.5;
        return {centre, (vertices[0] - centre).length()};
    }

    [[nodiscard]] Circle inscribedCircle() const {
        Vector v01 = Vector(vertices[0], vertices[1]);
        Vector v02 = Vector(vertices[0], vertices[2]);
        Vector v12 = Vector(vertices[1], vertices[2]);
        Vector v10 = -v01;
        Line bis0 = Line(vertices[0], Point(vertices[0] + bisector(v01, v02)));
        Line bis1 = Line(vertices[1], Point(vertices[1] + bisector(v12, v10)));
        Point centre = intersection(bis0, bis1);
        Line a = Line(vertices[0], vertices[1]);
        floating radius = a.distance(centre);
        return {centre, radius};
    }

    [[nodiscard]] Line EulerLine() const {
        return {orthocenter(), centroid()};
    }

    [[nodiscard]] Circle ninePointsCircle() const {
        Point a = midpoint(vertices[0], vertices[1]);
        Point b = midpoint(vertices[1], vertices[2]);
        Point c = midpoint(vertices[0], vertices[2]);
        return Triangle(a, b, c).circumscribedCircle();
    }
};

class Rectangle : public Polygon {
public:
    Rectangle() = default;
    Rectangle(const Point& A, const Point& C, floating k) {
        if (k < 1)
            k = 1 / k;
        Point B =
            Point(A + (
                ((C - A) * (1 / std::sqrt(1 + k * k)))
            ).rotate(std::atan(k)));
        Point D =
            Point(A + (((C - A) * (1 / std::sqrt(1 + k * k)))).rotate(std::atan(-k)));
        vertices = {A, B, C, D};
    }

    [[nodiscard]] Point center() const {
        return midpoint(vertices[0], vertices[2]);
    }

    [[nodiscard]] std::pair<Line, Line> diagonals() const {
        return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
    }
};

class Square : public Rectangle {
public:
    Square() = default;
    Square(const Point& a, const Point& c) : Rectangle(a, c, 1.0) {}

    [[nodiscard]] Circle circumscribedCircle() const {
        return {center(), (vertices[0] - vertices[2]).length() / 2};
    }

    Circle inscribedCircle() {
        return {center(), (vertices[0] - vertices[1]).length() / 2.0};
    }
};
























//padding