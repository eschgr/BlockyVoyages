#include <math.h>

namespace BlockyVoyages {
    // methods
    // constructors
    template<class T> inline Vector2<T>::Vector2(void) {}
    
    template<class T>
    inline Vector2<T>::Vector2(const Vector2<T>& other)
        : x(other.x), y(other.y)
    {}

    template<class T>
    inline Vector2<T>::Vector2(const Vector3<T>& other)
        : x(other.x), y(other.y)
    {}
    
    template<class T>
    inline Vector2<T>::Vector2(const Vector4<T>& other)
        : x(other.x), y(other.y)
    {}

    // set from floats
    template<class T>
    inline Vector2<T>::Vector2(T _x, T _y)
        : x(_x), y(_y)
    {}

    template<class T>
    inline Vector2<T>::Vector2(const T _v[2])
        : x(_v[0]), y(_v[1])
    {}

    template<class T>
    inline Vector2<T>::~Vector2(void) {}

    // setting methods
    template<class T>
    inline const Vector2<T>& Vector2<T>::operator=(const Vector2<T>& other) {
        if (&other != this) {
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    template<class T>
    inline void Vector2<T>::Set(const Vector2<T>& other) {
        x = other.x;
        y = other.y;
    }

    template<class T>
    inline void Vector2<T>::Set(const Vector3<T>& other) {
        x = other.x;
        y = other.y;
    }

    template<class T>
    inline void Vector2<T>::Set(const Vector4<T>& other) {
        x = other.x;
        y = other.y;
    }

    // set from floats
    template<class T>
    inline void Vector2<T>::Set(T _x, T _y) {
        x = _x;
        y = _y;
    }

    template<class T>
    inline void Vector2<T>::Set(const T _v[2]) {
        x = _v[0];
        y = _v[1];
    }

    // clean up
    template<class T>
    inline void Vector2<T>::Zero(void) {
        x = 0.0f;
        y = 0.0f;
    }

    template<class T>
    inline void Vector2<T>::Clean(void) {
        if (fabs(x) <EPSILON) {
            x = EPSILON;
        }

        if (fabs(y) <EPSILON) {
            y = EPSILON;
        }
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Add(const Vector2<T>& v1, const Vector2<T>& v2) {
        x = v1.x + v2.x;
        y = v1.y + v2.y;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Sub(const Vector2<T>& v1, const Vector2<T>& v2) {
        x = v1.x - v2.x;
        y = v1.y - v2.y;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Mul(const Vector2<T>& v1, T scalar) {
        x = v1.x * scalar;
        y = v1.y * scalar;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Mul(const Vector2<T>& v1, const Matrix22<T>& mat) {
        x = v1.x * mat.m00 + v1.y * mat.m01;
        y = v1.x * mat.m10 + v1.y * mat.m11;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Mul(const Matrix22<T>& mat, const Vector2<T>& v1) {
        x = v1.x * mat.m00 + v1.y * mat.m10;
        y = v1.x * mat.m01 + v1.y * mat.m11;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Div(const Vector2<T>& v1, T scalar) {
        x = v1.x / scalar;
        y = v1.y / scalar;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Build(const Vector2<T>& init, const Vector2<T>& term) {
        Sub(term, init);
        return *this;
    }

    template<class T>
    inline T Vector2<T>::Dot(const Vector2<T>& v2) const {                 // computes the dot product between this and v2
        return x * v2.x + y * v2.y;
    }

    template<class T>
    inline const Matrix22<T> Vector2<T>::Tensor(const Vector2<T>& v2) const {
        return Matrix22<T>(x * v2.x, x * v2.y,
                                  y * v2.x, y * v2.y);
    }

    template<class T>
    inline const Vector2<T> Vector2<T>::operator + (const Vector2<T> &rhs) const {        // Adds two vectors. use Add() or += instead if possible
        return Vector2<T>(x + rhs.x, y + rhs.y);
    }

    template<class T>
    inline const Vector2<T> Vector2<T>::operator - (const Vector2<T> &rhs) const {        // subtracts two vectors. use Sub() or -= instead if possible
        return Vector2<T>(x - rhs.x, y - rhs.y);
    }

    template<class T>
    inline T Vector2<T>::operator * (const Vector2<T> &rhs) const {        // computes the dot product between this and rhs
        return x * rhs.x + y * rhs.y;
    }

    template<class T>
    inline const Vector2<T> Vector2<T>::operator * (T scalar) const {             // multiplies this by a scalar. use Mul() or *= instead if possible
        return Vector2<T>(x * scalar, y * scalar);
    }

    template<class T>
    inline const Vector2<T> Vector2<T>::operator * (const Matrix22<T>& rhs) const {
        return Vector2<T>(x * rhs.m00 + y * rhs.m10,
                                 x * rhs.m01 + y * rhs.m11);
    }

    template<class T>
    inline const Vector2<T> Vector2<T>::operator / (T scalar) const {             // divides this by a scalar. use Div() or /= instead if possible
        return Vector2<T>(x / scalar, y / scalar);
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::operator += (const Vector2<T> &rhs) {     // does an in-place addition with rhs
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::operator -= (const Vector2<T> &rhs) {     // does an in-place subtraction with rhs
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::operator *= (T scalar) {          // does an in-place scalar multiplication
        x *= scalar;
        y *= scalar;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::operator *= (const Matrix22<T>& mat) {          // does an in-place matrix multiplication
        T tX = x * mat.m00 + y * mat.m10;
        y = x * mat.m01 + y * mat.m11;
        x = tX;
        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::operator /= (T scalar) {          // does an in-place scalar division
        x /= scalar;
        y /= scalar;
        return *this;
    }

    template<class T>
    inline T Vector2<T>::LengthSquared(void) const {
        return x * x + y * y;
    }

    template<class T>
    inline T Vector2<T>::Length(void) const {
        return sqrt(LengthSquared());
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::NegationOf(const Vector2<T>& other) {
        x = -other.x;
        y = -other.y;

        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Negate(void) {
        x *= -1;
        y *= -1;

        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::UnitOf(const Vector2<T>& other) {
        Set(other);
        Normalize();

        return *this;
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::Normalize(void) {
        T len = LengthSquared();
        if (len > TINY) {
            len = sqrt(len);
            x /= len;
            y /= len;
        }
        else {
            Zero();
        }

        return *this;
    }

    template<class T>
    inline T Vector2<T>::CosAngleBetween(const Vector2<T>& v2) const { // computes the cosine of the angle between this vector and vect
        T len1 = LengthSquared();
        T len2 = v2.LengthSquared();
        if (len1 > TINY && len2 > TINY) {
            return Dot(v2) / (sqrt(len1 * len2));
        }

        return 0.0f;
    }

    template<class T>
    inline T Vector2<T>::AngleBetween(const Vector2<T>& v2) const {     // computes the angle between this vector and vect
        return acos(CosAngleBetween(v2));
    }

    template<class T>
    inline T Vector2<T>::DistanceSquared(const Vector2<T>& v2) const { // finds the square of the distance between this and v2
        T dx = v2.x - x;
        T dy = v2.y - y;
        return dx * dx + dy * dy;
    }

    template<class T>
    inline T Vector2<T>::Distance(const Vector2<T>& v2) const {          // finds the distance between this and v2
        return sqrt(DistanceSquared(v2));
    }

    template<class T>
    inline Vector2<T>& Vector2<T>::PerpendicularOf(const Vector2<T> &v2) {            // sets this vector to the 2D perpendicular vector to v
        x = -v2.y;
        y =  v2.x;

        return *this;
    }

    template<class T>
    inline Vector2<T> Vector2<T>::Perpendicular(void) const {                     // returns the vector perpendicular to this
        return Vector2<T>(-y, x);
    }

    template<class T>
    inline bool Vector2<T>::operator == (const Vector2<T>& rhs) const {          // compares this to rhs and returns if they are equal
        return (this == &rhs || (x == rhs.x && y == rhs.y));
    }

    template<class T>
    inline bool Vector2<T>::operator != (const Vector2<T>& rhs) const {          // compares this to rhs and returns if they are not equal
        return !((*this) == rhs);
    }
}