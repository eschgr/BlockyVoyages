#include <math.h>

namespace BlockyVoyages {
    // methods
    // constructors
    template<class T> 
    inline Vector3<T>::Vector3(void) {}
    
    template<class T> inline Vector3<T>::Vector3(const Vector2<T>& other, T _z) 
        : x(other.x), y(other.y), z(_z)
    {}

    template<class T> inline Vector3<T>::Vector3(const Vector3<T>& other)
        : x(other.x), y(other.y), z(other.z)
    {}

    template<class T> inline Vector3<T>::Vector3(const Vector4<T>& other)
        : x(other.x), y(other.y), z(other.z)
    {}

    // set from floats
    template<class T> inline Vector3<T>::Vector3(T _x, T _y, T _z)
        : x(_x), y(_y), z(_z)
    {}

    template<class T> inline Vector3<T>::Vector3(const T _v[3])
        : x(_v[0]), y(_v[1]), z(_v[2])
    {}

    template<class T> inline Vector3<T>::~Vector3(void) 
    {}

    // setting methods
    template<class T> inline Vector3<T>& Vector3<T>::operator=(const Vector3<T>& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
            z = other.z;
        }

        return *this;
    }

    template<class T> inline void Vector3<T>::Set(const Vector2<T>& other, T _z) {
        x = other.x;
        y = other.y;
        z = _z;
    }

    template<class T> inline void Vector3<T>::Set(const Vector3<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }

    template<class T> inline void Vector3<T>::Set(const Vector4<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }

    // set from floats
    template<class T> inline void Vector3<T>::Set(T _x, T _y, T _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    template<class T> inline void Vector3<T>::Set(const T _v[3]) {
        x = _v[0];
        y = _v[1];
        z = _v[2];
    }

    // clean up and fast implementation
    template<class T> inline void Vector3<T>::Zero(void) {
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
    }

    template<class T> inline void Vector3<T>::Clean(void) {
        if (abs(x) <EPSILON) {
            x = EPSILON;
        }

        if (abs(y) <EPSILON) {
            y = EPSILON;
        }

        if (abs(z) <EPSILON) {
            z = EPSILON;
        }
    }

    template<class T> inline Vector3<T>& Vector3<T>::Add(const Vector3<T>&  lhs, const Vector3<T>& rhs) {
        x = lhs.x + rhs.x;
        y = lhs.y + rhs.y;
        z = lhs.z + rhs.z;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Sub(const Vector3<T>&  lhs, const Vector3<T>& rhs) {
        x = lhs.x - rhs.x;
        y = lhs.y - rhs.y;
        z = lhs.z - rhs.z;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Mul(const Vector3<T>&  lhs, T rhs) {
        x = lhs.x * rhs;
        y = lhs.y * rhs;
        z = lhs.z * rhs;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Mul(const Vector3<T>&  lhs, const Matrix33<T>& rhs) {
        x = lhs.x * rhs.m00 + lhs.y * rhs.m10 + lhs.z * rhs.m20;
        y = lhs.x * rhs.m01 + lhs.y * rhs.m11 + lhs.z * rhs.m21;
        z = lhs.x * rhs.m02 + lhs.y * rhs.m12 + lhs.z * rhs.m22;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Mul(const Matrix33<T>& lhs, const Vector3<T>& rhs) {
        x = lhs.m00 * rhs.x + lhs.m01 * rhs.y + lhs.m02 * rhs.z;
        y = lhs.m10 * rhs.x + lhs.m11 * rhs.y + lhs.m12 * rhs.z;
        z = lhs.m20 * rhs.x + lhs.m21 * rhs.y + lhs.m22 * rhs.z;
        return *this;
    }

    template<> inline Vector3<float32>& Vector3<float32>::Div(const Vector3<float32>&  lhs, float32 rhs) {
        float32 inv = 1.0f / rhs;
        x = lhs.x * inv;
        y = lhs.y * inv;
        z = lhs.z * inv;
        return *this;
    }

    template<> inline Vector3<float64>& Vector3<float64>::Div(const Vector3<float64>&  lhs, float64 rhs) {
        float64 inv = 1.0 / rhs;
        x = lhs.x * inv;
        y = lhs.y * inv;
        z = lhs.z * inv;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Div(const Vector3<T>&  lhs, T rhs) {
        x = lhs.x / rhs;
        y = lhs.y / rhs;
        z = lhs.z / rhs;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Build(const Vector3<T>& init, const Vector3<T>& term) {
        return Sub(term, init);
    }

    template<class T> inline T Vector3<T>::Dot(const Vector3<T>& rhs) const {                 // computes the dot product between this and v2
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Cross(const Vector3<T>& lhs, const Vector3<T>& rhs) { // sets this vector to the cross product of v1 and v2
        x = lhs.y * rhs.z - lhs.z * rhs.y;
        y = lhs.z * rhs.x - lhs.x * rhs.z;
        z = lhs.x * rhs.y - lhs.y * rhs.x;
        return *this;
    }

    template<class T> inline const Vector3<T> Vector3<T>::Cross(const Vector3<T>& rhs) const {                     // returns the cross product of this and v2
        return Vector3<T>(y * rhs.z - z * rhs.y,
                             z * rhs.x - x * rhs.z,
                             x * rhs.y - y * rhs.x);
    }

    template<class T> inline const Vector3<T> Vector3<T>::operator + (const Vector3<T> &rhs) const {        // Adds two vectors. use Add() or += instead if possible
        return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    template<class T> inline const Vector3<T> Vector3<T>::operator - (const Vector3<T> &rhs) const {        // subtracts two vectors. use Sub() or -= instead if possible
        return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    template<class T> inline         T Vector3<T>::operator * (const Vector3<T> &rhs) const {        // computes the dot product between this and rhs
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    template<class T> inline const Vector3<T> Vector3<T>::operator * (T rhs) const {             // multiplies this by a scalar. use Mul() or *= instead if possible
        return Vector3<T>(x * rhs, y * rhs, z * rhs);
    }

    template<class T> inline const Vector3<T> Vector3<T>::operator * (const Matrix33<T>& rhs) const {
        return Vector3<T>(x * rhs.m00 + y * rhs.m10 + z * rhs.m20,
                                    x * rhs.m01 + y * rhs.m11 + z * rhs.m21,
                                    x * rhs.m02 + y * rhs.m12 + z * rhs.m22);
    }
    template<class T> inline const Vector3<T> Vector3<T>::operator / (T rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        return Vector3<T>(x / rhs, y / rhs, z / rhs);
    }

    template<> inline const Vector3<float32> Vector3<float32>::operator / (float32 rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        rhs = 1.0f / rhs;
        return Vector3<float32>(x * rhs, y * rhs, z * rhs);
    }

    template<> inline const Vector3<float64> Vector3<float64>::operator / (float64 rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        rhs = 1.0 / rhs;
        return Vector3<float64>(x * rhs, y * rhs, z * rhs);
    }

    template<class T> inline Vector3<T>& Vector3<T>::operator += (const Vector3<T> &rhs) {     // does an in-place addition with rhs
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::operator -= (const Vector3<T> &rhs) {     // does an in-place subtraction with rhs
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::operator *= (T rhs) {          // does an in-place scalar multiplication
        x *= rhs;
        y *= rhs;
        z *= rhs;

        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::operator *= (const Matrix33<T>& mat) {          // does an in-place matrix multiplication
        T tx = x * mat.m00 + y * mat.m10 + z * mat.m20;
        T ty = x * mat.m01 + y * mat.m11 + z * mat.m21;
        z = x * mat.m02 + y * mat.m12 + z * mat.m22;
        x = tx;
        y = ty;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::operator /= (T rhs) {          // does an in-place scalar division
        x /= rhs;
        y /= rhs;
        z /= rhs;
        return *this;
    }

    template<> inline Vector3<float32>& Vector3<float32>::operator /= (float32 rhs) {          // does an in-place scalar division
        rhs = 1.0f / rhs;
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    template<> inline Vector3<float64>& Vector3<float64>::operator /= (float64 rhs) {          // does an in-place scalar division
        rhs = 1.0 / rhs;
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    template<class T> inline T Vector3<T>::LengthSquared(void) const {
        return x * x + y * y + z * z;
    }

    template<class T> inline T Vector3<T>::Length(void) const {
        return sqrt(LengthSquared());
    }

    template<class T> inline Vector3<T>& Vector3<T>::NegationOf(const Vector3<T>& other) {
        x = -other.x;
        y = -other.y;
        z = -other.z;
        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::Negate(void) {
        x *= -1.0f;
        y *= -1.0f;
        z *= -1.0f;
        return *this;
    }

    template<> inline Vector3<float64>& Vector3<float64>::Normalize(void) {
        float64 len = LengthSquared();
        if (len> TINY) {
            len = sqrt(len);
            x /= len;
            y /= len;
            z /= len;
        }
        else {
            Zero();
        }

        return *this;
    }

    template<> inline Vector3<float32>& Vector3<float32>::Normalize(void) {
        float32 len = LengthSquared();
        if (len> TINY) {
            len = 1.0f / sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        else {
            Zero();
        }

        return *this;
    }

    template<class T> inline Vector3<T>& Vector3<T>::UnitOf(const Vector3<T>& other) {
        Set(other);
        Normalize();
    }

    template<class float64> inline Vector3<float64>& Vector3<float64>::Normalize(void) {
        float64 len = LengthSquared();
        if (len> TINY) {
            len = 1.0 / sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        else {
            Zero();
        }

        return *this;
    }

    template<class T> inline const Matrix33<T> Vector3<T>::Tensor(const Vector3<T> &rhs) const {  // computes the outer product of this and rhs
        return Matrix33<T>(x * rhs.x, x * rhs.y, x * rhs.z,
                                     y * rhs.x, y * rhs.y, y * rhs.z,
                                     z * rhs.x, z * rhs.y, z * rhs.z);
    }

    template<class T> inline void Vector3<T>::GetOrthonormalBasis(Vector3<T>& right, Vector3<T>& up) const { // makes u and v into an orthonormal basis of this. Assumes this is unit length
        right.x = x;
		right.y = y;
		right.z = z;
        if (z <x && z <y) {
            right.z += 1.0f;
        }
        else if (y <x && y <z) {
            right.y += 1.0f;
        }
        else {
            right.x += 1.0f;
        }
	     
        up.Cross(*this, right);
        right.Cross(*this, up);

        right.Normalize();
        up.Normalize();
    }

    template<class T> inline T Vector3<T>::CosAngleBetween(const Vector3<T>& v2) const { // computes the cosine of the angle between this vector and vect
        T len1 = LengthSquared();
        T len2 = v2.LengthSquared();
        if (len1> TINY && len2> TINY) {
            return Dot(v2) / sqrt(len1 * len2);
        }
        return 0.0f;
    }

    template<class T> inline T Vector3<T>::AngleBetween(const Vector3<T>& v2) const {     // computes the angle between this vector and vect
        return acos(CosAngleBetween(v2));
    }

    template<class T> inline T Vector3<T>::DistanceSquared(const Vector3<T>& v2) const { // finds the square of the distance between this and v2
        T dx = x - v2.x;
        T dy = y - v2.y;
        T dz = z - v2.z;
        return dx * dx + dy * dy + dz * dz;
    }

    template<class T> inline T Vector3<T>::Distance(const Vector3<T>& v2) const {          // finds the distance between this and v2
        return sqrt(DistanceSquared(v2));
    }

    template<class T> inline bool Vector3<T>::operator == (const Vector3<T>& rhs) const {          // compares this to rhs and returns if they are equal
        return (this == &rhs) || (x == rhs.x && y == rhs.y && z == rhs.z);
    }

    template<class T> inline bool Vector3<T>::operator != (const Vector3<T>& rhs) const {          // compares this to rhs and returns if they are not equal
        return !(*this == rhs);
    }
};