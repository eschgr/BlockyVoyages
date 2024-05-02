#include <math.h>

namespace BlockyVoyages {
    template<class T> inline Vector4<T>::Vector4(void) {}

    template<class T> inline Vector4<T>::Vector4(const Vector4<T>& other)
        : x(other.x), y(other.y), z(other.z), w(other.w)
    {}

    template<class T> inline Vector4<T>::Vector4(const Vector2<T>& other, T _z, T _w)
        : x(other.x), y(other.y), z(_z), w(_w)
    {}

    template<class T> inline Vector4<T>::Vector4(const Vector2<T>& other, T _z)
        : x(other.x), y(other.y), z(_z), w(0)
    {}

    template<class T> inline Vector4<T>::Vector4(const Vector3<T>& other, T _w)
        : x(other.x), y(other.y), z(other.z), w(_w)
    {}

    template<class T> inline Vector4<T>::Vector4(const Vector3<T>& other)
        : x(other.x), y(other.y), z(other.z), w(0)
    {}

    // set from floats
    template<class T> inline Vector4<T>::Vector4(T _x, T _y, T _z, T _w)
        : x(_x), y(_y), z(_z), w(_w)
    {}

    template<class T> inline Vector4<T>::Vector4(T _x, T _y, T _z)
        : x(_x), y(_y), z(_z), w(0)
    {}

    template<class T> inline Vector4<T>::Vector4(const T _v[4])
        : x(_v[0]), y(_v[1]), z(_v[2]), w(_v[3])
    {}

    template<class T> inline Vector4<T>::~Vector4(void) {}

    // setting methods
    template<class T> inline Vector4<T>& Vector4<T>::operator=(const Vector4<T>& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
            z = other.z;
            w = other.w;
        }

        return *this;
    }

    template<class T> inline void Vector4<T>::Set(const Vector2<T>& other, T _z, T _w) {
        x = other.x;
        y = other.y;
        z = _z;
        w = _w;
    }

    template<class T> inline void Vector4<T>::Set(const Vector2<T>& other, T _z) {
        x = other.x;
        y = other.y;
        z = _z;
        w = 0;
    }

    template<class T> inline void Vector4<T>::Set(const Vector3<T>& other, T _w) {
        x = other.x;
        y = other.y;
        z = other.z;
        w = _w;
    }

    template<class T> inline void Vector4<T>::Set(const Vector3<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        w = 0;
    }

    template<class T> inline void Vector4<T>::Set(const Vector4<T>& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        w = other.w;
    }

    // set from floats
    template<class T> inline void Vector4<T>::Set(T _x, T _y, T _z, T _w) {
        x = _x;
        y = _y;
        z = _z;
        w = _w;
    }

    template<class T> inline void Vector4<T>::Set(T _x, T _y, T _z) {
        x = _x;
        y = _y;
        z = _z;
        w = 0;
    }

    template<class T> inline void Vector4<T>::Set(const T _v[4]) {
        x = _v[0];
        y = _v[1];
        z = _v[2];
        w = _v[3];
    }

    // clean up and fast implementation
    template<class T> inline void Vector4<T>::Zero(void) {
        x = 0;
        y = 0;
        z = 0;
        w = 0;
    }

    template<class T> inline void Vector4<T>::Clean(void) {
        if (abs(x) <EPSILON) {
            x = 0;
        }

        if (abs(y) <EPSILON) {
            y = 0;
        }

        if (abs(z) <EPSILON) {
            z = 0;
        }

        if (abs(w) <EPSILON) {
            w = 0;
        }
    }

    template<class T> inline Vector4<T>& Vector4<T>::Add(const Vector4<T>&  lhs, const Vector4<T>& rhs) {
        x = lhs.x + rhs.x;
        y = lhs.y + rhs.y;
        z = lhs.z + rhs.z;
        w = lhs.w + rhs.w;
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Sub(const Vector4<T>&  lhs, const Vector4<T>& rhs) {
        x = lhs.x - rhs.x;
        y = lhs.y - rhs.y;
        z = lhs.z - rhs.z;
        w = lhs.w - rhs.w;
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Mul(const Vector4<T>&  lhs, T rhs) {
        x = lhs.x * rhs;
        y = lhs.y * rhs;
        z = lhs.z * rhs;
        w = lhs.w * rhs;
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Mul(const Vector4<T>&  lhs, const Matrix44<T>& rhs) {
        x = lhs.x * rhs.m00 + lhs.y * rhs.m10 + lhs.z * rhs.m20 + lhs.w * rhs.m30;
        y = lhs.x * rhs.m01 + lhs.y * rhs.m11 + lhs.z * rhs.m21 + lhs.w * rhs.m31;
        z = lhs.x * rhs.m02 + lhs.y * rhs.m12 + lhs.z * rhs.m22 + lhs.w * rhs.m32;
        w = lhs.x * rhs.m03 + lhs.y * rhs.m13 + lhs.z * rhs.m23 + lhs.w * rhs.m33;

        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::Mul(const Vector4<float32>&  lhs, const Matrix44<float32>& rhs) {
        const float32* vin = lhs.pntr;
        const float32* matin = rhs.pntr;
                float32* vout = pntr;

        __asm {
            // move the data address into the main memory
            mov        ecx,    vin
            mov        edx,    matin

            movlps  xmm6, qword ptr [ecx]
            movlps  xmm0, qword ptr [edx]
            shufps  xmm6, xmm6, 0x44
            movhps  xmm0, qword ptr [edx+16]
            mulps    xmm0, xmm6
            movlps  xmm7, qword ptr [ecx+ 8]
            movlps  xmm2, qword ptr [edx+ 8]
            shufps  xmm7, xmm7, 0x44
            movhps  xmm2, qword ptr [edx+24]
            mulps    xmm2, xmm7
            movlps  xmm1, qword ptr [edx+32]
            movhps  xmm1, qword ptr [edx+48]
            mulps    xmm1, xmm6
            movlps  xmm3, qword ptr [edx+40]
            addps    xmm0, xmm2
            movhps  xmm3, qword ptr [edx+56]

            mov        eax,    vout

            mulps    xmm3, xmm7
            movaps  xmm4, xmm0
            addps    xmm1, xmm3
            shufps  xmm4, xmm1, 0x88
            shufps  xmm0, xmm1, 0xDD
            addps    xmm0, xmm4

            // store the answer into this
            movups    [eax],    xmm0
        };

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Mul(const Matrix44<T>& lhs, const Vector4<T>& rhs) {
        x = rhs.x * lhs.m00 + rhs.y * lhs.m01 + rhs.z * lhs.m02 + rhs.w * lhs.m03;
        y = rhs.x * lhs.m10 + rhs.y * lhs.m11 + rhs.z * lhs.m12 + rhs.w * lhs.m13;
        z = rhs.x * lhs.m20 + rhs.y * lhs.m21 + rhs.z * lhs.m22 + rhs.w * lhs.m23;
        w = rhs.x * lhs.m30 + rhs.y * lhs.m31 + rhs.z * lhs.m32 + rhs.w * lhs.m33;

        return *this;
    }

    //void MatrixMultiply3(Matrix4f &m, Vector4f *vin, Vector4f *vout) {
    template<> inline Vector4<float32>& Vector4<float32>::Mul(const Matrix44<float32>& lhs, const Vector4<float32>& rhs) {
        // Get a pointer to the elements of m
        const float *row0 = lhs.pntr;
        const float *vin = rhs.pntr;
                float* vout = this->pntr;

        __asm {
            mov            esi, vin         // in
            mov            edi, vout        // out

            // load columns of matrix into xmm4-7
            mov            edx, row0
            movups    xmm4, [edx]
            movups    xmm5, [edx+0x10]
            movups    xmm6, [edx+0x20]
            movups    xmm7, [edx+0x30]

            // load v into xmm0.
            movups    xmm0, [esi]

            // we'll store the final result in xmm2; initialize it
            // to zero
            xorps        xmm2, xmm2

            // broadcast x into xmm1, multiply it by the first
            // column of the matrix (xmm4), and add it to the total
            movups    xmm1, xmm0
            shufps    xmm1, xmm1, 0x00
            mulps        xmm1, xmm4
            addps        xmm2, xmm1

            // repeat the process for y, z and w
            movups    xmm1, xmm0
            shufps    xmm1, xmm1, 0x55
            mulps        xmm1, xmm5
            addps        xmm2, xmm1
            movups    xmm1, xmm0
            shufps    xmm1, xmm1, 0xAA
            mulps        xmm1, xmm6
            addps        xmm2, xmm1
            movups    xmm1, xmm0
            shufps    xmm1, xmm1, 0xFF
            mulps        xmm1, xmm7
            addps        xmm2, xmm1

            // write the results to vout
            movups    [edi], xmm2
        }
    }

    template<class T> inline Vector4<T>& Vector4<T>::Div(const Vector4<T>&  lhs, T rhs) {
        x = lhs.x / rhs;
        y = lhs.y / rhs;
        z = lhs.z / rhs;
        w = lhs.w / rhs;
        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::Div(const Vector4<float32>&  lhs, float32 rhs) {
        rhs = 1.0f / rhs;
        x = lhs.x * rhs;
        y = lhs.y * rhs;
        z = lhs.z * rhs;
        w = lhs.w * rhs;
        return *this;
    }

    template<> inline Vector4<float64>& Vector4<float64>::Div(const Vector4<float64>&  lhs, float64 rhs) {
        rhs = 1.0 / rhs;
        x = lhs.x * rhs;
        y = lhs.y * rhs;
        z = lhs.z * rhs;
        w = lhs.w * rhs;
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Build(const Vector4<T>& init, const Vector4<T>& term) {
        return Sub(term, init);
    }


    template<class T> inline T Vector4<T>::Dot(const Vector4<T>& rhs) const {                 // computes the dot product between this and v2
        return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Cross(const Vector4<T>& lhs, const Vector4<T>& rhs) { // sets this vector to the cross product of v1 and v2
        x = lhs.y * rhs.z - lhs.z * rhs.y;
        y = lhs.z * rhs.x - lhs.x * rhs.z;
        z = lhs.x * rhs.y - lhs.y * rhs.x;
        w = 0;

        return *this;
    }

    template<class T> inline const Vector4<T> Vector4<T>::Cross(const Vector4<T>& rhs) const {                     // returns the cross product of this and v2
        return Vector4<T>(y * rhs.z - z * rhs.y,
                             z * rhs.x - x * rhs.z,
                             x * rhs.y - y * rhs.x,
                             0);
    }

    template<class T> inline const Vector4<T> Vector4<T>::operator + (const Vector4<T> &rhs) const {        // Adds two vectors. use Add() or += instead if possible
        return Vector4<T>(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
    }

    template<class T> inline const Vector4<T> Vector4<T>::operator - (const Vector4<T> &rhs) const {        // subtracts two vectors. use Sub() or -= instead if possible
        return Vector4<T>(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
    }

    template<class T> inline         T Vector4<T>::operator * (const Vector4<T> &rhs) const {        // computes the dot product between this and rhs
        return x * rhs.x + y * rhs.y + z * rhs.z + w * rhs.w;
    }

    template<class T> inline const Vector4<T> Vector4<T>::operator * (T rhs) const {             // multiplies this by a scalar. use Mul() or *= instead if possible
        return Vector4<T>(x * rhs, y * rhs, z * rhs, w * rhs);
    }

    template<class T> inline const Vector4<T> Vector4<T>::operator * (const Matrix44<T>& rhs) const {
        return Vector4<T>(x * rhs.m00 + y * rhs.m10 + z * rhs.m20 + w * rhs.m30,
                                    x * rhs.m01 + y * rhs.m11 + z * rhs.m21 + w * rhs.m31,
                                    x * rhs.m02 + y * rhs.m12 + z * rhs.m22 + w * rhs.m32,
                                    x * rhs.m03 + y * rhs.m13 + z * rhs.m23 + w * rhs.m33);
    }

    template<> inline const Vector4<float32> Vector4<float32>::operator * (const Matrix44<float32>& rhs) const {
        Vector4<float32> ans;
                float32* vout = ans.pntr;
        const float32* vin = pntr;
        const float32* matIn = rhs.pntr;
        __asm {
            // move the data address into the main memory
            mov        ecx,    vin
            mov        edx,    matIn

            movlps  xmm6, qword ptr [ecx]
            movlps  xmm0, qword ptr [edx]
            shufps  xmm6, xmm6, 0x44
            movhps  xmm0, qword ptr [edx+16]
            mulps    xmm0, xmm6
            movlps  xmm7, qword ptr [ecx+ 8]
            movlps  xmm2, qword ptr [edx+ 8]
            shufps  xmm7, xmm7, 0x44
            movhps  xmm2, qword ptr [edx+24]
            mulps    xmm2, xmm7
            movlps  xmm1, qword ptr [edx+32]
            movhps  xmm1, qword ptr [edx+48]
            mulps    xmm1, xmm6
            movlps  xmm3, qword ptr [edx+40]
            addps    xmm0, xmm2
            movhps  xmm3, qword ptr [edx+56]

            mov        eax,    vout

            mulps    xmm3, xmm7
            movaps  xmm4, xmm0
            addps    xmm1, xmm3
            shufps  xmm4, xmm1, 0x88
            shufps  xmm0, xmm1, 0xDD
            addps    xmm0, xmm4

            //movaps xmmword ptr [eax], xmm0

            // store the answer into this
            movups    [eax],    xmm0
        };
        return ans;
    }

    template<class T> inline const Vector4<T> Vector4<T>::operator / (T rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        return Vector4<T>(x / rhs, y / rhs, z / rhs, w / rhs);
    }

    template<> inline const Vector4<float32> Vector4<float32>::operator / (float32 rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        rhs = 1.0f / rhs;
        return Vector4<float32>(x * rhs, y * rhs, z * rhs, w * rhs);
    }

    template<> inline const Vector4<float64> Vector4<float64>::operator / (float64 rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        rhs = 1.0 / rhs;
        return Vector4<float64>(x * rhs, y * rhs, z * rhs, w * rhs);
    }

    template<class T> inline Vector4<T>& Vector4<T>::operator += (const Vector4<T> &rhs) {     // does an in-place addition with rhs
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        w += rhs.w;

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::operator -= (const Vector4<T> &rhs) {     // does an in-place subtraction with rhs
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        w -= rhs.w;

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::operator *= (T rhs) {          // does an in-place scalar multiplication
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::operator *= (const Matrix44<T>& rhs) {          // does an in-place matrix multiplication
        T tx = x * rhs.m00 + y * rhs.m10 + z * rhs.m20 + w * rhs.m30;
        T ty = x * rhs.m01 + y * rhs.m11 + z * rhs.m21 + w * rhs.m31;
        T tz = x * rhs.m02 + y * rhs.m12 + z * rhs.m22 + w * rhs.m32;
        w = x * rhs.m03 + y * rhs.m13 + z * rhs.m23 + w * rhs.m33;

        x = tx;
        y = ty;
        z = tz;

        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::operator *= (const Matrix44<float32>& rhs) {          // does an in-place matrix multiplication
        Vector4<float32> ans;
                float32* vout = ans.pntr;
                float32* vin = pntr;
        const float32* matIn = rhs.pntr;
        __asm {
            mov        ecx,    vin
            mov        edx,    matIn
            mov        eax,    vout

            movlps  xmm6, qword ptr [ecx]
            movlps  xmm0, qword ptr [edx]
            shufps  xmm6, xmm6, 0x44
            movhps  xmm0, qword ptr [edx+16]
            mulps    xmm0, xmm6
            movlps  xmm7, qword ptr [ecx+ 8]
            movlps  xmm2, qword ptr [edx+ 8]
            shufps  xmm7, xmm7, 0x44
            movhps  xmm2, qword ptr [edx+24]
            mulps    xmm2, xmm7
            movlps  xmm1, qword ptr [edx+32]
            movhps  xmm1, qword ptr [edx+48]
            mulps    xmm1, xmm6
            movlps  xmm3, qword ptr [edx+40]
            addps    xmm0, xmm2
            movhps  xmm3, qword ptr [edx+56]

            mulps    xmm3, xmm7
            movaps  xmm4, xmm0
            addps    xmm1, xmm3
            shufps  xmm4, xmm1, 0x88
            shufps  xmm0, xmm1, 0xDD
            addps    xmm0, xmm4

            movups    [eax],    xmm0
        };

        memcpy(this, ans.pntr, sizeof(Vector4<float32>));
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::operator /= (T rhs) {          // does an in-place scalar division
        x /= rhs;
        y /= rhs;
        z /= rhs;
        w /= rhs;

        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::operator /= (float32 rhs) {          // does an in-place scalar division
        rhs = 1.0f / rhs;
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;

        return *this;
    }

    template<> inline Vector4<float64>& Vector4<float64>::operator /= (float64 rhs) {          // does an in-place scalar division
        rhs = 1.0 / rhs;
        x *= rhs;
        y *= rhs;
        z *= rhs;
        w *= rhs;

        return *this;
    }

    template<class T> inline T Vector4<T>::LengthSquared(void) const {
        return x * x + y * y + z * z + w * w;
    }

    template<class T> inline T Vector4<T>::Length(void) const {
        return sqrt(LengthSquared());
    }

    template<class T> inline Vector4<T>& Vector4<T>::NegationOf(const Vector4<T>& other) {
        x = -other.x;
        y = -other.y;
        z = -other.z;
        w = -other.w;

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Negate(void) {
        x *= -1;
        y *= -1;
        z *= -1;
        w *= -1;

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::UnitOf(const Vector4<T>& other) {
        Set(other);
        Normalize();
        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::Normalize(void) {
        T len = x * x + y * y + z * z + w * w;
        if (abs(len)> TINY) {
            len = sqrt(len);
            x /= len;
            y /= len;
            z /= len;
            w /= len;
        }

        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::Normalize(void) {
        float32 len = x * x + y * y + z * z + w * w;
        if (abs(len)> TINY) {
            len = 1.0f / sqrt(len);
            x *= len;
            y *= len;
            z *= len;
            w *= len;
        }

        return *this;
    }

    template<> inline Vector4<float64>& Vector4<float64>::Normalize(void) {
        float64 len = x * x + y * y + z * z + w * w;
        if (abs(len)> TINY) {
            len = 1.0 / sqrt(len);
            x *= len;
            y *= len;
            z *= len;
            w *= len;
        }

        return *this;
    }

    template<class T> inline const Matrix44<T> Vector4<T>::Tensor(const Vector4<T> &rhs) const {  // computes the outer product of this and rhs
        return Matrix44<T>(x * rhs.x, x * rhs.y, x * rhs.z, x * rhs.w,
                              y * rhs.x, y * rhs.y, y * rhs.z, y * rhs.w,
                              z * rhs.x, z * rhs.y, z * rhs.z, z * rhs.w,
                              w * rhs.x, w * rhs.y, w * rhs.z, w * rhs.w);
    }

    template<class T> inline void Vector4<T>::GetOrthonormalBasis(Vector4<T>& right, Vector4<T>& up) const { // makes u and v into an orthonormal basis of this. Assumes this is unit length
        right.x = x;
        right.y = y;
        right.z = z;
        if (z <x && z <y) {
            right.z += 1;
        }
        else if (y <x && y <z) {
            right.y += 1;
        }
        else {
            right.x += 1;
        }
         
        up.Cross(*this, right);
        right.Cross(*this, up);

        right.Normalize();
        up.Normalize();
    }

    template<class T> inline T Vector4<T>::CosAngleBetween(const Vector4<T>& vect) const { // computes the cosine of the angle between this vector and vect
        T len1 = vect.LengthSquared();
        T len2 = LengthSquared();

        if (len1> TINY && len2> TINY) {
            return Dot(vect) / sqrt(len1 * len2);
        }
        
        return 0;
    }

    template<class T> inline T Vector4<T>::AngleBetween(const Vector4<T>& vect) const {     // computes the angle between this vector and vect
        return acos(CosAngleBetween(vect));
    }

    template<class T> inline T Vector4<T>::DistanceSquared(const Vector4<T>& v2) const { // finds the square of the distance between this and v2
        T dx = v2.x - x;
        T dy = v2.y - y;
        T dz = v2.z - z;
        T dw = v2.w - w;
        return dx * dx + dy * dy + dz * dz + dw * dw;
    }

    template<class T> inline T Vector4<T>::Distance(const Vector4<T>& v2) const {          // finds the distance between this and v2
        return sqrt(DistanceSquared(v2));
    }

    template<class T> inline Vector4<T>& Vector4<T>::Homogenize() {
        if (abs(w)> EPSILON) {
            x /= inv;
            y /= inv;
            z /= inv;
            w = 1;
        }
        else {
            w = 0;
        }

        return *this;
    }

    template<> inline Vector4<float32>& Vector4<float32>::Homogenize() {
        if (abs(w)> EPSILON) {
            float32 inv = 1.0f / w;
            x *= inv;
            y *= inv;
            z *= inv;
            w = 1;
        }
        else {
            w = 0;
        }

        return *this;
    }

    template<> inline Vector4<float64>& Vector4<float64>::Homogenize() {
        if (abs(w)> EPSILON) {
            float64 inv = 1.0 / w;
            x *= inv;
            y *= inv;
            z *= inv;
            w = 1;
        }
        else {
            w = 0;
        }

        return *this;
    }

    template<class T> inline Vector4<T>& Vector4<T>::HomogenizedOf(const Vector4<T>& v2) {
        Set(v2);
        Homogenize();

        return *this;
    }

    template<class T> inline bool Vector4<T>::operator == (const Vector4<T>& rhs) const {          // compares this to rhs and returns if they are equal
        return (this == &rhs) || (x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w);
    }

    template<class T> inline bool Vector4<T>::operator != (const Vector4<T>& rhs) const {          // compares this to rhs and returns if they are not equal
        return !((*this) == rhs);
    }
}