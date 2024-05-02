#include <math.h>

namespace BlockyVoyages {
    // Constructors
    template<class T> inline Matrix44<T>::Matrix44(void) {}

    template<class T> inline Matrix44<T>::Matrix44(const Matrix44<T>& other)
        : m00(other.m00), m10(other.m10), m20(other.m20), m30(other.m30), 
          m01(other.m01), m11(other.m11), m21(other.m21), m31(other.m31), 
          m02(other.m02), m12(other.m12), m22(other.m22), m32(other.m32), 
          m03(other.m03), m13(other.m13), m23(other.m23), m33(other.m33)
    {}

    template<class T> inline Matrix44<T>::Matrix44(T v[16]) 
        : m00(v[ 0]), m10(v[ 1]), m20(v[ 2]), m30(v[ 3]), 
          m01(v[ 4]), m11(v[ 5]), m21(v[ 6]), m31(v[ 7]), 
          m02(v[ 8]), m12(v[ 9]), m22(v[10]), m32(v[11]), 
          m03(v[12]), m13(v[13]), m23(v[14]), m33(v[15])
    {}

    template<class T> inline Matrix44<T>::Matrix44(T _00, T _01, T _02, T _03,
                                                                         T _10, T _11, T _12, T _13, 
                                                                         T _20, T _21, T _22, T _23, 
                                                                         T _30, T _31, T _32, T _33) 
        : m00(_00), m10(_10), m20(_20), m30(_30), 
          m01(_01), m11(_11), m21(_21), m31(_31), 
          m02(_02), m12(_12), m22(_22), m32(_32), 
          m03(_03), m13(_13), m23(_23), m33(_33)
    {}

    template<class T> inline Matrix44<T>::Matrix44(const Matrix22<T>& other)
        : m00(other.m00), m10(other.m10), m20(0.0f), m30(0.0f), 
          m01(other.m01), m11(other.m11), m21(0.0f), m31(0.0f), 
          m02(0.0f),        m12(0.0f),        m22(0.0f), m32(0.0f), 
          m03(0.0f),        m13(0.0f),        m23(0.0f), m33(0.0f)
    {}

    template<class T> inline Matrix44<T>::Matrix44(const Matrix33<T>& other)
        : m00(other.m00), m10(other.m10), m20(other.m20), m30(0.0f), 
          m01(other.m01), m11(other.m11), m21(other.m21), m31(0.0f), 
          m02(other.m02), m12(other.m12), m22(other.m22), m32(0.0f), 
          m03(0.0f),        m13(0.0f),        m23(0.0f),        m33(1.0f)
    {}

    template<class T> inline Matrix44<T>::Matrix44(const Quaternion<T>& q) {
        T s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;

        s = 2.0f / (q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
        xs = s * q.x;    ys = s * q.y;    zs = s * q.z;
        wx = q.w * xs;    wy = q.w * ys;    wz = q.w * zs;
        xx = q.x * xs;    xy = q.x * ys;    xz = q.x * zs;
        yy = q.y * ys;    yz = q.y * zs;    zz = q.z * zs;
            
        m00 = 1.0f - (yy + zz); m10 = xy + wz;             m20 = xz - wy;
        m01 = xy - wz;             m11 = 1.0f - (xx + zz); m21 = yz + wx;
        m02 = xz + wy;             m12 = yz - wx;             m22 = 1.0f - (xx + yy);
    }

    template<> inline Matrix44<float32>::Matrix44(const Quaternion<float32>& q) {
        Matrix44<float32> mat1( q.w, -q.z,  q.y, q.x,
                                            q.z,  q.w, -q.x, q.y,
                                          -q.y,  q.x,  q.w, q.z,
                                          -q.x, -q.y, -q.z, q.w);

        Matrix44<float32> mat2( q.w, -q.z,  q.y, -q.x,
                                            q.z,  q.w, -q.x, -q.y,
                                          -q.y,  q.x,  q.w, -q.z,
                                            q.x,  q.y,  q.z,  q.w);

        Mul(mat1, mat2);

        // do some cleanup. If not done, numerical inaccuracy causes problems
        m03 = m13 = m23 = 0;
        m30 = m31 = m32 = 0;
        m33 = 1;
    }

    // Destructor
    template<class T> inline Matrix44<T>::~Matrix44(void) {}

    // Assignment methods
    template<class T> inline Matrix44<T>& Matrix44<T>::operator= (const Matrix44<T>& other) {
        if (this != &other) {
            m00 = other.m00; m10 = other.m10; m20 = other.m20; m30 = other.m30; 
            m01 = other.m01; m11 = other.m11; m21 = other.m21; m31 = other.m31; 
            m02 = other.m02; m12 = other.m12; m22 = other.m22; m32 = other.m32; 
            m03 = other.m03; m13 = other.m13; m23 = other.m23; m33 = other.m33;
        }

        return *this;
    }

    template<class T> inline void Matrix44<T>::Set(T v[16]) {
        m00 = v[ 0]; m10 = v[ 1]; m20 = v[ 2]; m30 = v[ 3]; 
        m01 = v[ 4]; m11 = v[ 5]; m21 = v[ 6]; m31 = v[ 7]; 
        m02 = v[ 8]; m12 = v[ 9]; m22 = v[10]; m32 = v[11]; 
        m03 = v[12]; m13 = v[13]; m23 = v[14]; m33 = v[15];
    }

    template<class T> inline void Matrix44<T>::Set(T _00, T _01, T _02, T _03,
                                                                         T _10, T _11, T _12, T _13, 
                                                                         T _20, T _21, T _22, T _23, 
                                                                         T _30, T _31, T _32, T _33) 
    {
        m00 = _00; m01 = _01; m02 = _02; m03 = _03; 
        m10 = _10; m11 = _11; m12 = _12; m13 = _13; 
        m20 = _20; m21 = _21; m22 = _22; m23 = _23; 
        m30 = _30; m31 = _31; m32 = _32; m33 = _33;
    }

    template<class T> inline void Matrix44<T>::Set(const Matrix44<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = other.m20; m30 = other.m30; 
        m01 = other.m01; m11 = other.m11; m21 = other.m21; m31 = other.m31; 
        m02 = other.m02; m12 = other.m12; m22 = other.m22; m32 = other.m32; 
        m03 = other.m03; m13 = other.m13; m23 = other.m23; m33 = other.m33;
    }

    template<class T> inline void Matrix44<T>::Set(const Matrix33<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = other.m20; m30 = 0.0f;
        m01 = other.m01; m11 = other.m11; m21 = other.m21; m31 = 0.0f;
        m02 = other.m02; m12 = other.m12; m22 = other.m22; m32 = 0.0f;
        m03 = 0.0f;        m13 = 0.0f;        m23 = 0.0f;        m33 = 1.0f;
    }

    template<class T> inline void Matrix44<T>::Set(const Matrix22<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = 0.0f; m30 = 0.0f; 
        m01 = other.m01; m11 = other.m11; m21 = 0.0f; m31 = 0.0f; 
        m02 = 0.0f;        m12 = 0.0f;        m22 = 0.0f; m32 = 0.0f; 
        m03 = 0.0f;        m13 = 0.0f;        m23 = 0.0f; m33 = 0.0f;
    }

    template<class T> inline void Matrix44<T>::Set(const Quaternion<T>& q) {
        T s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;

        s = 2.0f / (q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
        xs = s * q.x;    ys = s * q.y;    zs = s * q.z;
        wx = q.w * xs;    wy = q.w * ys;    wz = q.w * zs;
        xx = q.x * xs;    xy = q.x * ys;    xz = q.x * zs;
        yy = q.y * ys;    yz = q.y * zs;    zz = q.z * zs;
            
        m00 = 1.0f - (yy + zz); m10 = xy + wz;             m20 = xz - wy;
        m01 = xy - wz;             m11 = 1.0f - (xx + zz); m21 = yz + wx;
        m02 = xz + wy;             m12 = yz - wx;             m22 = 1.0f - (xx + yy);
    }

    template<> inline void Matrix44<float32>::Set(const Quaternion<float32>& q) {
        Matrix44<float32> mat1( q.w, -q.z,  q.y, q.x,
                                    q.z,  q.w, -q.x, q.y,
                                  -q.y,  q.x,  q.w, q.z,
                                  -q.x, -q.y, -q.z, q.w);

        Matrix44<float32> mat2( q.w, -q.z,  q.y, -q.x,
                                    q.z,  q.w, -q.x, -q.y,
                                  -q.y,  q.x,  q.w, -q.z,
                                    q.x,  q.y,  q.z,  q.w);

        // do some cleanup. If not done, numerical inaccuracy causes problems
        m03 = m13 = m23 = 0;
        m30 = m31 = m32 = 0;
        m33 = 1;

        Mul(mat1, mat2);
    }

    // General utility
    // Cleans this matrix. Sets all close to zero values to zero
    template<class T> inline void Matrix44<T>::Clean(void) {
        for(int32 ind = 0; ind <16; ++ind) {
            if (IsCloseToZero(pntr[ind])) {
                pntr[ind] = 0.0f;
            }
        }
    }

    template<class T> inline void Matrix44<T>::Identity(void) {
        m00 = m11 = m22 = m33 = 1;
        m10 = m20 = m30 =
        m01 = m21 = m31 =
        m02 = m12 = m32 =
        m03 = m13 = m23 = 0;
    }

    template<class T> inline void Matrix44<T>::Zero(void) {
        m00 = m10 = m20 = m30 =
        m01 = m11 = m21 = m31 =
        m02 = m12 = m22 = m32 =
        m03 = m13 = m23 = m33 = 0;
    }

    // Gets the row specified.
    template<class T> inline const Vector4<T> Matrix44<T>::GetRow(int32 row) const {
        return Vector4<T>(pntr[row], pntr[row + 4], pntr[row + 8], pntr[row + 12]);
    }

    template<class T> inline void Matrix44<T>::GetRows(Vector4<T>& row0, Vector4<T>& row1, Vector4<T>& row2, Vector4<T>& row3) const {
        row0.Set(pntr[0], pntr[4], pntr[ 8], pntr[12]);
        row1.Set(pntr[1], pntr[5], pntr[ 9], pntr[13]);
        row2.Set(pntr[2], pntr[6], pntr[10], pntr[14]);
        row3.Set(pntr[3], pntr[7], pntr[11], pntr[15]);
    }

    // Gets the column specified
    template<class T> inline const Vector4<T> Matrix44<T>::GetColumn(int32 col) const {
        return Vector4<T>(&pntr[col * 4]);
    }

    template<class T> inline void  Matrix44<T>::GetColumns(Vector4<T>& col0, Vector4<T>& col1, Vector4<T>& col2, Vector4<T>& col3) const {
        col0.Set(&pntr[ 0]);
        col1.Set(&pntr[ 4]);
        col2.Set(&pntr[ 8]);
        col3.Set(&pntr[12]);
    }

    template<class T> inline void Matrix44<T>::SetRow(const Vector4<T>& v, int32 row) {
        pntr[row       ] = v.x;
        pntr[row +  4] = v.y;
        pntr[row +  8] = v.z;
        pntr[row + 12] = v.w;
    }

    template<class T> inline void Matrix44<T>::SetRow(T f1, T f2, T f3, T f4, int32 row) {
        pntr[row       ] = f1;
        pntr[row +  4] = f2;
        pntr[row +  8] = f3;
        pntr[row + 12] = f4;
    }

    template<class T> inline void Matrix44<T>::SetRows(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3) {
        m00 = row0.x; m10 = row1.x; m20 = row2.x; m30 = row3.x;
        m01 = row0.y; m11 = row1.y; m21 = row2.y; m31 = row3.y;
        m02 = row0.z; m12 = row1.z; m22 = row2.z; m32 = row3.z;
        m03 = row0.w; m13 = row1.w; m23 = row2.w; m33 = row3.w;
    }

    template<class T> inline void Matrix44<T>::SetColumn(const Vector4<T>& v, int32 col) {
        pntr[col * 3     ] = v.x; 
        pntr[col * 3 + 1] = v.y; 
        pntr[col * 3 + 2] = v.z; 
        pntr[col * 3 + 3] = v.w;
    }

    template<class T> inline void Matrix44<T>::SetColumn(T f1, T f2, T f3, T f4, int32 col) {
        pntr[col * 3     ] = f1; 
        pntr[col * 3 + 1] = f2; 
        pntr[col * 3 + 2] = f3; 
        pntr[col * 3 + 3] = f4;
    }

    template<class T> inline void Matrix44<T>::SetColumns(const Vector4<T>& col0, const Vector4<T>& col1, const Vector4<T>& col2, const Vector4<T>& col3) {
        m00 = col0.x; m10 = col0.y; m20 = col0.z; m30 = col0.w;
        m01 = col1.x; m11 = col1.y; m21 = col1.z; m31 = col1.w;
        m02 = col2.x; m12 = col2.y; m22 = col2.z; m32 = col2.w;
        m03 = col3.x; m13 = col3.y; m23 = col3.z; m33 = col3.w;
    }

    template<class T> inline void Matrix44<T>::Translation(T x, T y, T z) {
        m10 = m20 = m30 =
        m01 = m21 = m31 = 
        m02 = m12 = m32 = 0;
        m03 = x; m13 = y; m23 = z; 
        m00 = m11 = m22 = m33 = 1;
    }

    template<class T> inline void Matrix44<T>::Rotation(T angle, T x, T y, T z) {
        throw "Can't use method Matrix44<T>::Rotation with the requested type.";
    }

    template<> inline void Matrix44<float32>::Rotation(float32 angle, float32 x, float32 y, float32 z) {
        float32 fCos = cos(angle);
        float32 fSin = sin(angle);
        float32 fSum = 1.0f - fCos;

        float32 vectLen = x * x + y * y + z * z;

        if (IsNotZero(vectLen)) {
            vectLen = 1 / sqrt(vectLen);
            x *= vectLen;
            y *= vectLen;
            z *= vectLen;
        }

        m00 = (x * x) * fSum + fCos;
        m01 = (x * y) * fSum - (z * fSin);
        m02 = (x * z) * fSum + (y * fSin);

        m10 = (y * x) * fSum + (z * fSin);
        m11 = (y * y) * fSum + fCos;
        m12 = (y * z) * fSum - (x * fSin);

        m20 = (z * x) * fSum - (y * fSin);
        m21 = (z * y) * fSum + (x * fSin);
        m22 = (z * z) * fSum + fCos;
    }

    template<> inline void Matrix44<float64>::Rotation(float64 angle, float64 x, float64 y, float64 z) {
        float64 fCos = cos(angle);
        float64 fSin = sin(angle);
        float64 fSum = 1.0f - fCos;

        float64 vectLen = x * x + y * y + z * z;

        if (IsNotZero(vectLen)) {
            vectLen = 1 / sqrt(vectLen);
            x *= vectLen;
            y *= vectLen;
            z *= vectLen;
        }

        m00 = (x * x) * fSum + fCos;
        m01 = (x * y) * fSum - (z * fSin);
        m02 = (x * z) * fSum + (y * fSin);

        m10 = (y * x) * fSum + (z * fSin);
        m11 = (y * y) * fSum + fCos;
        m12 = (y * z) * fSum - (x * fSin);

        m20 = (z * x) * fSum - (y * fSin);
        m21 = (z * y) * fSum + (x * fSin);
        m22 = (z * z) * fSum + fCos;
    }

    template<class T> inline void Matrix44<T>::Scale(T x, T y, T z) {
        m00 = x; m10 = m20 = m30 = m01 = 0;
        m11 = y; m21 = m31 = m02 = m12 = 0;
        m22 = z; m32 = m03 = m13 = m23 = 0;
        m33 = 1;
    }

    template<class T> inline void Matrix44<T>::Shear(const Vector4<T>& axis, const Vector4<T>& shear) {
        m00 = 1 + normal.x * shear.x;
        m01 =      normal.y * shear.x;
        m02 =      normal.z * shear.x;
        m03 = 0;

        m10 =      normal.x * shear.y;
        m11 = 1 + normal.y * shear.y;
        m12 =      normal.z * shear.y;
        m13 = 0;

        m20 =      normal.x * shear.z;
        m21 =      normal.y * shear.z;
        m22 = 1 + normal.z * shear.z;
        m23 = m30 = m31 = m32 = 0;
        m33 = 1;
    }

    template<class T> inline void Matrix44<T>::LookAtGL(const Vector4<T>& from, const Vector4<T>& at, const Vector4<T>& up) {
        Vector4<T> viewDir;
        viewDir.Build(from, at);
        viewDir.Normalize();

        Vector4<T>viewRight;
        viewRight.Cross(viewDir, up);
        viewRight.Normalize();

        Vector4<T> viewUp;
        viewUp.Cross(viewRight, viewDir);
        
        m00 =  viewRight.x; m01 =  viewRight.y; m02 =  viewRight.z;
        m10 =  viewUp.x;     m11 =  viewUp.y;     m12 =  viewUp.z;
        m20 = -viewDir.x;    m21 = -viewDir.y;    m22 = -viewDir.z;

        m03 = m13 = m23 = m30 = m31 = m32 = 0;
        m33 = 1;

        Vector4<T> eyeInv(-((*this) * from));
        m03 = eyeInv.x;
        m13 = eyeInv.y;
        m23 = eyeInv.z;
    }

    template<class T> inline void Matrix44<T>::LookAtD3D(const Vector4<T>& from, const Vector4<T>& at, const Vector4<T>& up) {
        Vector4<T> viewDir;
        viewDir.Build(at, from);
        viewDir.Normalize();

        Vector4<T> viewUp(viewDir);

        viewUp *= (up * viewDir);
        viewUp.Sub(up, viewUp);
        viewUp.Normalize();

        Vector4<T> viewRight;
        viewRight.Cross(viewDir, viewUp);
        
        Zero();
        m00 = viewRight.x; m01 = viewRight.y; m02 = viewRight.z;
        m10 = viewUp.x;     m11 = viewUp.y;     m12 = viewUp.z;
        m20 = viewDir.x;    m21 = viewDir.y;    m22 = viewDir.z;

        Vector4<T> eyeInv(-((*this) * from));
        m03 = eyeInv.x;
        m13 = eyeInv.y;
        m23 = eyeInv.z;
        m33 = 1.0f;
    }

    template<class T> inline bool Matrix44<T>::InverseOf(const Matrix44<T>& mat) {
        T      fTemp[12];
        T      fDet;

        fTemp[0]  = mat.m22 * mat.m33;
        fTemp[1]  = mat.m32 * mat.m23;
        fTemp[2]  = mat.m12 * mat.m33;
        fTemp[3]  = mat.m32 * mat.m13;
        fTemp[4]  = mat.m12 * mat.m23;
        fTemp[5]  = mat.m22 * mat.m13;
        fTemp[6]  = mat.m02 * mat.m33;
        fTemp[7]  = mat.m32 * mat.m03;
        fTemp[8]  = mat.m02 * mat.m23;
        fTemp[9]  = mat.m22 * mat.m03;
        fTemp[10]  = mat.m02 * mat.m13;
        fTemp[11]  = mat.m12 * mat.m03;

        m00  = fTemp[0]*mat.m11 + fTemp[3]*mat.m21 + fTemp[4] *mat.m31;
        m00 -= fTemp[1]*mat.m11 + fTemp[2]*mat.m21 + fTemp[5] *mat.m31;
        m01  = fTemp[1]*mat.m01 + fTemp[6]*mat.m21 + fTemp[9] *mat.m31;
        m01 -= fTemp[0]*mat.m01 + fTemp[7]*mat.m21 + fTemp[8] *mat.m31;
        m02  = fTemp[2]*mat.m01 + fTemp[7]*mat.m11 + fTemp[10]*mat.m31;
        m02 -= fTemp[3]*mat.m01 + fTemp[6]*mat.m11 + fTemp[11]*mat.m31;
        m03  = fTemp[5]*mat.m01 + fTemp[8]*mat.m11 + fTemp[11]*mat.m21;
        m03 -= fTemp[4]*mat.m01 + fTemp[9]*mat.m11 + fTemp[10]*mat.m21;
        m10  = fTemp[1]*mat.m10 + fTemp[2]*mat.m20 + fTemp[5] *mat.m30;
        m10 -= fTemp[0]*mat.m10 + fTemp[3]*mat.m20 + fTemp[4] *mat.m30;
        m11  = fTemp[0]*mat.m00 + fTemp[7]*mat.m20 + fTemp[8] *mat.m30;
        m11 -= fTemp[1]*mat.m00 + fTemp[6]*mat.m20 + fTemp[9] *mat.m30;
        m12  = fTemp[3]*mat.m00 + fTemp[6]*mat.m10 + fTemp[11]*mat.m30;
        m12 -= fTemp[2]*mat.m00 + fTemp[7]*mat.m10 + fTemp[10]*mat.m30;
        m13  = fTemp[4]*mat.m00 + fTemp[9]*mat.m10 + fTemp[10]*mat.m20;
        m13 -= fTemp[5]*mat.m00 + fTemp[8]*mat.m10 + fTemp[11]*mat.m20;

        fTemp[0]  = mat.m20 * mat.m31;
        fTemp[1]  = mat.m30 * mat.m21;
        fTemp[2]  = mat.m10 * mat.m31;
        fTemp[3]  = mat.m30 * mat.m11;
        fTemp[4]  = mat.m10 * mat.m21;
        fTemp[5]  = mat.m20 * mat.m11;
        fTemp[6]  = mat.m00 * mat.m31;
        fTemp[7]  = mat.m30 * mat.m01;
        fTemp[8]  = mat.m00 * mat.m21;
        fTemp[9]  = mat.m20 * mat.m01;
        fTemp[10]  = mat.m00 * mat.m11;
        fTemp[11]  = mat.m10 * mat.m01;

        m20  = fTemp[0] *mat.m13 + fTemp[3] *mat.m23 + fTemp[4] *mat.m33;
        m20 -= fTemp[1] *mat.m13 + fTemp[2] *mat.m23 + fTemp[5] *mat.m33;
        m21  = fTemp[1] *mat.m03 + fTemp[6] *mat.m23 + fTemp[9] *mat.m33;
        m21 -= fTemp[0] *mat.m03 + fTemp[7] *mat.m23 + fTemp[8] *mat.m33;
        m22  = fTemp[2] *mat.m03 + fTemp[7] *mat.m13 + fTemp[10]*mat.m33;
        m22 -= fTemp[3] *mat.m03 + fTemp[6] *mat.m13 + fTemp[11]*mat.m33;
        m23  = fTemp[5] *mat.m03 + fTemp[8] *mat.m13 + fTemp[11]*mat.m23;
        m23 -= fTemp[4] *mat.m03 + fTemp[9] *mat.m13 + fTemp[10]*mat.m23;
        m30  = fTemp[2] *mat.m22 + fTemp[5] *mat.m32 + fTemp[1] *mat.m12;
        m30 -= fTemp[4] *mat.m32 + fTemp[0] *mat.m12 + fTemp[3] *mat.m22;
        m31  = fTemp[8] *mat.m32 + fTemp[0] *mat.m02 + fTemp[7] *mat.m22;
        m31 -= fTemp[6] *mat.m22 + fTemp[9] *mat.m32 + fTemp[1] *mat.m02;
        m32  = fTemp[6] *mat.m12 + fTemp[11]*mat.m32 + fTemp[3] *mat.m02;
        m32 -= fTemp[10]*mat.m32 + fTemp[2] *mat.m02 + fTemp[7] *mat.m12;
        m33  = fTemp[10]*mat.m22 + fTemp[4] *mat.m02 + fTemp[9] *mat.m12;
        m33 -= fTemp[8] *mat.m12 + fTemp[11]*mat.m22 + fTemp[5] *mat.m02;

        fDet = mat.m00 * m00 + mat.m10 * m01 + mat.m20 * m02 + mat.m30 * m03;

        if (IsZero(fDet)) {
            return false;
        }

        fDet = 1 / fDet;

        m00 *= fDet;
        m01 *= fDet;
        m02 *= fDet;
        m03 *= fDet;

        m10 *= fDet;
        m11 *= fDet;
        m12 *= fDet;
        m13 *= fDet;

        m20 *= fDet;
        m21 *= fDet;
        m22 *= fDet;
        m23 *= fDet;

        m30 *= fDet;
        m31 *= fDet;
        m32 *= fDet;
        m33 *= fDet;

        return true;
    }

    template<class T> inline bool Matrix44<T>::Invert(void) {
        Matrix44 ans;
        if (ans.InverseOf(*this)) {
            Set(ans);
            return true;
        }
        return false;
    }

    template<class T> inline bool Matrix44<T>::AffineInverseOf(const Matrix44<T>& mat) {
        // compute determinant
        T cofactor0 = mat.m11 * mat.m22 - mat.m21 * mat.m12;
        T cofactor3 = mat.m20 * mat.m12 - mat.m10 * mat.m22;
        T cofactor6 = mat.m10 * mat.m21 - mat.m20 * mat.m11;
        T det = mat.m00 * cofactor0 + mat.m01 * cofactor3 + mat.m02 * cofactor6;
        if (IsZero(det)) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        m00 = cofactor0 / det;
        m10 = cofactor3 / det;
        m20 = cofactor6 / det;
        
        m01 = (mat.m21 * mat.m02 - mat.m01 * mat.m22) / det;
        m11 = (mat.m00 * mat.m22 - mat.m20 * mat.m02) / det;
        m21 = (mat.m20 * mat.m01 - mat.m00 * mat.m21) / det;

        m02 = (mat.m01 * mat.m12 - mat.m11 * mat.m02) / det;
        m12 = (mat.m10 * mat.m02 - mat.m00 * mat.m12) / det;
        m22 = (mat.m00 * mat.m11 - mat.m10 * mat.m01) / det;

        m03 = -m00 * mat.m03 - m01 * mat.m13 - m02 * mat.m23;
        m13 = -m10 * mat.m03 - m11 * mat.m13 - m12 * mat.m23;
        m23 = -m20 * mat.m03 - m21 * mat.m13 - m22 * mat.m23;

        return true;
    }

    template<> inline bool Matrix44<float32>::AffineInverseOf(const Matrix44<float32>& mat) {
        // compute determinant
        float32 cofactor0 = mat.m11 * mat.m22 - mat.m21 * mat.m12;
        float32 cofactor3 = mat.m20 * mat.m12 - mat.m10 * mat.m22;
        float32 cofactor6 = mat.m10 * mat.m21 - mat.m20 * mat.m11;
        float32 det = mat.m00 * cofactor0 + mat.m01 * cofactor3 + mat.m02 * cofactor6;
        if (IsZero(det)) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        det = 1.0f / det;
        m00 = det * cofactor0;
        m10 = det * cofactor3;
        m20 = det * cofactor6;
        
        m01 = det * (mat.m21 * mat.m02 - mat.m01 * mat.m22);
        m11 = det * (mat.m00 * mat.m22 - mat.m20 * mat.m02);
        m21 = det * (mat.m20 * mat.m01 - mat.m00 * mat.m21);

        m02 = det * (mat.m01 * mat.m12 - mat.m11 * mat.m02);
        m12 = det * (mat.m10 * mat.m02 - mat.m00 * mat.m12);
        m22 = det * (mat.m00 * mat.m11 - mat.m10 * mat.m01);

        m03 = -m00 * mat.m03 - m01 * mat.m13 - m02 * mat.m23;
        m13 = -m10 * mat.m03 - m11 * mat.m13 - m12 * mat.m23;
        m23 = -m20 * mat.m03 - m21 * mat.m13 - m22 * mat.m23;

        return true;
    }

    template<> inline bool Matrix44<float64>::AffineInverseOf(const Matrix44<float64>& mat) {
        // compute determinant
        float64 cofactor0 = mat.m11 * mat.m22 - mat.m21 * mat.m12;
        float64 cofactor3 = mat.m20 * mat.m12 - mat.m10 * mat.m22;
        float64 cofactor6 = mat.m10 * mat.m21 - mat.m20 * mat.m11;
        float64 det = mat.m00 * cofactor0 + mat.m01 * cofactor3 + mat.m02 * cofactor6;
        if (IsZero(det)) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        det = 1.0 / det;
        m00 = det * cofactor0;
        m10 = det * cofactor3;
        m20 = det * cofactor6;
        
        m01 = det * (mat.m21 * mat.m02 - mat.m01 * mat.m22);
        m11 = det * (mat.m00 * mat.m22 - mat.m20 * mat.m02);
        m21 = det * (mat.m20 * mat.m01 - mat.m00 * mat.m21);

        m02 = det * (mat.m01 * mat.m12 - mat.m11 * mat.m02);
        m12 = det * (mat.m10 * mat.m02 - mat.m00 * mat.m12);
        m22 = det * (mat.m00 * mat.m11 - mat.m10 * mat.m01);

        m03 = -m00 * mat.m03 - m01 * mat.m13 - m02 * mat.m23;
        m13 = -m10 * mat.m03 - m11 * mat.m13 - m12 * mat.m23;
        m23 = -m20 * mat.m03 - m21 * mat.m13 - m22 * mat.m23;

        return true;
    }

    template<class T> inline bool Matrix44<T>::AffineInvert(void) {
        Matrix44<T> ans;
        if (ans.AffineInverseOf(*this)) {
            Set(ans);
            return true;
        }
        return false;
    }

    template<class T> inline void Matrix44<T>::Transpose(void) {
        std::swap(m01, m10);
        std::swap(m02, m20);
        std::swap(m03, m30);
        std::swap(m12, m21);
        std::swap(m13, m31);
        std::swap(m23, m32);
    }

    template<class T> inline void Matrix44<T>::TransposeOf(const Matrix44<T> &mat) {
        m00 = mat.m00;
        m01 = mat.m10;
        m02 = mat.m20;
        m03 = mat.m30;

        m10 = mat.m01;
        m11 = mat.m11;
        m12 = mat.m21;
        m13 = mat.m31;

        m20 = mat.m02;
        m21 = mat.m12;
        m22 = mat.m22;
        m23 = mat.m32;

        m30 = mat.m03;
        m31 = mat.m13;
        m32 = mat.m23;
        m33 = mat.m33;
    }

    template<class T> inline void Matrix44<T>::PerspectiveGL(T fov, T aspect, T n, T f) {
        Zero();

        if (IsEqual(n, f))
            return;
        
        T d = 1 / tanf(fov * 0.5f);
        T invNearMinFar = 1 / (n - f);

        m00 = d / aspect;
        m11 = d;
        m22 = (n + f) * invNearMinFar;
        m23 = 2 * n * f * invNearMinFar;
        m32 = -1.0f;
    }

    template<class T> inline void Matrix44<T>::PerspectiveD3D(T fov, T aspect, T n, T f) {
        Zero();

        if (IsEqual(n, f))
            return;
        
        T d = 1 / tanf(fov * 0.5f);
        T invFarMinNear = 1 / (f - n);

        m00 = d / aspect;
        m11 = d;
        m22 = f * invFarMinNear;
        m23 = -n * f * invFarMinNear;
        m32 = 1.0f;
    }

    template<class T> inline void Matrix44<T>::ObliquePerspGL(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        T inv = 1 / (right - left);
        m00 = 2 * n * inv;
        m02 = (right + left) * inv;
        
        inv = 1 / (top - bottom);
        m11 = 2 * n * inv;
        m12 = (top + bottom) * inv;

        inv = 1 / (n - f);
        m22 = (n + f) * inv;
        m23 = 2 * n * f * inv;

        m32 = -1;
    }

    template<class T> inline void Matrix44<T>::ObliquePerspD3D(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        float32 inv = 1 / (right - left);
        m00 = 2 * n * inv;
        m02 = -(right + left) * inv;
        
        inv = 1 / (top - bottom);
        m11 = 2 * n * inv;
        m12 = -(top + bottom) * inv;

        inv = 1 / (f - n);
        m22 = f * inv;
        m23 = -n * f * inv;

        m32 = 1;
    }

    template<class T> inline void Matrix44<T>::OrthographicGL(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        T inv = 1 / (right - left);
        m00 = 2 * inv;
        m03 = -(right + left) * inv;

        inv = 1 / (top - bottom);
        m11 = 2 * inv;
        m13 = -(top + bottom) * inv;

        inv = -1 / (f - n);
        m22 = 2 * inv;
        m23 = (f + n) * inv;

        m33 = 1.0f;
    }

    template<class T> inline void Matrix44<T>::OrthographicD3D(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        T inv = 1 / (right - left);
        m00 = 2 * inv;
        m03 = -(right + left) * inv;

        inv = 1 / (top - bottom);
        m11 = 2 * inv;
        m13 = -(top + bottom) * inv;

        m22 = 1 / (f - n);
        m23 = -(n) * m22;

        m33 = 1.0f;
    }

    template<class T> inline void Matrix44<T>::ObliqueParallelGL(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        m02 = 1 / (right - left);;
        m00 = 2 * m02;
        m03 = -(right + left - n) * m02;

        m12 = 1 / (top - bottom);
        m11 = 2 * m12;
        m13 = -(top + bottom - n) * m12;

        T inv = -1 / (f - n);
        m22 = 2 * inv;
        m23 = (f + n) * inv;

        m33 = 1.0f;
    }

    template<class T> inline void Matrix44<T>::ObliqueParallelD3D(T left, T right, T bottom, T top, T n, T f) {
        Zero();
        if (IsEqual(top, bottom) || IsEqual(right, left) || IsEqual(n, f))
            return;

        m02 = 1 / (right - left);;
        m00 = 2 * m02;
        m03 = -(right + left - n) * m02;

        m12 = 1 / (top - bottom);
        m11 = 2 * m12;
        m13 = -(top + bottom - n) * m12;

        m22 = 1 / (f - n);
        m23 = -n * m22;

        m33 = 1.0f;
    }

    template<class T> inline T Matrix44<T>::Determinant(void) const {
        return m00 * DetPart(m11, m21, m31, m12, m22, m32, m13, m23, m33) -
                 m01 * DetPart(m10, m20, m30, m12, m22, m32, m13, m23, m33) +
                 m02 * DetPart(m10, m20, m30, m11, m21, m31, m13, m23, m33) -
                 m03 * DetPart(m10, m20, m30, m11, m21, m31, m12, m22, m32);
    }

    template<class T> inline T Matrix44<T>::Trace(void) const {
        return m00 + m11 + m22 + m33;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::Tensor(const Vector4<T>& v1, const Vector4<T>& v2) {
        m00 = v1.x * v2.x;
        m10 = v1.y * v2.x;
        m20 = v1.z * v2.x;
        m30 = v1.w * v2.x;

        m01 = v1.x * v2.y;
        m11 = v1.y * v2.y;
        m21 = v1.z * v2.y;
        m31 = v1.w * v2.y;

        m02 = v1.x * v2.z;
        m12 = v1.y * v2.z;
        m22 = v1.z * v2.z;
        m32 = v1.w * v2.z;

        m03 = v1.x * v2.w;
        m13 = v1.y * v2.w;
        m23 = v1.z * v2.w;
        m33 = v1.w * v2.w;

        return *this;
    }

    template<class T> inline const Matrix44<T> Matrix44<T>::operator * (const Matrix44<T> &rhs) const {
        Matrix44d mResult;
        mResult.m00 = m00 * rhs.m00 + m01 * rhs.m10 + m02 * rhs.m20 + m03 * rhs.m30;
        mResult.m01 = m00 * rhs.m01 + m01 * rhs.m11 + m02 * rhs.m21 + m03 * rhs.m31;
        mResult.m02 = m00 * rhs.m02 + m01 * rhs.m12 + m02 * rhs.m22 + m03 * rhs.m32;
        mResult.m03 = m00 * rhs.m03 + m01 * rhs.m13 + m02 * rhs.m23 + m03 * rhs.m33;

        mResult.m10 = m10 * rhs.m00 + m11 * rhs.m10 + m12 * rhs.m20 + m13 * rhs.m30;
        mResult.m11 = m10 * rhs.m01 + m11 * rhs.m11 + m12 * rhs.m21 + m13 * rhs.m31;
        mResult.m12 = m10 * rhs.m02 + m11 * rhs.m12 + m12 * rhs.m22 + m13 * rhs.m32;
        mResult.m13 = m10 * rhs.m03 + m11 * rhs.m13 + m12 * rhs.m23 + m13 * rhs.m33;

        mResult.m20 = m20 * rhs.m00 + m21 * rhs.m10 + m22 * rhs.m20 + m23 * rhs.m30;
        mResult.m21 = m20 * rhs.m01 + m21 * rhs.m11 + m22 * rhs.m21 + m23 * rhs.m31;
        mResult.m22 = m20 * rhs.m02 + m21 * rhs.m12 + m22 * rhs.m22 + m23 * rhs.m32;
        mResult.m23 = m20 * rhs.m03 + m21 * rhs.m13 + m22 * rhs.m23 + m23 * rhs.m33;

        mResult.m30 = m30 * rhs.m00 + m31 * rhs.m10 + m32 * rhs.m20 + m33 * rhs.m30;
        mResult.m31 = m30 * rhs.m01 + m31 * rhs.m11 + m32 * rhs.m21 + m33 * rhs.m31;
        mResult.m32 = m30 * rhs.m02 + m31 * rhs.m12 + m32 * rhs.m22 + m33 * rhs.m32;
        mResult.m33 = m30 * rhs.m03 + m31 * rhs.m13 + m32 * rhs.m23 + m33 * rhs.m33;
        
        return mResult;
    }

    template<> inline const Matrix44<float32> Matrix44<float32>::operator *(const Matrix44<float32> &rhs) const {
        Matrix44<float32> ans;
        float32* mOut = ans.pntr;
        const float32* mlhs = pntr;
        const float32* mrhs = rhs.pntr;

        __asm 
        {
            mov      edx,         mrhs            ; src1
            mov      eax,         mOut          ; dst
            mov      ecx,         mlhs          ; src2
            movss    xmm0,        [edx]
            movups  xmm1,        [ecx]
            shufps  xmm0,        xmm0,        0
            movss    xmm2,        [edx+4]
            mulps    xmm0,        xmm1
            shufps  xmm2,        xmm2,        0
            movups  xmm3,        [ecx+10h]
            movss    xmm7,        [edx+8]
            mulps    xmm2,        xmm3
            shufps  xmm7,        xmm7,        0
            addps    xmm0,        xmm2
            movups  xmm4,        [ecx+20h]
            movss    xmm2,        [edx+0Ch]
            mulps    xmm7,        xmm4
            shufps  xmm2,        xmm2,        0
            addps    xmm0,        xmm7
            movups  xmm5,        [ecx+30h]
            movss    xmm6,        [edx+10h]
            mulps    xmm2,        xmm5
            movss    xmm7,        [edx+14h]
            shufps  xmm6,        xmm6,        0
            addps    xmm0,        xmm2
            shufps  xmm7,        xmm7,        0
            movlps  [eax],      xmm0
            movhps  [eax+8],    xmm0
            mulps    xmm7,        xmm3
            movss    xmm0,        [edx+18h]
            mulps    xmm6,        xmm1
            shufps  xmm0,        xmm0,        0
            addps    xmm6,        xmm7
            mulps    xmm0,        xmm4
            movss    xmm2,        [edx+24h]
            addps    xmm6,        xmm0
            movss    xmm0,        [edx+1Ch]
            movss    xmm7,        [edx+20h]
            shufps  xmm0,        xmm0,        0
            shufps  xmm7,        xmm7,        0
            mulps    xmm0,        xmm5
            mulps    xmm7,        xmm1
            addps    xmm6,        xmm0
            shufps  xmm2,        xmm2,        0
            movlps  [eax+10h], xmm6
            movhps  [eax+18h], xmm6
            mulps    xmm2,        xmm3
            movss    xmm6,        [edx+28h]
            addps    xmm7,        xmm2
            shufps  xmm6,        xmm6,        0
            movss    xmm2,        [edx+2Ch]
            mulps    xmm6,        xmm4
            shufps  xmm2,        xmm2,        0
            addps    xmm7,        xmm6
            mulps    xmm2,        xmm5
            movss    xmm0,        [edx+34h]
            addps    xmm7,        xmm2
            shufps  xmm0,        xmm0,        0
            movlps  [eax+20h], xmm7
            movss    xmm2,        [edx+30h]
            movhps  [eax+28h], xmm7
            mulps    xmm0,        xmm3
            shufps  xmm2,        xmm2,        0
            movss    xmm6,        [edx+38h]
            mulps    xmm2,        xmm1
            shufps  xmm6,        xmm6,        0
            addps    xmm2,        xmm0
            mulps    xmm6,        xmm4
            movss    xmm7,        [edx+3Ch]
            shufps  xmm7,        xmm7,        0
            addps    xmm2,        xmm6
            mulps    xmm7,        xmm5
            addps    xmm2,        xmm7
            movups  [eax+30h], xmm2
        } // asm
        return ans;
    }

    template<class T> inline const Matrix44<T> Matrix44<T>::operator *(      T    rhs) const {
        return(Matrix44<T>(*this) *= rhs);
    }

    template<class T> inline const Matrix44<T> Matrix44<T>::operator /(      T    rhs) const {
        return(Matrix44<T>(*this) /= rhs);
    }

    template<class T> inline const Matrix44<T> Matrix44<T>::operator +(const Matrix44<T> &rhs) const {
        return(Matrix44<T>(*this) += rhs);
    }

    template<class T> inline const Matrix44<T> Matrix44<T>::operator -(const Matrix44<T> &rhs) const {
        return(Matrix44<T>(*this) -= rhs);
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::operator *=(const Matrix44<T> &rhs) {
        return ((*this) = (*this) * rhs);
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::operator *=(      T    rhs) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] *= rhs;
        }
        return *this;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::operator /=(      T    rhs) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] /= rhs;
        }
        return *this;
    }

    template<> inline Matrix44<float32>& Matrix44<float32>::operator /=(float32 rhs) {
        rhs = 1.0f / rhs;
        return (*this) *= rhs;
    }

    template<> inline Matrix44<float64>& Matrix44<float64>::operator /=(float64 rhs) {
        rhs = 1.0 / rhs;
        return (*this) *= rhs;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::operator +=(const Matrix44<T> &rhs) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] += rhs.pntr[ind];
        }
        return *this;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::operator -=(const Matrix44<T> &rhs) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] -= rhs.pntr[ind];
        }
        return *this;
    }

    // Methods used to replace the = sign. Use as ans = m1 op m2
    template<class T> inline Matrix44<T>& Matrix44<T>::Add(const Matrix44<T>& m1, const Matrix44<T>& m2) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] = m1.pntr[ind] + m2.pntr[ind];
        }
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::Sub(const Matrix44<T>& m1, const Matrix44<T>& m2) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] = m1.pntr[ind] - m2.pntr[ind];
        }
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::Mul(const Matrix44<T>& m1, const Matrix44<T>& m2) {
        float32* mOut = pntr;
        const float32* mlhs = m1.pntr;
        const float32* mrhs = m2.pntr;
        __asm 
        {
            mov      edx,         mrhs             ; src1
            mov      eax,         mOut          ; dst
            mov      ecx,         mlhs             ; src2
            movss    xmm0,        [edx]
            movups  xmm1,        [ecx]
            shufps  xmm0,        xmm0,        0
            movss    xmm2,        [edx+4]
            mulps    xmm0,        xmm1
            shufps  xmm2,        xmm2,        0
            movups  xmm3,        [ecx+10h]
            movss    xmm7,        [edx+8]
            mulps    xmm2,        xmm3
            shufps  xmm7,        xmm7,        0
            addps    xmm0,        xmm2
            movups  xmm4,        [ecx+20h]
            movss    xmm2,        [edx+0Ch]
            mulps    xmm7,        xmm4
            shufps  xmm2,        xmm2,        0
            addps    xmm0,        xmm7
            movups  xmm5,        [ecx+30h]
            movss    xmm6,        [edx+10h]
            mulps    xmm2,        xmm5
            movss    xmm7,        [edx+14h]
            shufps  xmm6,        xmm6,        0
            addps    xmm0,        xmm2
            shufps  xmm7,        xmm7,        0
            movlps  [eax],      xmm0
            movhps  [eax+8],    xmm0
            mulps    xmm7,        xmm3
            movss    xmm0,        [edx+18h]
            mulps    xmm6,        xmm1
            shufps  xmm0,        xmm0,        0
            addps    xmm6,        xmm7
            mulps    xmm0,        xmm4
            movss    xmm2,        [edx+24h]
            addps    xmm6,        xmm0
            movss    xmm0,        [edx+1Ch]
            movss    xmm7,        [edx+20h]
            shufps  xmm0,        xmm0,        0
            shufps  xmm7,        xmm7,        0
            mulps    xmm0,        xmm5
            mulps    xmm7,        xmm1
            addps    xmm6,        xmm0
            shufps  xmm2,        xmm2,        0
            movlps  [eax+10h], xmm6
            movhps  [eax+18h], xmm6
            mulps    xmm2,        xmm3
            movss    xmm6,        [edx+28h]
            addps    xmm7,        xmm2
            shufps  xmm6,        xmm6,        0
            movss    xmm2,        [edx+2Ch]
            mulps    xmm6,        xmm4
            shufps  xmm2,        xmm2,        0
            addps    xmm7,        xmm6
            mulps    xmm2,        xmm5
            movss    xmm0,        [edx+34h]
            addps    xmm7,        xmm2
            shufps  xmm0,        xmm0,        0
            movlps  [eax+20h], xmm7
            movss    xmm2,        [edx+30h]
            movhps  [eax+28h], xmm7
            mulps    xmm0,        xmm3
            shufps  xmm2,        xmm2,        0
            movss    xmm6,        [edx+38h]
            mulps    xmm2,        xmm1
            shufps  xmm6,        xmm6,        0
            addps    xmm2,        xmm0
            mulps    xmm6,        xmm4
            movss    xmm7,        [edx+3Ch]
            shufps  xmm7,        xmm7,        0
            addps    xmm2,        xmm6
            mulps    xmm7,        xmm5
            addps    xmm2,        xmm7
            movups  [eax+30h], xmm2
        } // asm
        return *this;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::Mul(const Matrix44<T>& m1,         T    f) {
        for(int32 ind = 0; ind <16; ++ind) {
            pntr[ind] = m1.pntr[ind] * f;
        }

        return *this;
    }

    template<class T> inline Matrix44<T>& Matrix44<T>::Div(const Matrix44<T>& m1,         T    f) {
        return (*this = m1) /= f;
    }

    template<class T> inline const Vector4<T> Matrix44<T>::operator * (const Vector4<T>& rhs) const {
        Vector4<T> vResult;
        vResult.x = m00 * rhs.x + m01 * rhs.y + m02 * rhs.z + m03 * rhs.w;
        vResult.y = m10 * rhs.x + m11 * rhs.y + m12 * rhs.z + m13 * rhs.w;
        vResult.z = m20 * rhs.x + m21 * rhs.y + m22 * rhs.z + m23 * rhs.w;
        vResult.w = m30 * rhs.x + m31 * rhs.y + m32 * rhs.z + m33 * rhs.w;
        
        return vResult;
    }

    template<> inline const Vector4<float32> Matrix44<float32>::operator * (const Vector4<float32>& rhs) const {
        Vector4<float32> ans;
                float32* vout = ans.pntr;
        const float32* vin = rhs.pntr;
        const float32* min = pntr;
        __asm {
            mov      ecx,     vin
            mov      edx,     min
            mov      eax,     vout
            movss    xmm0,    [ecx]
            movss    xmm1,    [ecx+4]
            shufps  xmm0,    xmm0,        0
            movss    xmm2,    [ecx+8]
            movups  xmm4,    [edx]
            mulps    xmm0,    xmm4
            shufps  xmm1,    xmm1,        0
            movss    xmm3,    [ecx+12]
            movups  xmm4,    [edx+16]
            mulps    xmm1,    xmm4
            shufps  xmm2,    xmm2,        0
            shufps  xmm3,    xmm3,        0
            movups  xmm4,    [edx+32]
            mulps    xmm2,    xmm4
            addps    xmm0,    xmm1
            movups  xmm4,    [edx+48]
            mulps    xmm3,    xmm4
            addps    xmm2,    xmm3
            addps    xmm0,    xmm2
            
            movups    [eax],    xmm0
        };
        return ans;
    }

    template<class T> inline bool Matrix44<T>::operator==(const Matrix44<T>& other) {
        return (this == &other) || (m00 == other.m00 && m10 == other.m10 && m20 == other.m20 && m30 == other.m30 &&
                                              m01 == other.m01 && m11 == other.m11 && m21 == other.m21 && m31 == other.m31 &&
                                              m02 == other.m02 && m12 == other.m12 && m22 == other.m22 && m32 == other.m32 &&
                                              m03 == other.m03 && m13 == other.m13 && m23 == other.m23 && m33 == other.m33);
    }

    template<class T> inline bool Matrix44<T>::operator!=(const Matrix44<T>& other) {
        return !((*this) == other);
    }

    template<class T> inline bool Matrix44<T>::GetSymEigenInfo(Matrix44<T>& ev, Vector4<T>& ew) const {
        ev = *this;
        return ev.SymEigenDecompose(ew);
    }

    template<class T> inline bool Matrix44<T>::SymEigenDecompose(Vector4<T>& ew) {
        throw "Can't use method Matrix44<T>::SymEigenDecompose with the requested type.";
        return false;
    }

    template<> inline bool Matrix44<float32>::SymEigenDecompose(Vector4<float32>& ew) {
        Vector4<float32> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<> inline bool Matrix44<float64>::SymEigenDecompose(Vector4<float64>& ew) {
        Vector4<float64> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<class T> inline void Matrix44<T>::Tridiagonalize(Vector4<T>& diag, Vector4<T>& subd) {
        int32 i0, i1, i2, i3;

        T h, scale;
        for(i0 = 4 - 1, i3 = 4 - 2; i0>= 2; --i0, --i3)
        {
            h = 0.0f;
            scale = 0.0f;
            for(i2 = 0; i2 <= i3; ++i2)
                scale += fabs(pntr[i0 + (i2 <<2)]);

            if (scale == 0)
                subd[i0] = pntr[i0 + (i3 <<2)];
            else
            {
                T invScale = 1.0f / scale;
                for (i2 = 0; i2 <= i3; ++i2)
                {
                    pntr[i0 + (i2 <<2)] *= invScale;
                    h += pntr[i0 + (i2 <<2)] * pntr[i0 + (i2 <<2)];
                }

                T f = pntr[i0 + (i3 <<2)];
                T g = sqrt(h);
                if (f> 0.0f)
                    g = -g;

                subd[i0] = scale * g;
                h -= f * g;
                pntr[i0 + (i3 <<2)] = f - g;
                f = 0.0f;
                T invH = 1.0f / h;
                for (i1 = 0; i1 <= i3; ++i1)
                {
                    pntr[i1 + (i0 <<2)] = pntr[i0 + (i1 <<2)] * invH;
                    g = 0.0f;
                    for (i2 = 0; i2 <= i1; ++i2)
                        g += pntr[i1 + (i2 <<2)] * pntr[i0 + (i2 <<2)];

                    for (i2 = i1+1; i2 <= i3; ++i2)
                        g += pntr[i2 + (i1 <<2)] * pntr[i0 + (i2 <<2)];

                    subd[i1] = g * invH;
                    f += subd[i1] * pntr[i0 + (i1 <<2)];
                }

                T halfdivH = 0.5f * f * invH;
                for (i1 = 0; i1 <= i3; i1++)
                {
                    f = pntr[i0 + (i1 <<2)];
                    g = subd[i1] - halfdivH*f;
                    subd[i1] = g;
                    for (i2 = 0; i2 <= i1; ++i2)
                        pntr[i1 + (i2 <<2)] -= f * subd[i2] + g * pntr[i0 + (i2 <<2)];
                }
            }

            diag[i0] = h;
        }

        subd[1] = pntr[1];
        diag[1] = 0.0f;

        diag[0] = subd[0] = 0;
        for (i0 = 0, i3 = 0; i0 <4; ++i0, ++i3)
        {
            if (diag[i0] != 0)
            {
                for (i1 = 0; i1 <i3; ++i1)
                {
                    T sum = 0.0f;
                    for (i2 = 0; i2 <i3; ++i2)
                        sum += pntr[i0 + (i2 <<2)] * pntr[i2 + (i1 <<2)];

                    for (i2 = 0; i2 <i3; ++i2)
                        pntr[i2 + (i1 <<2)] -= sum * pntr[i2 + (i0 <<2)];
                }
            }

            diag[i0] = pntr[i0 * 5];
            pntr[i0 * 5] = 1;
            for (i1 = 0; i1 <i3; ++i1)
                pntr[i1 + (i0 <<2)] = pntr[i0 + (i1 <<2)] = 0.0f;
        }

        // re-ordering if MgcEigen::QLAlgorithm is used subsequently
        for (i0 = 1, i3 = 0; i0 <4; i0++, i3++)
            subd[i3] = subd[i0];

        subd[3] = 0;
    }

    template<class T> inline bool Matrix44<T>::QLFactorize(Vector4<T>& diag, Vector4<T>& subd) {
	    const int32 iMaxIter = 32;

		int32 i0, i1, i2, i3, i4;
		for (i0 = 0; i0 <4; i0++)
		{
			for (i1 = 0; i1 <iMaxIter; i1++)
			{
				i2 = i0;
				while(i2 <3 && IsNotZero(subd[i2]))
					++i2;

				if (i2 == i0)
					break;

				T g = (diag[i0 + 1] - diag[i0]) / (2.0f * subd[i0]);
				T r = sqrt(g * g + 1.0f);

				if (g <0.0f)
					g = diag[i2] - diag[i0] + subd[i0] / (g - r);
				else
					g = diag[i2] - diag[i0] + subd[i0] / (g + r);
				
				T fSin = 1.0f, fCos = 1.0f, p = 0.0f;
				for (i3 = i2; i3> i0; --i3)
				{
					T f = fSin * subd[i3 - 1];
					T b = fCos * subd[i3 - 1];
					if (fabs(f)>= fabs(g))
					{
						fCos = g / f;
						r = sqrt(fCos * fCos + 1.0f);
						subd[i3] = f * r;
						fSin = 1.0f / r;
						fCos *= fSin;
					}
					else
					{
						fSin = f / g;
						r = sqrt(fSin * fSin + 1.0f);
						subd[i3] = g * r;
						fCos = 1.0f / r;
						fSin *= fCos;
					}
					g = diag[i3] - p;
					r = (diag[i3 - 1] - g) * fSin + 2.0f * b * fCos;
					p = fSin * r;
					diag[i3] = g + p;
					g = fCos * r - b;

					for (i4 = 0; i4 <4; i4++)
					{
						f = pntr[i4 + (i3 <<2)];
						pntr[i4 +    (i3         <<2)] = fSin * pntr[i4 + ((i3 - 1) <<2)] + fCos * f;
						pntr[i4 + ((i3 - 1) <<2)] = fCos * pntr[i4 + ((i3 - 1) <<2)] - fSin * f;
					}
				}
				diag[i0] -= p;
				subd[i0] = g;
				subd[i2] = 0.0f;
			}
			if (i1 == iMaxIter)
				return false;
		}
		return true;
    }

    template<class T> inline void Matrix44<T>::SortEigenInfo(Vector4<T>& diag) {
        int32 i0, i1, i2;
		for(i0 = 0; i0 <3; ++i0) {
			// locate maximum eigenvalue
			i1 = i0;
			T max = fabs(diag[i1]);
			for(i2 = i0 + 1; i2 <4; ++i2) {
				if (fabs(diag[i2])> fabs(max)) {
					i1 = i2;
					max = diag[i1];
				}
			}

			if (i1 != i0) {
				// swap eigenvalues
				diag[i1] = diag[i0];
				diag[i0] = max;

				// swap eigenvectors
                for(i2 = 0; i2 <4; ++i2) {
                    std::swap(pntr[i2 + (i0 <<2)], pntr[i2 + (i1 <<2)]);
                }
			}
		}
    }


    template<class T> inline T Matrix44<T>::DetPart(T c0, T c1, T c2, T c3, T c4, T c5, T c6, T c7, T c8) const {
        return  (c0 * (c4 * c8 - c5 * c7) +
                     c1 * (c5 * c6 - c3 * c8) +
                     c2 * (c3 * c7 - c4 * c6));
    }

}