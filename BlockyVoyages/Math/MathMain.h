#pragma once

#include "../types.h"

#include <iostream>
#include <math.h>

// put all types in here then link to *.hpp

namespace BlockyVoyages {
    // floating point types
    template<class T = float32> class Vector2;
    template<class T = float32> class Vector3;
    template<class T = float32> class Vector4;
    template<class T = float32> class VectorN;

    template<class T = float32> class Matrix22;
    template<class T = float32> class Matrix33;
    template<class T = float32> class Matrix44;
    template<class T = float32> class MatrixMN;

    template<class T = float32> class Quaternion;

    // global constants
    extern const float32 PI;          // Pi
    extern const float32 PIDIV2;     // Pi / 2
    extern const float32 PI2;         // Pi * 2
    extern const float32 PIDIV4;     // Pi / 4
    extern const float32 PIDIV16;    // Pi / 16
    extern const float32 PIINV;      // 1 / Pi
    extern const float32 G_FT;        // Acceleration of an object due to gravity in ft / s^2
    extern const float32 G_M;         // Acceleration of an object due to gravity in m / s^2
    extern const float32 EPSILON;    // some small constant 1e-10 in this instance
    extern const float32 TINY;        // an extremely constant 1e-20 in this instance

    // the actual class declarations
    // floating point types
    template<class T> class Vector2 {
    public:
        // data
        union {
            struct {
                T x, y;
            };

            struct {
                T s, t;
            };

            struct {
                T u, v;
            };

            T pntr[2];
        };

        // methods
        // constructors
        Vector2(void);

        Vector2(const Vector2<T>& other);
        explicit Vector2(const Vector3<T>& other);
        explicit Vector2(const Vector4<T>& other);

        // set from floats
        Vector2(T _x, T _y);
        Vector2(const T _v[2]);

        ~Vector2(void);

        // setting methods
        const Vector2<T>& operator=(const Vector2<T>& other);

        void Set(const Vector2<T>& other);
        void Set(const Vector3<T>& other);
        void Set(const Vector4<T>& other);

        // set from floats
        void Set(T _x, T _y);
        void Set(const T _v[2]);

        // clean up and fast implementation
        void Zero(void);
        void Clean(void);

        Vector2<T>& Add  (const Vector2<T>&  v1,    const Vector2<T>& v2     );
        Vector2<T>& Sub  (const Vector2<T>&  v1,    const Vector2<T>& v2     );
        Vector2<T>& Mul  (const Vector2<T>&  v1,    T scalar                     );
        Vector2<T>& Mul  (const Vector2<T>&  v1,    const Matrix22<T>& mat);
        Vector2<T>& Mul  (const Matrix22<T>& mat,  const Vector2<T>& v1  );
        Vector2<T>& Div  (const Vector2<T>&  v1,    T scalar                     );
        Vector2<T>& Build(const Vector2<T>&  init, const Vector2<T>& term);

        T Dot(const Vector2<T>& v2) const;                 // computes the dot product between this and v2

        const Matrix22<T> Tensor(const Vector2<T>& v2) const;

        const Vector2<T> operator + (const Vector2<T>& rhs) const;        // Adds two vectors. use Add() or += instead if possible
        const Vector2<T> operator - (const Vector2<T>& rhs) const;        // subtracts two vectors. use Sub() or -= instead if possible
                T                operator * (const Vector2<T>& rhs) const;        // computes the dot product between this and rhs
        const Vector2<T> operator * (T scalar) const;                            // multiplies this by a scalar. use Mul() or *= instead if possible
        const Vector2<T> operator * (const Matrix22<T>& rhs) const;
        const Vector2<T> operator / (T scalar) const;                            // divides this by a scalar. use Div() or /= instead if possible

        Vector2<T>& operator += (const Vector2<T>& rhs);     // does an in-place addition with rhs
        Vector2<T>& operator -= (const Vector2<T>& rhs);     // does an in-place subtraction with rhs
        Vector2<T>& operator *= (T scalar);          // does an in-place scalar multiplication
        Vector2<T>& operator *= (const Matrix22<T>& mat);          // does an in-place matrix multiplication
        Vector2<T>& operator /= (T scalar);          // does an in-place scalar division

        T LengthSquared(void) const;
        T Length(void) const;

        Vector2<T>& NegationOf(const Vector2<T>& other);
        Vector2<T>& Negate(void);

        Vector2& UnitOf(const Vector2<T>& other);
        Vector2& Normalize(void);

        T CosAngleBetween(const Vector2<T>& vect) const; // computes the cosine of the angle between this vector and vect
        T AngleBetween(const Vector2<T>& vect) const;     // computes the angle between this vector and vect

        T DistanceSquared(const Vector2<T>& v2) const; // finds the square of the distance between this and v2
        T Distance(const Vector2<T>& v2) const;          // finds the distance between this and v2

        Vector2<T>& PerpendicularOf(const Vector2<T>& v);            // sets this vector to the 2D perpendicular vector to v
        Vector2<T> Perpendicular(void) const;                     // returns the vector perpendicular to this

        bool operator == (const Vector2<T>& rhs) const;          // compares this to rhs and returns if they are equal
        bool operator != (const Vector2<T>& rhs) const;          // compares this to rhs and returns if they are not equal

        // Accessors
        T& operator [] (int32 i)                { return pntr[i]; }
        T  operator [] (int32 i) const        { return pntr[i]; }

        friend std::ostream& operator<<(std::ostream& out, const Vector2<T>& source) {     // outputs the vector to the streams
            out <<"<" <<source.x <<", " <<source.y <<"> ";
            return out;
        }

        friend const Vector2<T> operator * (T scalar, const Vector2<T>& rhs) {        // returns a scalar multiplied by this vector 
            return Vector2<T>(rhs.x * scalar, rhs.y * scalar);
        }

        friend const Vector2<T> operator - (const Vector2<T>& v) {                                // returns the negation of this vector
            return Vector2<T>(-v.x, -v.y);
        }
    };

    template<class T> class Vector3 {
    public:
        // data
        union {
            struct {
                T x, y, z;
            };

            struct {
                T r, s, t;
            };

            struct {
                T r, g, b;
            };

            struct {
                T u, v, w;
            };

            T pntr[3];
        };

        // methods
        // constructors
        Vector3(void);

        explicit Vector3(const Vector2<T>& other, T _z);
        Vector3(const Vector3<T>& other);
        explicit Vector3(const Vector4<T>& other);

        // set from floats
        Vector3(T _x, T _y, T _z);
        Vector3(const T _v[3]);

        ~Vector3(void);

        // setting methods
        Vector3<T>& operator=(const Vector3<T>& other);

        void Set(const Vector2<T>& other, T _z);
        void Set(const Vector3<T>& other);
        void Set(const Vector4<T>& other);

        // set from floats
        void Set(T _x, T _y, T _z);
        void Set(const T _v[3]);

        // clean up and fast implementation
        void Zero(void);
        void Clean(void);

        Vector3<T>& Add(const Vector3<T>&  lhs, const Vector3<T>& rhs);
        Vector3<T>& Sub(const Vector3<T>&  lhs, const Vector3<T>& rhs);
        Vector3<T>& Mul(const Vector3<T>&  lhs, T rhs);
        Vector3<T>& Mul(const Vector3<T>&  lhs, const Matrix33<T>& rhs);
        Vector3<T>& Mul(const Matrix33<T>& lhs, const Vector3<T>& rhs);
        Vector3<T>& Div(const Vector3<T>&  lhs, T rhs);
        Vector3<T>& Build(const Vector3<T>& init, const Vector3<T>& term);

        T Dot(const Vector3<T>& rhs) const;                 // computes the dot product between this and v2

        Vector3<T>& Cross(const Vector3<T>& lhs, const Vector3<T>& rhs); // sets this vector to the cross product of v1 and v2
        const Vector3<T> Cross(const Vector3<T>& rhs) const;                     // returns the cross product of this and v2

        const Vector3<T> operator + (const Vector3<T> &rhs) const;        // Adds two vectors. use Add() or += instead if possible
        const Vector3<T> operator - (const Vector3<T> &rhs) const;        // subtracts two vectors. use Sub() or -= instead if possible
        T operator * (const Vector3<T> &rhs) const;        // computes the dot product between this and rhs
        const Vector3<T> operator * (T rhs) const;             // multiplies this by a scalar. use Mul() or *= instead if possible
        const Vector3<T> operator * (const Matrix33<T>& rhs) const;
        const Vector3<T> operator / (T rhs) const;             // divides this by a scalar. use Div() or /= instead if possible

        Vector3<T>& operator += (const Vector3<T> &rhs);     // does an in-place addition with rhs
        Vector3<T>& operator -= (const Vector3<T> &rhs);     // does an in-place subtraction with rhs
        Vector3<T>& operator *= (T rhs);          // does an in-place scalar multiplication
        Vector3<T>& operator *= (const Matrix33<T>& mat);          // does an in-place matrix multiplication
        Vector3<T>& operator /= (T rhs);          // does an in-place scalar division

        T LengthSquared(void) const;
        T Length(void) const;

        Vector3<T>& NegationOf(const Vector3<T>& other);
        Vector3<T>& Negate(void);

        Vector3<T>& UnitOf(const Vector3<T>& other);
        Vector3<T>& Normalize(void);

        const Matrix33<T> Tensor(const Vector3<T> &rhs) const;  // computes the outer product of this and rhs
        void GetOrthonormalBasis(Vector3<T>& right, Vector3<T>& up) const; // makes u and v into an orthonormal basis of this. Assumes this is unit length

        T CosAngleBetween(const Vector3<T>& vect) const; // computes the cosine of the angle between this vector and vect
        T AngleBetween(const Vector3<T>& vect) const;     // computes the angle between this vector and vect

        T DistanceSquared(const Vector3<T>& v2) const; // finds the square of the distance between this and v2
        T Distance(const Vector3<T>& v2) const;          // finds the distance between this and v2

        bool operator == (const Vector3<T>& rhs) const;          // compares this to rhs and returns if they are equal
        bool operator != (const Vector3<T>& rhs) const;          // compares this to rhs and returns if they are not equal

        // Accessors
        T& operator [] (int32 i)                { return pntr[i]; }
        T  operator [] (int32 i) const        { return pntr[i]; }

        friend std::ostream& operator <<(std::ostream& out, const Vector3<T>& source) {     // outputs the vector to the streams
            out <<"<" <<source.x <<", " <<source.y <<", " <<source.z <<"> ";
            return out;
        }

        friend const Vector3<T> operator * (T scalar, const Vector3<T>& rhs) {         // returns a scalar multiplied by this vector 
            return Vector3<T>(scalar * rhs.x, scalar * rhs.y, scalar * rhs.z);
        }

        friend const Vector3<T> operator - (const Vector3<T>& v) {                                 // returns the negation of this vector
            return Vector3<T>(-v.x, -v.y, -v.z);
        }
    };

    template<class T> class Vector4 {
    public:
        // data
        union {
            struct {
                T x, y, z, w;
            };

            struct {
                T r, g, b, a;
            };

            T pntr[4];
        };

        // methods
        // constructors
        Vector4(void);
        Vector4(const Vector4<T>& other);

        explicit Vector4(const Vector2<T>& other, T _z, T _w);
        explicit Vector4(const Vector2<T>& other, T _z);
        explicit Vector4(const Vector3<T>& other, T _w);
        explicit Vector4(const Vector3<T>& other);

        // set from floats
        Vector4(T _x, T _y, T _z, T _w);
        Vector4(T _x, T _y, T _z);
        Vector4(const T _v[4]);

        ~Vector4(void);

        // setting methods
        Vector4<T>& operator=(const Vector4<T>& other);

        void Set(const Vector2<T>& other, T _z, T _w);
        void Set(const Vector2<T>& other, T _z);
        void Set(const Vector3<T>& other, T _w);
        void Set(const Vector3<T>& other);
        void Set(const Vector4<T>& other);

        // set from floats
        void Set(T _x, T _y, T _z, T _w);
        void Set(T _x, T _y, T _z);
        void Set(const T _v[4]);

        // clean up
        void Zero(void);
        void Clean(void);

        Vector4<T>& Add(const Vector4<T>&  lhs, const Vector4<T>& rhs);
        Vector4<T>& Sub(const Vector4<T>&  lhs, const Vector4<T>& rhs);
        Vector4<T>& Mul(const Vector4<T>&  lhs, T rhs);
        Vector4<T>& Mul(const Vector4<T>&  lhs, const Matrix44<T>& rhs);
        Vector4<T>& Mul(const Matrix44<T>& lhs, const Vector4<T>& rhs);
        Vector4<T>& Div(const Vector4<T>&  lhs, T rhs);
        Vector4<T>& Build(const Vector4<T>& init, const Vector4<T>& term);

        T                 Dot  (const Vector4<T>& rhs) const;                 // computes the dot product between this and v2

        Vector4<T>& Cross(const Vector4<T>& lhs, const Vector4<T>& rhs); // sets this vector to the cross product of v1 and v2
        const Vector4<T>  Cross(const Vector4<T>& rhs) const;                     // returns the cross product of this and v2

        const Vector4<T>  operator + (const Vector4<T>  &rhs) const;        // Adds two vectors. use Add() or += instead if possible
        const Vector4<T>  operator - (const Vector4<T>  &rhs) const;        // subtracts two vectors. use Sub() or -= instead if possible
                T                 operator * (const Vector4<T>  &rhs) const;        // computes the dot product between this and rhs
        const Vector4<T>  operator * (        T                  rhs) const;             // multiplies this by a scalar. use Mul() or *= instead if possible
        const Vector4<T>  operator * (const Matrix44<T> &rhs) const;
        const Vector4<T>  operator / (        T                  rhs) const;             // divides this by a scalar. use Div() or /= instead if possible

        Vector4<T>& operator += (const Vector4<T>  &rhs);     // does an in-place addition with rhs
        Vector4<T>& operator -= (const Vector4<T>  &rhs);     // does an in-place subtraction with rhs
        Vector4<T>& operator *= (        T                  rhs);          // does an in-place scalar multiplication
        Vector4<T>& operator *= (const Matrix44<T> &mat);          // does an in-place matrix multiplication
        Vector4<T>& operator /= (        T                  rhs);          // does an in-place scalar division

        T                 LengthSquared(void) const;
        T                 Length         (void) const;

        Vector4<T>& NegationOf(const Vector4<T>& other);
        Vector4<T>& Negate(void);

        Vector4<T>& UnitOf(const Vector4<T>& other);
        Vector4<T>& Normalize(void);

        const Matrix44<T> Tensor(const Vector4<T> &rhs) const;  // computes the outer product of this and rhs
        void GetOrthonormalBasis(Vector4<T>& right, Vector4<T>& up) const; // makes u and v into an orthonormal basis of this. Assumes this is unit length

        T CosAngleBetween(const Vector4<T>& vect) const; // computes the cosine of the angle between this vector and vect
        T AngleBetween(const Vector4<T>& vect) const;     // computes the angle between this vector and vect

        T DistanceSquared(const Vector4<T>& v2) const; // finds the square of the distance between this and v2
        T Distance(const Vector4<T>& v2) const;          // finds the distance between this and v2

        Vector4<T>& Homogenize(void);
        Vector4<T>& HomogenizedOf(const Vector4<T>& v2);

        bool operator == (const Vector4<T>& rhs) const;          // compares this to rhs and returns if they are equal
        bool operator != (const Vector4<T>& rhs) const;          // compares this to rhs and returns if they are not equal

        // Accessors
        T& operator [] (int32 i)                { return pntr[i]; }
        T  operator [] (int32 i) const        { return pntr[i]; }

        friend std::ostream& operator <<(std::ostream& out, const Vector4<T>& source) {     // outputs the vector to the streams
            out <<"<" <<source.x <<", " <<source.y <<", " <<source.z <<", " <<source.w <<"> ";
            return out;
        }

        friend const Vector4<T> operator * (T scalar, const Vector4<T>& rhs) {         // returns a scalar multiplied by this vector 
            return Vector4<T>(rhs) *= scalar;
        }

        friend const Vector4<T> operator - (const Vector4<T>& v){                                 // returns the negation of this vector
            return Vector4<T>(-v.x, -v.y, -v.z, -v.w);
        }
    };

    template<class T> class VectorN {
    public:
        T* pntr;
        int32 n;

        // methods
        // constructors
        VectorN();
        VectorN(int32 _n, T value);
        VectorN(const VectorN<T>& other);

        explicit VectorN(const Vector2<T>& other);
        explicit VectorN(const Vector3<T>& other);
        explicit VectorN(const Vector4<T>& other);

        // set from floats
        VectorN(const T* _v, int32 _n);
        explicit VectorN(int32 _n);

        ~VectorN(void);

        // size setting methods
        void SetSize(int32 nSlots); // sets the size of the array
        void SetSize(int32 nSlots, T value); // sets the size of the array
        void Clear(void);              // clears the memory used by this vector

        // setting methods
        VectorN<T>& operator=(const VectorN<T>& other);

        void Set(const Vector2<T>& other);
        void Set(const Vector3<T>& other);
        void Set(const Vector4<T>& other);
        void Set(const VectorN<T>& other);

        // set from floats
        void Set(const T* _v, int32 n);

        // clean up and fast implementation
        void Zero(void);
        void Clean(void);

        VectorN<T>& Add(const VectorN<T>&  lhs, const VectorN<T>& rhs);
        VectorN<T>& Sub(const VectorN<T>&  lhs, const VectorN<T>& rhs);
        VectorN<T>& Mul(const VectorN<T>&  lhs, T rhs);
        VectorN<T>& Mul(const VectorN<T>&  lhs, const MatrixMN<T>& rhs);
        VectorN<T>& Mul(const MatrixMN<T>& lhs, const VectorN<T>& rhs);
        VectorN<T>& Div(const VectorN<T>&  lhs, T rhs);
        VectorN<T>& Build(const VectorN<T>& init, const VectorN<T>& term);

        T Dot(const VectorN<T>& rhs) const;                 // computes the dot product between this and v2

        const VectorN<T> operator + (const VectorN<T>& rhs) const;        // Adds two vectors. use Add() or += instead if possible
        const VectorN<T> operator - (const VectorN<T>& rhs) const;        // subtracts two vectors. use Sub() or -= instead if possible
        T operator * (const VectorN<T> &rhs) const;        // computes the dot product between this and rhs
        const VectorN<T> operator * (T rhs) const;             // multiplies this by a scalar. use Mul() or *= instead if possible
        const VectorN<T> operator * (const MatrixMN<T>& rhs) const;
        const VectorN<T> operator / (T rhs) const;             // divides this by a scalar. use Div() or /= instead if possible

        VectorN<T>& operator += (const VectorN<T>& rhs);     // does an in-place addition with rhs
        VectorN<T>& operator -= (const VectorN<T>& rhs);     // does an in-place subtraction with rhs
        VectorN<T>& operator *= (T rhs);          // does an in-place scalar multiplication
        VectorN<T>& operator *= (const MatrixMN<T>& mat);          // does an in-place matrix multiplication
        VectorN<T>& operator /= (T rhs);          // does an in-place scalar division

        T MagnitudeSquared(void) const; // returns the squared 2-norm of this vector
        T Magnitude(void) const;

        T MaxElement(void) const; // returns the value of the maximum element
        int32    MaxIndex(void) const; // returns the index of the maximum element

        VectorN<T>& NegationOf(const VectorN<T>& other);
        VectorN<T>& Negate(void);

        VectorN<T>& UnitOf(const VectorN<T>& other);
        VectorN<T>& Normalize(void);

        const MatrixMN<T> Tensor(const VectorN<T> &rhs) const;  // computes the outer product of this and rhs

        T CosAngleBetween(const VectorN<T>& vect) const; // computes the cosine of the angle between this vector and vect
        T AngleBetween(const VectorN<T>& vect) const;     // computes the angle between this vector and vect

        T DistanceSquared(const VectorN<T>& v2) const; // finds the square of the distance between this and v2
        T Distance(const VectorN<T>& v2) const;          // finds the distance between this and v2

        bool operator == (const VectorN<T>& rhs) const;          // compares this to rhs and returns if they are equal
        bool operator != (const VectorN<T>& rhs) const;          // compares this to rhs and returns if they are not equal

        // Accessors
        T& operator [] (int32 i)                { return pntr[i]; }
        T  operator [] (int32 i) const        { return pntr[i]; }

        friend std::ostream& operator <<(std::ostream& out, const VectorN<T>& source) {     // outputs the vector to the streams
            out <<"<";
            if (source.n> 0) {
                out <<source.pntr[0];
                for(int32 ind = 1; ind <source.n; ++ind) {
                    out <<", " <<source.pntr[ind];
                }
            }

            out <<">";

            return out;
        }

        friend const VectorN<T> operator * (T scalar, const VectorN<T>& rhs) {         // returns a scalar multiplied by this vector 
            return VectorN<T>(rhs) *= scalar;
        }

        friend const VectorN<T> operator - (const VectorN<T>& pntr) {                                 // returns the negation of this vector
            VectorN<T> ans(pntr.n);
            for(int32 ind = 0; ind <pntr.n; ++ind) {
                ans.pntr[ind] = -pntr.pntr[ind];
            }

            return ans;
        }
    };

    template<class T> class Matrix22 {
    public:
        union {
            struct {
                T m00, m10;
                T m01, m11;
            };
            T mat[2][2];
            T pntr[4];
        };

        // methods
        // constructors
        Matrix22(void);
        Matrix22(const Matrix22<T>& other);
        explicit Matrix22(const Matrix33<T>& other);
        explicit Matrix22(const Matrix44<T>& other);

        Matrix22(const T v[4]); // Set the matrix to the values in v. v must be in column major form
        Matrix22(T _00, T _01,
            T _10, T _11);

        ~Matrix22(void);

        Matrix22<T>& operator = (const Matrix22<T>& other);

        // set the matrix to a type other than Matrix22<T>
        void Set(const Matrix22<T>& other);
        void Set(const Matrix33<T>& other);
        void Set(const Matrix44<T>& other);
        void Set(const T v[4]);
        void Set(T _00, T _01,
            T _10, T _11);

        void SetRow(const Vector2<T>& vals, int32 row);
        void SetRows(const Vector2<T>& row0, const Vector2<T>& row1);

        void SetColumn(const Vector2<T>& vals, int32 col);
        void SetColumns(const Vector2<T>& col0, const Vector2<T>& col1);

        Vector2<T> GetRow(int32) const;
        void GetRows(Vector2<T>& row0, Vector2<T>& row1) const;

        Vector2<T> GetColumn(int32) const;
        void GetColumns(Vector2<T>& col0, Vector2<T>& col1) const;

        void Clean(void);
        void Zero(void);
        void Identity(void);

        // Creates a rotation matrix which rotates about the z axis
        void Rotation(T ang);
        void Reflection(T ang);

        void TransposeOf(const Matrix22<T>& mat);
        void Transpose(void);

        bool InverseOf(const Matrix22<T>& mat);
        bool Invert(void);

        T Determinant(void) const;
        T Trace(void) const;

        Matrix22<T>& Tensor(const Vector2<T>& v1, const Vector2<T>& v2);

        const Matrix22<T> operator * (const Matrix22<T> &rhs) const;
        const Matrix22<T> operator * (        T    rhs) const;
        const Matrix22<T> operator / (        T    rhs) const;
        const Matrix22<T> operator + (const Matrix22<T> &rhs) const;
        const Matrix22<T> operator - (const Matrix22<T> &rhs) const;

        Matrix22<T>& operator *= (const Matrix22<T> &rhs);
        Matrix22<T>& operator *= (        T    rhs);
        Matrix22<T>& operator /= (        T    rhs);
        Matrix22<T>& operator += (const Matrix22<T> &rhs);
        Matrix22<T>& operator -= (const Matrix22<T> &rhs);

        // Methods used to replace the = sign. Use as ans = m1 op m2
        Matrix22<T>& Add(const Matrix22<T>& m1, const Matrix22<T>& m2);
        Matrix22<T>& Sub(const Matrix22<T>& m1, const Matrix22<T>& m2);
        Matrix22<T>& Mul(const Matrix22<T>& m1, const Matrix22<T>& m2);
        Matrix22<T>& Mul(const Matrix22<T>& m1,         T f);
        Matrix22<T>& Div(const Matrix22<T>& m1,         T f);

        Vector2<T>  operator *(const Vector2<T> &rhs) const;

        T& operator()(int32 row, int32 col)       { return pntr[row + (col << 1)]; }
        T  operator()(int32 row, int32 col) const { return pntr[row + (col << 1)]; }

        bool operator == (const Matrix22<T>& other) const;
        bool operator != (const Matrix22<T>& other) const;

        bool GetSymEigenInfo(Matrix22<T>& vects, Vector2<T>& val) const; // Get the eigenvalues and eigenvectors of the matrix. Will not work for integer types
        bool SymEigenDecompose(Vector2<T>& ew);

        friend std::ostream& operator<<(std::ostream& out, const Matrix22<T>& m) {
            out <<"{ { " <<m.m00 <<", " <<m.m01 <<" }, " \
                "{ " <<m.m10 <<", " <<m.m11 <<" } }";

            return out;
        }

        friend const Matrix22<T> operator * (T scalar, const Matrix22<T>& mat) {
            return Matrix22<T>(mat) *= scalar;
        }

        friend const Matrix22<T> operator - (const Matrix22<T>& mat) {
            return Matrix22<T>(-mat.m00, -mat.m01,
                                         -mat.m10, -mat.m11);
        }

    private:
        void Tridiagonalize(Vector2<T>& diag, Vector2<T>& subd);
        bool QLFactorize(Vector2<T>& diag, Vector2<T>& subd);
        void SortEigenInfo(Vector2<T>& diag);
    };

    template<class T> class Matrix33 {
    public:
        union {
            struct {
                T m00, m10, m20;
                T m01, m11, m21;
                T m02, m12, m22;
            };
            T mat[3][3];
            T pntr[9];
        };

        Matrix33(void);
        Matrix33(const Matrix33<T>& other);

        // Set the matrix to the values in v. v must be in column major form
        Matrix33(const T v[9]);
        Matrix33(T _00, T _01, T _02,
            T _10, T _11, T _12,
            T _20, T _21, T _22);

        explicit Matrix33(const Matrix22<T>& other);
        explicit Matrix33(const Matrix44<T>& other);
        explicit Matrix33(const Quaternion<T>& q);
        ~Matrix33(void);

        Matrix33<T>& operator = (const Matrix33<T>& other);

        // Set the matrix to the values in v. v must be in column major form
        void Set(const T v[9]);
        void Set(T _00, T _01, T _02,
            T _10, T _11, T _12, 
            T _20, T _21, T _22);

        void Set(const Matrix22<T>& other);
        void Set(const Matrix33<T>& other);
        void Set(const Matrix44<T>& other);
        void Set(const Quaternion<T>& q);

        void SetRow(const Vector3<T>& vals, int32 row);
        void SetRows(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2);

        void SetColumn(const Vector3<T>& vals, int32 col);
        void SetColumns(const Vector3<T>& col0, const Vector3<T>& col1, const Vector3<T>& col2);

        Vector3<T> GetRow(int32) const;
        void GetRows(Vector3<T>& row0, Vector3<T>& row1, Vector3<T>& row2) const;

        Vector3<T> GetColumn(int32) const;
        void GetColumns(Vector3<T>& col0, Vector3<T>& col1, Vector3<T>& col2) const;

        void Clean(void);
        void Zero(void);
        void Identity(void);

        void Rotation(T angle, T x, T y, T z);
        void Scale(T x, T y, T z);
        void Shear(const Vector3<T>& normal, const Vector3<T>& shear);

        void TransposeOf(const Matrix33<T>& mat);
        void Transpose(void);

        bool InverseOf(const Matrix33<T>& mat);
        bool Invert(void);

        T Determinant(void) const;
        T Trace(void) const;

        Matrix33<T>& Tensor(const Vector3<T>& v1, const Vector3<T>& v2);

        const Matrix33<T> operator *(const Matrix33<T> &rhs) const;
        const Matrix33<T> operator *(        T    rhs) const;
        const Matrix33<T> operator /(        T    rhs) const;
        const Matrix33<T> operator +(const Matrix33<T> &rhs) const;
        const Matrix33<T> operator -(const Matrix33<T> &rhs) const;

        Matrix33<T>& operator *=(const Matrix33<T> &rhs);
        Matrix33<T>& operator *=(        T    rhs);
        Matrix33<T>& operator /=(        T    rhs);
        Matrix33<T>& operator +=(const Matrix33<T> &rhs);
        Matrix33<T>& operator -=(const Matrix33<T> &rhs);

        // Methods used to replace the = sign. Use as ans = m1 op m2
        Matrix33<T>& Add(const Matrix33<T>& m1, const Matrix33<T>& m2);
        Matrix33<T>& Sub(const Matrix33<T>& m1, const Matrix33<T>& m2);
        Matrix33<T>& Mul(const Matrix33<T>& m1, const Matrix33<T>& m2);
        Matrix33<T>& Mul(const Matrix33<T>& m1,         T    f);
        Matrix33<T>& Div(const Matrix33<T>& m1,         T    f);

        const Vector3<T>  operator * (const Vector3<T>  &rhs) const;

        T& operator()(int32 row, int32 col)         { return pntr[row + col * 3]; }
        T  operator()(int32 row, int32 col) const { return pntr[row + col * 3]; }

        bool operator==(const Matrix33<T>& other) const;
        bool operator!=(const Matrix33<T>& other) const;

        bool GetSymEigenInfo(Matrix33<T>& eVects, Vector3<T>& eValues) const;
        bool SymEigenDecompose(Vector3<T>& eValues);

        friend std::ostream& operator<<(std::ostream& out, const Matrix33<T>& m) {
            out <<"{ { " <<m.m00 <<", " <<m.m01 <<", " <<m.m02 <<" }, " \
                "{ " <<m.m10 <<", " <<m.m11 <<", " <<m.m12 <<" }, " \
                "{ " <<m.m20 <<", " <<m.m21 <<", " <<m.m22 <<" } }";
            return out;
        }

        friend const Matrix33<T> operator * (T scalar, const Matrix33<T>& mat) {
            return Matrix33<T>(mat) *= scalar;
        }

        friend const Matrix33<T> operator - (const Matrix33<T>& mat) {
            return Matrix33<T>(mat) *= -1;
        }

    private:
        void Tridiagonalize(Vector3<T>& diag, Vector3<T>& subd);
        bool QLFactorize(Vector3<T>& diag, Vector3<T>& subd);
        void SortEigenInfo(Vector3<T>& diag);
    };

    template<class T> class Matrix44 {
    public:
        union {
            struct {
                T m00, m10, m20, m30;
                T m01, m11, m21, m31;
                T m02, m12, m22, m32;
                T m03, m13, m23, m33;
            };
            T mat[4][4];
            T pntr[16];
        };

        // Constructors
        Matrix44(void);

        Matrix44(T v[16]);
        Matrix44(T _00, T _01, T _02, T _03,
            T _10, T _11, T _12, T _13, 
            T _20, T _21, T _22, T _23, 
            T _30, T _31, T _32, T _33);


        explicit Matrix44(const Matrix22<T>& other);
        explicit Matrix44(const Matrix33<T>& other);
        Matrix44(const Matrix44<T>& other);
        explicit Matrix44(const Quaternion<T>& q);

        // Destructor
        ~Matrix44(void);

        // Assignment methods
        Matrix44<T>& operator = (const Matrix44<T>& other);

        void Set(T v[16]);
        void Set(T _00, T _01, T _02, T _03,
            T _10, T _11, T _12, T _13, 
            T _20, T _21, T _22, T _23, 
            T _30, T _31, T _32, T _33);

        void Set(const Matrix22<T>& other);
        void Set(const Matrix33<T>& other);
        void Set(const Matrix44<T>& other);
        void Set(const Quaternion<T>& q);

        // General utility
        // Cleans this matrix. Sets all close to zero values to zero
        void Clean(void);
        void Identity(void);
        void Zero(void);

        // Gets the row specified.
        const Vector4<T> GetRow(int32 row) const;
        void     GetRows(Vector4<T>& row0, Vector4<T>& row1, Vector4<T>& row2, Vector4<T>& row3) const;
        // Gets the column specified
        const Vector4<T> GetColumn(int32 col) const;
        void     GetColumns(Vector4<T>& col0, Vector4<T>& col1, Vector4<T>& col2, Vector4<T>& col3) const;

        void SetRow(const Vector4<T>& v, int32 row);
        void SetRow(T f1, T f2, T f3, T f4, int32 row);
        void SetRows(const Vector4<T>& row0, const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3);

        void SetColumn(const Vector4<T>& v, int32 col);
        void SetColumn(T f1, T f2, T f3, T f4, int32 col);
        void SetColumns(const Vector4<T>& col0, const Vector4<T>& col1, const Vector4<T>& col2, const Vector4<T>& col3);

        void Translation(T x, T y, T z);

        void Rotation(T angle, T x, T y, T z);
        void Rotation(const Quaternion<T>& q)                                    { Set(q); }
        void Rotation(const Matrix33<T>& mat)                                    { Set(mat); }

        void Scale(T x, T y, T z);
        void Shear(const Vector4<T>& axis, const Vector4<T>& shear);

        void LookAtGL(const Vector4<T>& from, const Vector4<T>& at, const Vector4<T>& up);
        void LookAtD3D(const Vector4<T>& from, const Vector4<T>& at, const Vector4<T>& up);
        
        bool InverseOf(const Matrix44<T>& mat);
        bool Invert(void);

        bool AffineInverseOf(const Matrix44<T>& mat);
        bool AffineInvert(void);

        void Transpose(void);
        void TransposeOf(const Matrix44<T> &mat);

        void PerspectiveGL(T fov, T aspect, T n, T f);
        void PerspectiveD3D(T fov, T aspect, T n, T f);

        void ObliquePerspGL(T left, T right, T bottom, T top, T n, T f);
        void ObliquePerspD3D(T left, T right, T bottom, T top, T n, T f);

        void OrthographicGL(T left, T right, T bottom, T top, T n, T f);
        void OrthographicD3D(T left, T right, T bottom, T top, T n, T f);

        void ObliqueParallelGL(T left, T right, T bottom, T top, T n, T f);
        void ObliqueParallelD3D(T left, T right, T bottom, T top, T n, T f);

        T Determinant(void) const;
        T Trace(void) const;

        Matrix44<T>& Tensor(const Vector4<T>& v1, const Vector4<T>& v2);

        const Matrix44<T> operator *(const Matrix44<T> &rhs) const;
        const Matrix44<T> operator *(         T    rhs) const;
        const Matrix44<T> operator /(         T    rhs) const;
        const Matrix44<T> operator +(const Matrix44<T> &rhs) const;
        const Matrix44<T> operator -(const Matrix44<T> &rhs) const;

        Matrix44<T>& operator *=(const Matrix44<T> &rhs);
        Matrix44<T>& operator *=(      T    rhs);
        Matrix44<T>& operator /=(      T    rhs);
        Matrix44<T>& operator +=(const Matrix44<T> &rhs);
        Matrix44<T>& operator -=(const Matrix44<T> &rhs);

        // Methods used to replace the = sign. Use as ans = m1 op m2
        Matrix44<T>& Add(const Matrix44<T>& m1, const Matrix44<T>& m2);
        Matrix44<T>& Sub(const Matrix44<T>& m1, const Matrix44<T>& m2);
        Matrix44<T>& Mul(const Matrix44<T>& m1, const Matrix44<T>& m2);
        Matrix44<T>& Mul(const Matrix44<T>& m1,         T    f);
        Matrix44<T>& Div(const Matrix44<T>& m1,         T    f);

        const Vector4<T>    operator * (const Vector4<T> &rhs) const;

        bool operator==(const Matrix44<T>& other);
        bool operator!=(const Matrix44<T>& other);

        T& operator()(int32 row, int32 col)       { return pntr[row + (col << 2)]; }
        T  operator()(int32 row, int32 col) const { return pntr[row + (col << 2)]; }

        bool GetSymEigenInfo(Matrix44<T>& eVects, Vector4<T>& eValues) const;
        bool SymEigenDecompose(Vector4<T>& ew);

        friend std::ostream& operator<<(std::ostream& out, const Matrix44<T>& m) {
            out <<"{ { " <<m.m00 <<", " <<m.m01 <<", " <<m.m02 <<", " <<m.m03 <<" }, " \
                "{ " <<m.m10 <<", " <<m.m11 <<", " <<m.m12 <<", " <<m.m13 <<" }, " \
                "{ " <<m.m20 <<", " <<m.m21 <<", " <<m.m22 <<", " <<m.m23 <<" }, " \
                "{ " <<m.m30 <<", " <<m.m31 <<", " <<m.m32 <<", " <<m.m33 <<" } }";
            return out;
        }

        friend const Matrix44<T> operator * (T scalar, const Matrix44<T>& mat) {
            return Matrix44<T>(mat) *= scalar;
        }

        friend const Matrix44<T> operator - (const Matrix44<T>& mat) {
            return Matrix44<T>(mat) *= -1;
        }

    private:
        void Tridiagonalize(Vector4<T>& diag, Vector4<T>& subd);
        bool QLFactorize(Vector4<T>& diag, Vector4<T>& subd);
        void SortEigenInfo(Vector4<T>& diag);
        T     DetPart(T c0, T c1, T c2, T c3, T c4, T c5, T c6, T c7, T c8) const;
    };

    // this is a general MxN ROW major matrix
    // no optimizations for smaller matrices
    template<class T> class MatrixMN {
    public:
        // Constructors
        MatrixMN(void);
        MatrixMN(int32 m, int32 n);

        MatrixMN(T*  v, int32 m, int32 n);
        MatrixMN(T** v, int32 m, int32 n);

        explicit MatrixMN(const Matrix22<T>& other);
        explicit MatrixMN(const Matrix33<T>& other);
        explicit MatrixMN(const Matrix44<T>& other);
        MatrixMN(const MatrixMN<T>& other);

        // Destructor
        ~MatrixMN(void);

        // Assignment methods
        MatrixMN<T>& operator = (const MatrixMN<T>& other);

        void SetSize(int32 _m, int32 _n);
        void Clear(void);

        void Set(T*  v, int32 _m, int32 _n);
        void Set(T** v, int32 _m, int32 _n);

        void Set(const Matrix22<T>& other);
        void Set(const Matrix33<T>& other);
        void Set(const Matrix44<T>& other);
        void Set(const MatrixMN<T>& other);

        // General utility
        // Cleans this matrix. Sets all close to zero values to zero
        void Clean(void);
        void Identity(void);
        void Zero(void);

        // Gets the row specified.
        const VectorN<T> GetRow(int32 row) const;

        // Gets the column specified
        const VectorN<T> GetColumn(int32 col) const;

        void SetRow(const VectorN<T>& v, int32 row);

        void SetColumn(const VectorN<T>& v, int32 col);

        bool InverseOf(const MatrixMN<T>& mat);
        bool Invert(void);

        void Transpose(void);
        void TransposeOf(const MatrixMN<T> &mat);

        T Determinant(void) const;
        T Trace(void) const;

        MatrixMN<T>& Tensor(const VectorN<T>& v1, const VectorN<T>& v2);

        const MatrixMN<T> operator *(const MatrixMN<T> &rhs) const;
        const MatrixMN<T> operator *(         T    rhs) const;
        const MatrixMN<T> operator +(const MatrixMN<T> &rhs) const;
        const MatrixMN<T> operator -(const MatrixMN<T> &rhs) const;

        MatrixMN<T>& operator *=(const MatrixMN<T> &rhs);
        MatrixMN<T>& operator *=(      T    rhs);
        MatrixMN<T>& operator +=(const MatrixMN<T> &rhs);
        MatrixMN<T>& operator -=(const MatrixMN<T> &rhs);

        // Methods used to replace the = sign. Use as ans = m1 op m2
        MatrixMN<T>& Add(const MatrixMN<T>& m1, const MatrixMN<T>& m2);
        MatrixMN<T>& Sub(const MatrixMN<T>& m1, const MatrixMN<T>& m2);
        MatrixMN<T>& Mul(const MatrixMN<T>& m1, const MatrixMN<T>& m2);
        MatrixMN<T>& Mul(const MatrixMN<T>& m1,         T    f);

        const VectorN<T>    operator * (const VectorN<T> &rhs) const;

        bool operator==(const MatrixMN<T>& other);
        bool operator!=(const MatrixMN<T>& other);

        T& operator()(int32 row, int32 col)         { return pntr[row * n + col]; }
        T  operator()(int32 row, int32 col) const { return pntr[row * n + col]; }

        bool GetSymEigenInfo(MatrixMN<T>& eVects, VectorN<T>& eValues) const;
        bool SymEigenDecompose(VectorN<T>& eValues);

        // linear solvers

        // LU Decomposition
        // places the LU decomposiiton of mat into this
        // rowPerms: a vector used to how the rows are swapped during the
        //  decomposition
        // dir: +/-1 when succesful; +1 for even number of row swaps,
        //  -1 for odd number of row swaps
        bool LUDecompositionOf(const MatrixMN& mat, int32* rowPerms, T& dir);

        // makes this matrix into its LU Decomposition
        // adapted from numerical recipes in C
        // rowPerms: a vector used to how the rows are swapped during the
        //  decomposition
        // dir: +/-1 when succesful; +1 for even number of row swaps,
        //  -1 for odd number of row swaps
        bool LUDecompose(int32* rowPerms, T& dir);

        // Uses LU Backsubstitution to solve for a linear set of equations
        // Assumes this is a LU Decompoosed matrix
        // b: vector on the right side of the equations
        // x: vector to store the solutions
        // rowPerms is the row swap vector from the LU Decomposition
        // x is reset if no solution exists
        bool LUBacksubstitute(const VectorN<T>& b, VectorN<T>& x, const int32* rowPerms) const;

        // Solves a linear set of equations using LU Decomposition
        // Assums nothing about this matrix
        // b: vector on the right side of the equations
        // x: vector to store the solution if one exists
        // x is reset if no solution exists
        bool LULinearSolver(const VectorN<T>& b, VectorN<T>& x) const;

        // need to add
        // Decomposes this using singular value decomposition
        // sigma: vector to hold the singular values
        // v: the v part of the matrix
        bool SVDecompose(VectorN<T>& sigma, MatrixMN<T>& v);

        // Uses Singular Value Backsubstitution to solve for a linear set of equations
        // Assumes this is a Singular Value Decomposed matrix
        // b: vector on the right side of the equations
        // x: vector to store the solutions
        // x is reset if no solution exists
        void SVBacksubstitute(const VectorN<T>& sigma, const MatrixMN<T>& v, const VectorN<T>& b, VectorN<T>& x) const;

        // Solves a linear set of equations using LU Decomposition
        // Assums nothing about this matrix
        // b: vector on the right side of the equations
        // x: vector to store the solution if one exists
        // x is reset if no solution exists
        bool SVDLinearSolver(const VectorN<T>& b, VectorN<T>& x) const;

        friend std::ostream& operator<<(std::ostream& out, const MatrixMN<T>& m) {
            out <<"{ ";
            for(int32 row = 0; row <m.m - 1; ++row) {
                out <<"{ " <<m(row, 0);
                for(int32 col = 1; col <m.n; ++col) {
                    out <<", " <<m(row, col);
                }
                out <<" }, ";
            }

            out <<"{ " <<m(row, 0);
            for(int32 col = 1; col <m.n; ++col) {
                out <<", " <<m(row, col);
            }
            out <<" } }";
            return out;
        }

        friend const MatrixMN<T> operator * (T scalar, const MatrixMN<T>& mat) {
            return MatrixMN<T>(mat) *= scalar;
        }

        friend const MatrixMN<T> operator - (const MatrixMN<T>& mat) {
            return MatrixMN<T>(mat) *= -1;
        }

    private:
        T* pntr;
        int32 m;
        int32 n;

        // Eigenvalue and Eigenvector methods
        void Tridiagonalize(VectorN<T>& diag, VectorN<T>& subd);
        bool QLFactorize(VectorN<T>& diag, VectorN<T>& subd);
        void SortEigenInfo(VectorN<T>& diag);
        T     Pythag(T a, T b) const;
    };

    template<class T> class SparseMatrix {
    };

    template<class T> class Quaternion {
    public:
        union {
            struct {
                T w, x, y, z;
            };

            struct {
                T w;
                Vector3<T> v;
            };
            T pntr[4];
        };

        Quaternion(void) {}    // default constructor. Does nothing
        Quaternion(const Quaternion<T>& other);                // copies other into this

        Quaternion(const T _v[4]); // set the quaternion to the values given
        Quaternion(T _w, T _x, T _y, T _z); // set the quaternion to the values given

        explicit Quaternion(const Matrix33<T>& rotMat);        // creates a quaternion using the given 3x3 rotation matrix
        explicit Quaternion(const Matrix44<T>& rotMat);        // creates a quaternion using the rotation part of the given 4x4 matrix
        explicit Quaternion(const Vector4<T>& vect);            // uses the full w, x, y, z components of the vector for the quaternion

        explicit Quaternion(const Vector3<T>& axis, T angle = 0);      // computes the quaternion using axis-angle representation
        Quaternion(const Vector3<T>& v1, const Vector3<T>& v2);    // computes the quaternion to rotate from v1 to v2

        ~Quaternion(void) {}                                        // destroys the quaternion

        Quaternion<T>& operator = (const Quaternion<T>& other);

        void Set(const T _v[4]); // set the quaternion to the values given
        void Set(T _w, T _x, T _y, T _z); // set the quaternion to the values given

        void Set(const Matrix33<T>& rotMat);        // creates a quaternion using the given 3x3 rotation matrix
        void Set(const Matrix44<T>& rotMat);        // creates a quaternion using the rotation part of the given 4x4 matrix
        void Set(const Vector4<T>& vect);            // uses the full w, x, y, z components of the vector for the quaternion

        void Set(const Vector3<T>& axis, T angle = 0);        // computes the quaternion using axis-angle representation
        void Set(const Vector3<T>& v1, const Vector3<T>& v2);    // computes the quaternion to rotate from v1 to v2

        void Clean(void);
        void Zero(void);
        void Identity(void);

        void MakeFromEuler(T fPitch, T fYaw, T fRoll); // Makes this quaternion from xyz order of the given angles
        void MakeFromAxes(const Vector3<T>& axis0, const Vector3<T>& axis1, const Vector3<T>& axis2);
        void Normalize(void);                                    // makes this quaternion unit length
        T     Selection(void) const                { return w; }    // returns the w component of the quaternion

        Quaternion<T>& Conjugate(void);                                // Returns the complex conjugate of this quaternion
        Quaternion<T>& ConjugateOf(const Quaternion<T>& q);            // turns this quaternion into the complex conjugate of q

        Quaternion<T>& Invert(void);                          // Returns the inverse of this quaternion
        Quaternion<T>& InverseOf(const Quaternion& q2);    // makes this quaternion the inverse of q2

        // find the axis angle representation that is equivalent to this quaternion
        void GetAxisAngle(Vector3<T>& axis, T& angle) const;
        // Gets the matrix representation of this quaternion
        Matrix44<T> GetMatrix(void) const;
        // Gets the rotate frame axes
        void GetAxes(Vector3<T>& axis0, Vector3<T>& axis1, Vector3<T>& axis2) const;

        // Gets the magnitude of this quaternion
        T Magnitude(void) const;

        T Dot(const Quaternion<T>& q2) const;                // finds the dot product between this and q2
        T Norm(void) const;                                    // returns the norm of this quaternion

        Vector4<T> Rotate(const Vector4<T>& v) const;            // Uses this quaternion to rotate v
        Vector3<T> Rotate(const Vector3<T>& v) const;            // Uses this quaternion to rotate v

        Quaternion<T> operator + (const Quaternion<T>& q) const;    // adds two quaternions
        Quaternion<T> operator - (const Quaternion<T>& q) const;    // subtracts two quaternion
        Quaternion<T> operator * (const Quaternion<T>& q) const;    // multiplies two quaternions
        Quaternion<T> operator * (        T                     f) const;                // multiplies a quaternion with a floating point value
        Vector4<T> operator * (const Vector4<T>&     v) const;    // multiplies a 4D vector by this quaternion
        Quaternion<T> operator / (        T                     f) const;                // divides this quaternion by a Float32

        Quaternion<T>& operator += (const Quaternion<T>& q);        // adds q to this quaternion
        Quaternion<T>& operator -= (const Quaternion<T>& q);        // subtracts q to this quaternion
        Quaternion<T>& operator *= (const Quaternion<T>& q);        // multiplies this quaternion by q
        Quaternion<T>& operator *= (T f);                    // multiplies this quaternion by f
        Quaternion<T>& operator /= (T f);                    // divides this quaternion by f

        // Methods used to replace the = sign. Use as ans = m1 op m2
        Quaternion<T>& Add(const Quaternion<T>& q1, const Quaternion<T>& q2);
        Quaternion<T>& Sub(const Quaternion<T>& q1, const Quaternion<T>& q2);
        Quaternion<T>& Mul(const Quaternion<T>& q1, const Quaternion<T>& q2);
        Quaternion<T>& Mul(const Quaternion<T>& q1,         T    f);
        Quaternion<T>& Div(const Quaternion<T>& q1,         T    f);

        const Vector3<T>  operator * (const Vector3<T>  &rhs) const;

        // accessors
        T& operator[](uint32 i)            { return pntr[i]; }    // acceses the quaternion uses an index
        T operator[](uint32 i) const     { return pntr[i]; }    // acceses the quaternion uses an index

        // comparators
        bool operator==(const Quaternion<T>& other) const;
        bool operator!=(const Quaternion<T>& other) const;

        bool IsZero(void) const;                                // Determines if this quaternion is zero
        bool IsUnit(void) const;                                // Determines if this quaternion has unit length
        bool IsIdentity(void) const;                            // Determines if this quaternion is the identity quaternion

        // interpolation
		Quaternion<T>& Lerp(const Quaternion<T>& start, const Quaternion<T>& end, T t);
		Quaternion<T>& Slerp(const Quaternion<T>& start, const Quaternion<T>& end, T t);
		Quaternion<T>& ApproxSlerp(const Quaternion<T>& start, const Quaternion<T>& end, T t);

        friend Quaternion<T> operator* (T f, const Quaternion<T>& q2) {    // multiplies q2 by f
            return Quaternion<T>(q2.w * f, q2.x * f, q2.y * f, q2.z * f);
        }
        
        friend std::ostream& operator <<(std::ostream& out, const Quaternion<T>& q) {    // outputs this quaternion to the console
            out <<"<" <<q.w <<" " <<q.v <<">";
            return out;
        }

        friend Quaternion<T> operator -(const Quaternion<T>& q) {    // negates the given quaternion
            return Quaternion<T>(-q.w, -q.x, -q.y, -q.z);
        }
    };

    // function prototypes
    template<class T> inline T DegToRad(T angle) { return angle * PI / 180.0f; }
    template<class T> inline T RadToDeg(T angle) { return angle * 180.0f * PIINV; }

    // returns if val is close to zero
    template<class T> inline bool IsZero(T val)  { return fabs(val) <= TINY; }
    template<class T> inline bool IsCloseToZero(T val) { return fabs(val) <= EPSILON; }
    template<class T> inline bool IsNotZero(T val) { return fabs(val) >  TINY; }

    template<class T> inline bool IsEqual(T f, T value)          { return fabs(f - value) <= TINY; }
    template<class T> inline bool IsNotEqual(T f, T value)      { return fabs(f - value) >  TINY; }

    template<class T> inline bool IsLessEqual(T f, T value)        { return f - value <= TINY; }
    template<class T> inline bool IsGreaterEqual(T f, T value)    { return value - f <= TINY; }

    template<class T> inline T     sqr(T val)                                         { return val * val; }
    template<class T> inline T     sign(T a, T b)                                    { return (b >= 0.0 ? fabs(a) : -fabs(a)); }

    template<class T> inline T fastfloor(T val) { return static_cast<float32>(static_cast<int32>(val < 0 ? val - 1 : val)); }

    // Floating point information
    inline int32 GetFPBits(float32 val)                            { return *(reinterpret_cast<const int32*>(&val)); }
    inline int32 GetFPExponent(float32 val)                      { return ((GetFPBits(val) & 0x7F800000)>> 23) - 127; }
    inline int32 GetFPSign(float32 val)                            { return (GetFPBits(val)>> 31); }
    inline int32 GetFPMantissa(float32 val)                      { return (GetFPBits(val) & 0x007FFFFF); }

    template<class T> inline T GaussianDistribution(T x, T y, T rho) {
        T g = 1.0f / sqrt(2.0f * PI * rho * rho);
        g *= expf(-(x * x + y * y) / (2 * rho * rho));

        return g;
    }

    inline uint32 intLog2(uint32 val) {
        int32 log = 0;
        while(val> 1) {
            ++log;
            val>>= 1;
        }
        return log;
    };

    typedef Vector2<uint8   > Vector2b;
    typedef Vector2<uint32 > Vector2u;
    typedef Vector2<int32   > Vector2i;
    typedef Vector2<float32> Vector2f;
    typedef Vector2<float64> Vector2d;

    typedef Vector3<uint8   > Vector3b;
    typedef Vector3<uint32 > Vector3u;
    typedef Vector3<int32   > Vector3i;
    typedef Vector3<float32> Vector3f;
    typedef Vector3<float64> Vector3d;

    typedef Vector4<uint8   > Vector4b;
    typedef Vector4<uint32 > Vector4u;
    typedef Vector4<int32   > Vector4i;
    typedef Vector4<float32> Vector4f;
    typedef Vector4<float64> Vector4d;

    typedef VectorN<uint8   > VectorNb;
    typedef VectorN<uint32 > VectorNu;
    typedef VectorN<int32   > VectorNi;
    typedef VectorN<float32> VectorNf;
    typedef VectorN<float64> VectorNd;

    typedef Matrix22<uint8   > Matrix22b;
    typedef Matrix22<uint32 > Matrix22u;
    typedef Matrix22<int32   > Matrix22i;
    typedef Matrix22<float32> Matrix22f;
    typedef Matrix22<float64> Matrix22d;

    typedef Matrix33<uint8   > Matrix33b;
    typedef Matrix33<uint32 > Matrix33u;
    typedef Matrix33<int32   > Matrix33i;
    typedef Matrix33<float32> Matrix33f;
    typedef Matrix33<float64> Matrix33d;

    typedef Matrix44<uint8   > Matrix44b;
    typedef Matrix44<uint32 > Matrix44u;
    typedef Matrix44<int32   > Matrix44i;
    typedef Matrix44<float32> Matrix44f;
    typedef Matrix44<float64> Matrix44d;

    typedef MatrixMN<uint8   > MatrixMNb;
    typedef MatrixMN<uint32 > MatrixMNu;
    typedef MatrixMN<int32   > MatrixMNi;
    typedef MatrixMN<float32> MatrixMNf;
    typedef MatrixMN<float64> MatrixMNd;

    typedef Quaternion<uint8   > Quaternionb;
    typedef Quaternion<uint32 > Quaternionu;
    typedef Quaternion<int32   > Quaternioni;
    typedef Quaternion<float32> Quaternionf;
    typedef Quaternion<float64> Quaterniond;
};

#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Vector4.hpp"
#include "VectorN.hpp"
#include "Matrix22.hpp"
#include "Matrix33.hpp"
#include "Matrix44.hpp"
#include "MatrixMN.hpp"
#include "Quaternion.hpp"