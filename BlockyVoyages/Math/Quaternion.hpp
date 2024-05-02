#include <math.h>

namespace BlockyVoyages {
    template<class T> inline Quaternion<T>::Quaternion(const Quaternion<T>& other) 
        : w(other.w), v(other.v)
    {}

    template<class T> inline Quaternion<T>::Quaternion(const T _v[4]) 
        : w(_v[0]), x(_v[1]), y(_v[2]), z(_v[3])
    {}

    template<class T> inline Quaternion<T>::Quaternion(T _w, T _x, T _y, T _z) 
        : w(_w), x(_x), y(_y), z(_z)
    {}

    template<class T> inline Quaternion<T>::Quaternion(const Matrix33<T>& rotMat) {
        Set(rotmat);
    }

    template<class T> inline Quaternion<T>::Quaternion(const Matrix44<T>& rotMat) {        // creates a quaternion using the rotation part of the given 4x4 matrix
        Set(rotmat);
    }

    template<class T> inline Quaternion<T>::Quaternion(const Vector4<T>& vect) 
        : w(vect.w), x(vect.x), y(vect.y), z(vect.z)
    {}

    template<class T> inline Quaternion<T>::Quaternion(const Vector3<T>& axis, T angle) 
        : w(cos(angle / 2)), v(axis)
    {
        v *= sin(angle / 2);
    }

    template<class T> inline Quaternion<T>::Quaternion(const Vector3<T>& v1, const Vector3<T>& v2) {    // computes the quaternion to rotate from v1 to v2
        Set(v1, v2);
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator=(const Quaternion<T>& other) {
        if (this != &other) {
            w = other.w;
            v = other.v;
        }

        return *this;
    }

    template<class T> inline void Quaternion<T>::Set(const T _v[4]) { // set the quaternion to the values given
        w = _v[0];
        v.Set(&_v[1]);
    }

    template<class T> inline void Quaternion<T>::Set(T _w, T _x, T _y, T _z) { // set the quaternion to the values given
        w = _w;
        v.Set(_x, _y, _z);
    }

    template<class T> inline void Quaternion<T>::Set(const Matrix33<T>& rotMat) {        // creates a quaternion using the given 3x3 rotation matrix
        T trace = rotMat.Trace();
        if (trace> 0.0f) {
            T s = sqrt(trace + 1.0f);
            w = s * 0.5;
            T recip = 0.5 / s;
            x = (rotMat.m12 - rotMat.m21) * recip;
            y = (rotMat.m20 - rotMat.m02) * recip;
            z = (rotMat.m01 - rotMat.m10) * recip;
        }
        else
        {
            uint32 i = 0;
            if (rotMat.m11> rotMat.m00) {
                i = 1;
            }

            if (rotMat.m22> rotMat.mv[i <<2]) {
                i = 2;
            }

            uint32 j = (i + 1) % 3;
            uint32 k = (j + 1) % 3;
            
            T s = sqrt(rotMat.pntr[i <<2] - rotMat.pntr[j <<2] - rotMat.pntr[k <<2] + 1.0f);
            T recip = 0.5f / s;
            w        = (rotMat(j, k) - rotMat(k, j)) * recip;
            v[i] = s * 0.5f;
            v[j] = (rotMat(i, j) + rotMat(j, i)) * recip;
            v[k] = (rotMat(i, k) + rotMat(k, i)) * recip;
        }
    }

    template<class T> inline void Quaternion<T>::Set(const Matrix44<T>& rotMat) {        // creates a quaternion using the rotation part of the given 4x4 matrix
        T trace = rotMat.m00 + rotMat.m11 + rotMat.m22;
		if (trace> 0.0f) {
			T s = sqrt(trace + 1.0f);
			w = s * 0.5f;
			T recip = 0.5f / s;
			x = (rotMat.m12 - rotMat.m21) * recip;
			y = (rotMat.m20 - rotMat.m02) * recip;
			z = (rotMat.m01 - rotMat.m10) * recip;
		}
		else {
			uint32 i = 0;
            if (rotMat.m11> rotMat.m00) {
				i = 1;
            }

            if (rotMat.m22> rotMat.pntr[i * 5]) {
				i = 2;
            }

			uint32 j = (i + 1) % 3;
			uint32 k = (j + 1) % 3;
			
			T s = sqrt(1.0f + rotMat.pntr[i * 5] - rotMat.pntr[j * 5] - rotMat.pntr[k * 5]);
            T recip = 0.5f / s;
			w =        (rotMat(j, k) - rotMat(k, j)) * recip;
            v[i] = 0.5f * s;
			v[j] = (rotMat(i, j) + rotMat(j, i)) * recip;
			v[k] = (rotMat(i, k) + rotMat(k, i)) * recip;
		}
    }

    template<class T> inline void Quaternion<T>::Set(const Vector4<T>& vect) {            // uses the full w, x, y, z components of the vector for the quaternion
        w = vect.w;
        x = vect.x;
        y = vect.y;
        z = vect.z;
    }

    template<class T> inline void Quaternion<T>::Set(const Vector3<T>& axis, T angle) {        // computes the quaternion using axis-angle representation
        w = cos(angle / 2);
        v.Mul(axis, sin(angle / 2));
    }

    template<class T> inline void Quaternion<T>::Set(const Vector3<T>& v1, const Vector3<T>& v2) {    // computes the quaternion to rotate from v1 to v2
        v.Cross(v1, v2);
        T s = sqrt(2.0f * (1.0f + v1 * v2));
        w = s / 2;
        v /= s;
    }

    template<class T> inline void Quaternion<T>::Clean(void) {
        if (IsCloseToZero(w)) {
            w = 0;
        }

        v.Clean();
    }

    template<class T> inline void Quaternion<T>::Zero(void) {
        w = x = y = z = 0;
    }

    template<class T> inline void Quaternion<T>::Identity(void) {
        w = 1;
        x = y = z = 0;
    }

    template<class T> inline void Quaternion<T>::MakeFromEuler(T fPitch, T fYaw, T fRoll) { // Makes this quaternion from xyz order of the given angles
        T cosx = cos(fPitch * 0.5f);
        T cosy = cos(fYaw * 0.5f);
        T cosz = cos(fRoll * 0.5f);

        T sinx = sin(fPitch * 0.5f);
        T siny = sin(fYaw * 0.5f);
        T sinz = sin(fRoll * 0.5f);

        w = cosx * cosy * cosz - sinx * siny * sinz;
        x = sinx * cosy * cosz + cosx * siny * sinz;
        y = cosx * siny * cosz - sinx * cosy * sinz;
        z = cosx * cosy * sinz + sinx * siny * cosz;
    }

    template<class T> inline void Quaternion<T>::MakeFromAxes(const Vector3<T>& axis0, const Vector3<T>& axis1, const Vector3<T>& axis2) {
        Matrix33 rotMat;
        rotMat.SetColumns(axis0, axis1, axis2);
        Set(rotMat);
    }

    template<class T> inline void Quaternion<T>::Normalize(void) {                                    // makes this quaternion unit length
        T len = Norm();
        if (IsNotZero(len)) {
            len = 1.0f / sqrt(len);
            w *= len;
            v *= len;
        }
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Conjugate(void) {                                // Returns the complex conjugate of this quaternion
        v *= -1.0f;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::ConjugateOf(const Quaternion<T>& q) {            // turns this quaternion into the complex conjugate of q
        w = other.w;
        v.Mul(other.v, -1.0f);
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Invert(void) {                          // Returns the inverse of this quaternion
        T norm = Norm();
        if (Voyages::IsZero(len)) {
            return *this;
        }

        norm = 1.0f / norm;
        w *= norm;
        v *= -norm;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::InverseOf(const Quaternion& q2) {
        *this = q2;
        return Invert();
    }

    template<class T> inline void Quaternion<T>::GetAxisAngle(Vector3<T>& axis, T& angle) const {    // find the axis angle representation that is equivalent to this quaternion
        if (IsLessEqual(w, static_cast<T>(1))) {
            angle = sqrt(1.0f - w * w);
            axis.Set(x, y, z);
            axis /= angle;

            angle = 2.0f * acos(w);
        }
        else {
            axis.Set(0, 0, 1);
            angle = 0.0f;
        }
    }

    template<class T> inline Matrix44<T> Quaternion<T>::GetMatrix(void) const {                        // Gets the matrix representation of this quaternion
        Matrix44<T> mat1( w, -z,  y, x,
                                    z,  w, -x, y,
                                  -y,  x,  w, z,
                                  -x, -y, -z, w);

        Matrix44<T> mat2( w, -z,  y, -x,
                                    z,  w, -x, -y,
                                  -y,  x,  w, -z,
                                    x,  y,  z,  w);
        Matrix44<T> res;

        res.Mul(mat1, mat2);

        // do some cleanup. If not done, numerical inaccuracy causes problems
        res.m03 = res.m13 = res.m23 = 0;
        res.m30 = res.m31 = res.m32 = 0;
        res.m33 = 1;

        return res;
    }

    template<class T> inline void Quaternion<T>::GetAxes(Vector3<T>& axis0, Vector3<T>& axis1, Vector3<T>& axis2) const {
        Matrix44<T> mat1( w, -z,  y, x,
                                    z,  w, -x, y,
                                  -y,  x,  w, z,
                                  -x, -y, -z, w);

        Matrix44<T> mat2( w, -z,  y, -x,
                                    z,  w, -x, -y,
                                  -y,  x,  w, -z,
                                    x,  y,  z,  w);
        Matrix44<T> res;
        res.Mul(mat1, mat2);

        // do some cleanup. If not done, numerical inaccuracy causes problems
        res.m03 = res.m13 = res.m23 = 0;
        res.m30 = res.m31 = res.m32 = 0;
        res.m33 = 1;

        axis0.Set(res.m00, res.m10, res.m20);
        axis1.Set(res.m01, res.m11, res.m21);
        axis2.Set(res.m02, res.m12, res.m22);
    }

    template<class T> inline T Quaternion<T>::Magnitude(void) const {                            // Gets the magnitude of this quaternion
        return sqrt(norm);
    }

    template<class T> inline T Quaternion<T>::Dot(const Quaternion<T>& q2) const {                // finds the dot product between this and q2
        return w * q2.w + x * q2.x + y * q2.y + z * q2.z;
    }

    template<class T> inline T Quaternion<T>::Norm(void) const {                                    // returns the norm of this quaternion
        return w * w + x * x + y * y + z * z;
    }

    template<class T> inline Vector4<T> Quaternion<T>::Rotate(const Vector4<T>& v) const {            // Uses this quaternion to rotate v
        T pMult = 1.0f - w * w;
        T vMult = 2.0f * (x * v.x + y * v.y + z * v.z);
        T crossMult = 2.0f * w;

        return Vector4<T>(pMult * v.x + vMult * x + crossMult * (y * v.z - z * v.y),
                                    pMult * v.y + vMult * y + crossMult * (z * v.x - x * v.z),
                                    pMult * v.z + vMult * z + crossMult * (x * v.y - y * v.x),
                                    v.w);
    }

    template<class T> inline Vector3<T> Quaternion<T>::Rotate(const Vector3<T>& v) const {            // Uses this quaternion to rotate v
        T pMult = 1.0f - w * w;
        T vMult = 2.0f * (x * v.x + y * v.y + z * v.z);
        T crossMult = 2.0f * w;

        return Vector3<T>(pMult * v.x + vMult * x + crossMult * (y * v.z - z * v.y),
                                    pMult * v.y + vMult * y + crossMult * (z * v.x - x * v.z),
                                    pMult * v.z + vMult * z + crossMult * (x * v.y - y * v.x));
    }

    template<class T> inline Quaternion<T> Quaternion<T>::operator + (const Quaternion<T>& q) const {    // adds two quaternions
        return Quaternion<T>(*this) += q;
    }

    template<class T> inline Quaternion<T> Quaternion<T>::operator - (const Quaternion<T>& q) const {    // subtracts two quaternions
        return Quaternion<T>(*this) -= q;
    }

    template<class T> inline Quaternion<T> Quaternion<T>::operator * (const Quaternion<T>& q) const {    // multiplies two quaternions
        return Quaternion<T>(w * q.w - x * q.x - y * q.y - z * q.z,
						                y * q.z - z * q.y + w * q.x + x * q.w,
						                z * q.x - x * q.z + w * q.y + y * q.w,
						                x * q.y - y * q.x + w * q.z + z * q.w);
    }

    template<class T> inline Quaternion<T> Quaternion<T>::operator * (        T                     f) const {                // multiplies a quaternion with a floating point value
        return Quaternion<T>(*this) *= f;
    }

    template<class T> inline Vector4<T>     Quaternion<T>::operator * (const Vector4<T>&     v) const {    // multiplies a 4D vector by this quaternion
        return Vector4<T>(y * v.z - z * v.y + w * v.x,
						            z * v.x - x * v.z + w * v.y,
						            x * v.y - y * v.x + w * v.z,
						            v.w);
    }

    template<class T> inline const Vector3<T>  Quaternion<T>::operator * (const Vector3<T>  &rhs) const {
        return Vector3<T>(y * v.z - z * v.y + w * v.x,
						            z * v.x - x * v.z + w * v.y,
						            x * v.y - y * v.x + w * v.z);
    }

    template<class T> inline Quaternion<T> Quaternion<T>::operator / (        T                     f) const {                // divides this quaternion by a float32
        return Quaternion<T>(*this) /= f;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator += (const Quaternion<T>& q) {        // adds q to this quaternion
        w += q.w;
        v += q.v;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator -= (const Quaternion<T>& q) {        // subtracts q to this quaternion
        w -= q.w;
        v -= q.v;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator *= (const Quaternion<T>& q) {        // multiplies this quaternion by q
        return (*this) = (*this) * q;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator *= (T f) {                    // multiplies this quaternion by f
        w *= f;
        v *= f;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::operator /= (T f) {                    // divides this quaternion by f
        T inv = 1.0f / f;
        w *= inv;
        v *= inv;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Add(const Quaternion<T>& q1, const Quaternion<T>& q2) {
        w = q1.w + q2.w;
        x = q1.x + q2.x;
        y = q1.y + q2.y;
        z = q1.z + q2.z;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Sub(const Quaternion<T>& q1, const Quaternion<T>& q2) {
        w = q1.w - q2.w;
        x = q1.x - q2.x;
        y = q1.y - q2.y;
        z = q1.z - q2.z;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Mul(const Quaternion<T>& q1, const Quaternion<T>& q2) {
        w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
        x = q1.y * q2.z - q1.z * q2.y + q1.w * q2.x + q1.x * q2.w;
		y = q1.z * q2.x - q1.x * q2.z + q1.w * q2.y + q1.y * q2.w;
		z = q1.x * q2.y - q1.y * q2.x + q1.w * q2.z + q1.z * q2.w;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Mul(const Quaternion<T>& q1,         T    f) {
        w = q1.w * f;
        x = q1.x * f;
        y = q1.y * f;
        z = q1.z * f;
        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Div(const Quaternion<T>& q1,         T    f) {
        T inv = 1.0f / f;
        w = q1.w * inv;
        x = q1.x * inv;
        y = q1.y * inv;
        z = q1.z * inv;
        return *this;
    }

    template<class T> inline bool Quaternion<T>::operator==(const Quaternion<T>& other) const {
        return this == &other || (w == other.w && x == other.x && y == other.y && z == other.z);
    }

    template<class T> inline bool Quaternion<T>::operator!=(const Quaternion<T>& other) const {
        return !(*this == other);
    }

    template<class T> inline bool Quaternion<T>::IsZero(void) const {
        return BlockyVoyages::IsZero(w) &&
                 BlockyVoyages::IsZero(x) &&
                 BlockyVoyages::IsZero(y) &&
                 BlockyVoyages::IsZero(z);
    }

    template<class T> inline bool Quaternion<T>::IsUnit(void) const {
        return IsEqual(Norm(), static_cast<T>(1));
    }

    template<class T> inline bool Quaternion<T>::IsIdentity(void) const {
        return BlockyVoyages::IsEqual(w, static_cast<T>(1)) &&
                 BlockyVoyages::IsZero(x) &&
                 BlockyVoyages::IsZero(y) &&
                 BlockyVoyages::IsZero(z);
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Lerp(const Quaternion<T>& start, const Quaternion<T>& end, T t) {
        // get cos of "angle" between quaternions
        T cosTheta = start.Dot(end);

        // initialize result
        Quaternion<T> temp(start);

        w = end.w * t;
        x = end.x * t;
        y = end.y * t;
        z = end.z * t;

        // if "angle" between quaternions is less than 90 degrees
        if (cosTheta>= EPSILON)
        {
            // use standard interpolation
            temp *= (1.0f - t);
        }
        else
        {
            // otherwise, take the shorter path
            temp *= (t - 1.0f);
        }

        w += temp.w;
        x += temp.x;
        y += temp.y;
        z += temp.z;

        return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::Slerp(const Quaternion<T>& start, const Quaternion<T>& end, T t) {
        // get cosine of "angle" between quaternions
		T cosTheta = start.Dot(end);
		T startInterp, endInterp;

		// if "angle" between quaternions is less than 90 degrees
		if (cosTheta>= EPSILON)
		{
			// if angle is greater than zero
			if ((1.0f - cosTheta)> EPSILON)
			{
				// use standard slerp
				T theta = acosf(cosTheta);
				T recipSinTheta = 1.0f / sinf(theta);

				startInterp = sinf((1.0f - t) * theta) * recipSinTheta;
				endInterp = sinf(t * theta) * recipSinTheta;
			}
			// angle is close to zero
			else
			{
				// use linear interpolation
				startInterp = 1.0f - t;
				endInterp = t;
			}
		}
		// otherwise, take the shorter route
		else
		{
			// if angle is less than 180 degrees
			if ((1.0f + cosTheta)> EPSILON)
			{
				// use slerp w/negation of start quaternion
				T theta = acosf(-cosTheta);
				T recipSinTheta = 1.0f/sinf(theta);

				startInterp = sinf((t - 1.0f) * theta) * recipSinTheta;
				endInterp = sinf(t * theta) * recipSinTheta;
			}
			// angle is close to 180 degrees
			else
			{
				// use lerp w/negation of start quaternion
				startInterp = t - 1.0f;
				endInterp = t;
			}
		}
	     
		w = startInterp * start.w + endInterp * end.w;
		x = startInterp * start.x + endInterp * end.x;
		y = startInterp * start.y + endInterp * end.y;
		z = startInterp * start.z + endInterp * end.z;

		return *this;
    }

    template<class T> inline Quaternion<T>& Quaternion<T>::ApproxSlerp(const Quaternion<T>& start, const Quaternion<T>& end, T t) {
        T cosTheta = start.Dot(end);

		// correct time by using cosine of angle between quaternions
		T factor = 1.0f - 0.7878088f * cosTheta;
		T k = 0.5069269f;
		factor *= factor;
		k *= factor;

		T b = 2 * k;
		T c = -3 * k;
		T d = 1 + k;

		t = t * (b * t + c) + d;

		// initialize result
		w = end.w * t;
		x = end.x * t;
		y = end.y * t;
		z = end.z * t;

		// if "angle" between quaternions is less than 90 degrees
		if (cosTheta>= EPSILON)
		{
			// use standard interpolation
			w += (1.0f - t) * start.w;
			x += (1.0f - t) * start.x;
			y += (1.0f - t) * start.y;
			z += (1.0f - t) * start.z;
		}
		else
		{
			// otherwise, take the shorter path
			w += (t - 1.0f) * start.w;
			x += (t - 1.0f) * start.x;
			y += (t - 1.0f) * start.y;
			z += (t - 1.0f) * start.z;
		}

		return *this;
    }
}