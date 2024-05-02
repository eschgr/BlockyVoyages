#include "stdafx.h"

#include "Camera.h"

namespace BlockyVoyages {
namespace World {
Camera::Camera(void) 
    : m_projMat(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1),
      m_viewMat(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1),
      m_viewCenter(0, 0, 0),
      m_pos(0, 0, 0),
      m_rot(1, 0, 0, 0)
{}

void Camera::resetView(const Vector3f& pos, const Quaternionf& rot) {
    m_pos = pos;
    m_rot = rot;

    // rebuild the matrix
    BuildViewMatrix();
}

Vector3f Camera::getLookAtDirection(void) const {
    Vector4f lookAtDir = m_viewMat.GetRow(2);
    return Vector3f(-lookAtDir.x, -lookAtDir.y, -lookAtDir.z);
}

void Camera::moveBy(const Vector3f& dir) {
    // set up the xform matrix
    m_viewMat.m03 -= dir.x;
    m_viewMat.m13 -= dir.y;
    m_viewMat.m23 -= dir.z;
        
    // find the position of the camera
    Matrix44f xform;
    xform.AffineInverseOf(m_viewMat);
    Vector4f xDir;
    xDir.Mul(xform, Vector4f(dir));
        
    m_pos.x += xDir.x;
    m_pos.y += xDir.y;
    m_pos.z += xDir.z;
}

void Camera::rotateBy(const Quaternionf& rot) {
    m_rot *= rot;

    Quaternionf tRot(rot);
    tRot.w *= -1;
    Matrix44f rotMat(tRot);
    m_viewMat = rotMat * m_viewMat;
}

void Camera::globalRotateBy(const Quaternionf& rot) {
    m_rot = rot * m_rot;

    // rebuild the matrix
    BuildViewMatrix();
}

void Camera::setPosition(const Vector3f& pos) {
    m_pos = pos;

    // rebuild the matrix
    BuildViewMatrix();
}

void Camera::setRotation(const Quaternionf& rot) {
    m_rot = rot;

    // rebuild the matrix
    BuildViewMatrix();
}
// assumes up is <0,1,0>. if (from - at) * <0,1,0> == 0, then it uses <1,0,0>
void Camera::lookAt(const Vector3f& from, const Vector3f& at) {
    Vector3f up(0, 1, 0);
    Vector3f dir;
    dir.Sub(from, at);

    if (IsZero(dir.x) && IsZero(dir.z)) {
        up.Set(1, 0, 0);
    }
    m_viewMat.LookAtGL(Vector4f(from, 1),
                     Vector4f(at, 1),
                     Vector4f(up));

    // get the information on the position and rotation
    m_pos = from;
    m_rot.Set(m_viewMat);
}

// returns the 6 planes that make up the view frustum
std::vector<Vector4f> Camera::getViewFrustumPlanes() const {
    // create the view projection matrix
    Matrix44f viewProjMat;
    viewProjMat.Mul(m_projMat, m_viewMat);
    std::vector<Vector4f> planes(6);

    int32 i = 0;
    
    // Near Frustum Plane
    // Add the third column to the fourth column to get the near plane.
	planes[i].x = viewProjMat.m20 + viewProjMat.m30;
	planes[i].y = viewProjMat.m21 + viewProjMat.m31;
	planes[i].z = viewProjMat.m22 + viewProjMat.m32;
	planes[i++].w = viewProjMat.m23 + viewProjMat.m33;

    // Left Frustum Plane
    // Add first column of the matrix to the fourth column
    planes[i].x = viewProjMat.m30 + viewProjMat.m00; 
	planes[i].y = viewProjMat.m31 + viewProjMat.m01;
	planes[i].z = viewProjMat.m32 + viewProjMat.m02;
	planes[i++].w = viewProjMat.m33 + viewProjMat.m03;

	// Right Frustum Plane
    // Subtract first column of matrix from the fourth column
	planes[i].x = viewProjMat.m30 - viewProjMat.m00; 
	planes[i].y = viewProjMat.m31 - viewProjMat.m01;
	planes[i].z = viewProjMat.m32 - viewProjMat.m02;
	planes[i++].w = viewProjMat.m33 - viewProjMat.m03;

	// Top Frustum Plane
    // Subtract second column of matrix from the fourth column
	planes[i].x = viewProjMat.m30 - viewProjMat.m10; 
	planes[i].y = viewProjMat.m31 - viewProjMat.m11;
	planes[i].z = viewProjMat.m32 - viewProjMat.m12;
	planes[i++].w = viewProjMat.m33 - viewProjMat.m13;

	// Bottom Frustum Plane
    // Add second column of the matrix to the fourth column
	planes[i].x = viewProjMat.m30 + viewProjMat.m10;
	planes[i].y = viewProjMat.m31 + viewProjMat.m11;
	planes[i].z = viewProjMat.m32 + viewProjMat.m12;
	planes[i++].w = viewProjMat.m33 + viewProjMat.m13;

	// Far Frustum Plane
    // Subtract third column of matrix from the fourth column
	planes[i].x = viewProjMat.m30 - viewProjMat.m20; 
	planes[i].y = viewProjMat.m31 - viewProjMat.m21;
	planes[i].z = viewProjMat.m32 - viewProjMat.m22;
	planes[i++].w = viewProjMat.m33 - viewProjMat.m23;

    // normalize the plane normals
    for (int i = 0; i < 6; ++i) {
        float32 len = sqrt(planes[i].x * planes[i].x +
                           planes[i].y * planes[i].y +
                           planes[i].z * planes[i].z);
        planes[i] /= len;
    }

    return planes;
}

void Camera::BuildViewMatrix(void) {
    m_viewMat.Translation(-m_pos.x, -m_pos.y, -m_pos.z);
    Quaternionf tRot(m_rot);
    tRot.w *= -1;
    m_viewMat = tRot.GetMatrix() * m_viewMat;
}
}
}