#pragma once

#include "../Math/MathMain.h"

#include <vector>

namespace BlockyVoyages {
namespace World {
// the camera doesn't know if it is an orthographic of a perspective camera
// Extend this class to fill in the projection matrix
class Camera {
public:
    virtual ~Camera() = 0 {}

    void resetView(const Vector3f& pos, const Quaternionf& rot);

    const Matrix44f& getProjectionMatrix() const { return m_projMat; }
    const Matrix44f& getViewMatrix() const { return m_viewMat; }

    Matrix44f getViewProjectionMatrix() const { return m_projMat * m_viewMat; }

    const Vector3f& getPosition() const { return m_pos; }
    const Quaternionf& getRotation() const { return m_rot; }
    Vector3f getLookAtDirection() const;

    void moveBy(const Vector3f& dir);
    void rotateBy(const Quaternionf& rot);
    void globalRotateBy(const Quaternionf& rot);

    void setPosition(const Vector3f& pos);
    void setRotation(const Quaternionf& rot);
        
    // assumes up is <0,1,0>. if (from - at) * <0,1,0> == 0, then it uses <1,0,0>
    void lookAt(const Vector3f& from, const Vector3f& at);
    // uses the position of the camera and looks to "at"
    void lookAt(const Vector3f& at) { lookAt(m_pos, at); }

    // returns the 6 planes that make up the view frustum
    std::vector<Vector4f> getViewFrustumPlanes() const;

protected:
    Camera();
    Camera& operator=(const Camera& other);

    void BuildViewMatrix();

    Matrix44f m_projMat;
    Matrix44f m_viewMat;

    Vector3f m_viewCenter;
    Vector3f m_pos;
    Quaternionf m_rot;
};
}
}
