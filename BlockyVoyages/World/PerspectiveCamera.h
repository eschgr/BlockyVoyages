#pragma once

#include "Camera.h"

namespace BlockyVoyages {
namespace World {
// the camera doesn't know if it is an orthographic of a perspective camera
// Extend this class to fill in the projection matrix
class PerspectiveCamera : public Camera {
public:
    PerspectiveCamera(float32 fov, float32 aspect, float32 near, float32 far);
    ~PerspectiveCamera(void) {}

    void setAspect(float32 aspect);
    void setFOV(float32 fov);
    void setNearClipPlane(float32 near);
    void setFarClipPlane(float32 far);

private:
    float32 m_FOV;
    float32 m_aspect;
    float32 m_near;
    float32 m_far;
};
}
}