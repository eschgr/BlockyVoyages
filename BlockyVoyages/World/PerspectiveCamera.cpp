#include "stdafx.h"

#include "PerspectiveCamera.h"

namespace BlockyVoyages {
namespace World {
PerspectiveCamera::PerspectiveCamera(float32 fov, float32 aspect, float32 nearP, float32 farP) 
    : m_FOV(DegToRad(fov)),
      m_aspect(aspect),
      m_near(nearP),
      m_far(farP)
{
    m_projMat.PerspectiveGL(m_FOV, m_aspect, m_near, m_far);
}

void PerspectiveCamera::setAspect(float32 m_aspect) {
    m_aspect = m_aspect;

    m_projMat.PerspectiveGL(m_FOV, m_aspect, m_near, m_far);
}

void PerspectiveCamera::setFOV(float32 fov) {
    m_FOV = DegToRad(fov);

    m_projMat.PerspectiveGL(m_FOV, m_aspect, m_near, m_far);
}

void PerspectiveCamera::setNearClipPlane(float32 nearP) {
    m_near = nearP;

    m_projMat.PerspectiveGL(m_FOV, m_aspect, m_near, m_far);
}

void PerspectiveCamera::setFarClipPlane(float32 farP) {
    m_far = farP;

    m_projMat.PerspectiveGL(m_FOV, m_aspect, m_near, m_far);
}
}
}