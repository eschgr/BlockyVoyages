#pragma once

#include "Controller.h"
#include "../World/Camera.h"

namespace BlockyVoyages {
namespace Input {

class CameraController : public Controller {
public:
    CameraController(World::Camera* camera);
    virtual ~CameraController() {}

    // set the direction of travel
    virtual void setTravelDirection(const Vector3f& dir);
    
    // rotate the character
    virtual void turnRightLeft(float32 angle);
    virtual void lookUpDown(float32 angle);

    // move the character to the specified position
    virtual void setPosition(const Vector3f& pos);
    virtual void setRotation(float32 heading, float32 headPitch);

    // toggle running
    virtual void toggleRun();

    virtual void jump();
    virtual void toggleCrouch();

    // Cameras have no items or inventory, so they can't be used

    // use right hand
    virtual void useRightHand() {}
    // use left hand
    virtual void useLeftHand() {}

    // set active item. Doesn't know what item would be switched to.
    virtual void setActiveItem(int32 slot) {}

    void setTimeDelta(float32 delta) { m_dt = delta; }

private:
    CameraController();
    World::Camera* m_camera;
    bool m_running;
    float32 m_dt;

    static const float32 moveSpeed;
};

}
}