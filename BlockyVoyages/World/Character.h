#pragma once

#include "Map.h"
#include "../Graphics/VertexBuffer.h"
#include "../Input/CameraController.h"
#include "../Input/CharacterController.h"
#include "../Geometry/AABB.h"
#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace World {

// a character is an in game object that can be controlled. It has a controller
// associated with it that can be passed onto the input or AI engine.
struct CharacterStatistics {
    float32 runSpeed;
    float32 walkSpeed;
    float32 accelRate;
    float32 jumpHeight;
};

enum CharacterState {
    IDLE = 0x0,
    IN_AIR = 0x1,
    JUMPING = 0x3,
    // if not crouching or running, then walking
    CROUCHING = 0x4,
    RUNNING = 0x8
};

class Character {
public:
    Character(Map* Map);
    ~Character();

    void setPosition(const Vector3f& pos);
    void setRotation(float32 heading, float32 headPitch);

    const Vector3f& getPosition() const { return m_worldPos; }
    float32 getPitch() const { return m_pitch; }
    float32 getHeading() const { return m_heading; }

    // set the direction of travel
    void setTravelDirection(const Vector3f& dir);
    
    // rotate the character
    void turnRightLeft(float32 angle);
    void lookUpDown(float32 angle);

    // toggle running
    void toggleRun();

    // jumping and crouching
    void activateJump();
    void toggleCrouch();

    // use right hand
    void useRightHand();
    // use left hand
    void useLeftHand();

    // set active item. Doesn't know what item would be switched to.
    void setActiveItem(int32 slot);

    // attaches a camera to the player. Pass a nullptr to detach the camera
    // passing the controller allows the character to not need to know about
    // the movement model of the camera
    void attachCamera(Input::CameraController* camController);

    Input::CameraController* getAttachedCamera() { return m_attachedCamera; }

    Input::Controller* getController() { return &m_controller; }

    // update the position of the player.
    void update(float32 dt);

    // draw the character on screen
    void draw();
private:
    // CharacterController
    Input::CharacterController m_controller;

    // if a camera is attached to the player, pass all commands to the camera
    // as well. Only the movement commands will be passed as the camera only
    // contains positional information.
    Input::CameraController* m_attachedCamera;

    // Model
    Vector3f m_dims;
    Geometry::AABB m_aabb;

    // position information
    Vector3f m_worldPos;
    Vector3f m_vel;
    float32 m_heading;
    float32 m_pitch;

    // inventory

    // Character statistics. This is for the non-changing statistics. For the
    // changing status, see m_status.
    CharacterStatistics m_stats;

    // a bit field of all the states that the character may be in
    int32 m_stateFlags;

    // Information about the world
    Map* m_worldMap;
};

}
}