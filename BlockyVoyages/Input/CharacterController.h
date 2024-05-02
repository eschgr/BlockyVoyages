#pragma once

#include "Controller.h"

namespace BlockyVoyages {
namespace World {
    class Character;
}

namespace Input {

class CharacterController : public Controller {
public:
    CharacterController(World::Character* character);
    virtual ~CharacterController() {}
    
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

    // use right hand
    virtual void useRightHand();
    // use left hand
    virtual void useLeftHand();

    // set active item. Doesn't know what item would be switched to.
    virtual void setActiveItem(int32 slot);

private:
    World::Character* m_character;
};

}
}