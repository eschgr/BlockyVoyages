#pragma once

#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Input {

// A controller is an interface class that translates commands from the input
// processor and turns it into something the controlled object can understand.
class Controller {
public:
    virtual ~Controller() = 0 {}

    // set the direction of travel
    virtual void setTravelDirection(const Vector3f& dir) = 0;
    
    // rotate the character
    virtual void turnRightLeft(float32 angle) = 0;
    virtual void lookUpDown(float32 angle) = 0;

    // move the character to the specified position
    virtual void setPosition(const Vector3f& pos) = 0;
    virtual void setRotation(float32 heading, float32 headPitch) = 0;

    // toggle running
    virtual void toggleRun() = 0;

    virtual void jump() = 0;
    virtual void toggleCrouch() = 0;

    // use right hand
    virtual void useRightHand() = 0;
    // use left hand
    virtual void useLeftHand() = 0;

    // set active item. Doesn't know what item would be switched to.
    virtual void setActiveItem(int32 slot) = 0;
};

}
}