#include "stdafx.h"

#include "CharacterController.h"

#include "../World/Character.h"

namespace BlockyVoyages {
namespace Input {

CharacterController::CharacterController(World::Character* character)
    : m_character(character) {
}

// set the direction of travel
void CharacterController::setTravelDirection(const Vector3f& dir) {
    m_character->setTravelDirection(dir);
}
    
// rotate the character
void CharacterController::turnRightLeft(float32 angle) {
    m_character->turnRightLeft(angle);
}

void CharacterController::lookUpDown(float32 angle) {
    m_character->lookUpDown(angle);
}

// toggle running
void CharacterController::toggleRun() {
    m_character->toggleRun();
}

// use right hand
void CharacterController::useRightHand() {
    m_character->useRightHand();
}

// use left hand
void CharacterController::useLeftHand() {
    m_character->useLeftHand();
}

// set active item. Doesn't know what item would be switched to.
void CharacterController::setActiveItem(int32 slot) {
    m_character->setActiveItem(slot);
}

void CharacterController::setPosition(const Vector3f& pos) {
    m_character->setPosition(pos);
}

void CharacterController::setRotation(float32 heading, float32 headPitch) {
    m_character->setRotation(heading, headPitch);
}

void CharacterController::jump() {
    m_character->activateJump();
}

void CharacterController::toggleCrouch() {
    m_character->toggleCrouch();
}

}
}