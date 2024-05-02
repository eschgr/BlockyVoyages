#include "stdafx.h"

#include "CameraController.h"

namespace BlockyVoyages {
namespace Input {

const float32 CameraController::moveSpeed = 2.0f;

CameraController::CameraController(World::Camera* camera)
    : m_camera(camera),
      m_running(false),
      m_dt(0.0f)
{}

// set the direction of travel
void CameraController::setTravelDirection(const Vector3f& dir) {
    Vector3f scaledDir;
    scaledDir.Mul(dir, moveSpeed * (m_running ? 10.0f : 1.0f) * m_dt);
    m_camera->moveBy(scaledDir);
}

// rotate the character
void CameraController::turnRightLeft(float32 angle) {
    m_camera->globalRotateBy(Quaternionf(Vector3f(0, 1, 0), angle));
}

void CameraController::lookUpDown(float32 angle) {
    m_camera->rotateBy(Quaternionf(Vector3f(1.0f, 0.0f, 0.0f), angle));
}

void CameraController::setPosition(const Vector3f& pos) {
    m_camera->setPosition(pos);
}

void CameraController::setRotation(float32 heading, float32 headPitch) {
    Quaternionf headingRotQuat, pitchRotQuat, rotQuat;
    headingRotQuat.MakeFromEuler(0.0f, heading, 0.0f);
    pitchRotQuat.MakeFromEuler(headPitch, 0.0f, 0.0f);
    rotQuat.Mul(headingRotQuat, pitchRotQuat);
    m_camera->setRotation(rotQuat);
}

// toggle running
void CameraController::toggleRun() {
    m_running = !m_running;
}

void CameraController::jump() {
    Vector3f scaledDir(0.0f, (m_running ? 10.0f : 1.0f) * m_dt, 0.0f);
    m_camera->moveBy(scaledDir);
}

void CameraController::toggleCrouch() {
}

}
}