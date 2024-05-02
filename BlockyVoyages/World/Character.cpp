#include "stdafx.h"

#include "Character.h"

#include "Map.h"
#include "../Geometry/AABB.h"

namespace BlockyVoyages {
namespace World {

// these really should be elsewhere.
const float gravity = 9.81f; // m/s^2
const float terminalVel = 54.0f; // terminal velocity of a skydiver (m/s)
const float radius = 0.5f; // m
const float height = 1.75f; // m
const float stepHeigh = 0.3f; // not used
const float jump_speed = 5.5f; // m/s
const float walk_speed = 1.7f; // m/s
const float run_speed = 5.4f; // m/s

Character::Character(Map* worldMap)
    : m_controller(this),
      m_attachedCamera(nullptr),
      m_worldPos(0.0f, 0.0f, 0.0f),
      m_vel(0.0f, 0.0f, 0.0f),
      m_heading(0.0f),
      m_pitch(0.0f),
      m_stateFlags(0),
      m_worldMap(worldMap),
      m_dims(radius, height, radius)
{
    m_stats.walkSpeed = walk_speed;
    m_stats.runSpeed = run_speed;
    m_aabb.fromCenterDims(m_worldPos, m_dims);
}

Character::~Character() {}

void Character::setPosition(const Vector3f& pos) {
    m_worldPos = pos;
    m_aabb.fromCenterDims(m_worldPos, m_dims);
}

void Character::setRotation(float32 heading, float32 headPitch) {
    m_heading = heading;
    m_pitch = headPitch;
}

// set the direction of travel
void Character::setTravelDirection(const Vector3f& dir) {
    // if in the air, the travel direction of the character can't change
    
    // A character can only really travel in 2 dimensions.
    // Up and down don't make sense.
    Vector2f travelDir(dir.x, dir.z), tmp;
    travelDir.Normalize();
    Matrix22f rotMat;
    rotMat.Rotation(m_heading);
    tmp.Mul(rotMat, travelDir);
    // again, only transfer over the x and z components
    float32 speed = (m_stateFlags & RUNNING) ? m_stats.runSpeed : m_stats.walkSpeed;
    m_vel.x = tmp.x * speed;
    m_vel.z = tmp.y * speed;
}
    
// rotate the character
void Character::turnRightLeft(float32 angle) {
    m_heading += angle;
}

void Character::lookUpDown(float32 angle) {
    m_pitch += angle;
}

// toggle running
void Character::toggleRun() {
    m_stateFlags ^= RUNNING;
}

void Character::activateJump() {
    // if we're not already in the air and not already jumping, activate the jump
    if (!(m_stateFlags & IN_AIR) && !(m_stateFlags & JUMPING)) {
        m_vel.y = jump_speed;
        m_stateFlags |= JUMPING;
    }
}

void Character::toggleCrouch() {
    m_stateFlags ^= CROUCHING;
}

// use right hand
void Character::useRightHand() {
}

// use left hand
void Character::useLeftHand() {
}

// set active item. Doesn't know what item would be switched to.
void Character::setActiveItem(int32 slot) {
}

void Character::attachCamera(Input::CameraController* camController) {
    m_attachedCamera = camController;

    if (nullptr == camController) {
        m_vel.Zero();
    }
}

void Character::update(float32 dt) {
    // if the character is in the air, then add acceleration due to gravity to
    // the velocity, otherwise, just add it to the minimum coordinate to check
    // if the character is still on ground.
    if (m_stateFlags & IN_AIR) {
        m_vel.y -= gravity * dt;
        if (m_vel.y < -terminalVel) {
            m_vel.y = -terminalVel;
        }
    }

    // calculate the large AABB around the movement area
    Vector3f minCoord(m_aabb.getMinimum());
    Vector3f maxCoord(m_aabb.getMaximum());

    if (m_vel.x > 0) {
        maxCoord.x += m_vel.x * dt;
    } else {
        minCoord.x += m_vel.x * dt;
    }

    if (m_vel.y > 0) {
        maxCoord.y += m_vel.y * dt;
    } else {
        minCoord.y += m_vel.y * dt;
    }

    if (m_vel.z > 0) {
        maxCoord.z += m_vel.z * dt;
    } else {
        minCoord.z += m_vel.z * dt;
    }

    // get those blocks which can be used to walk up automatically.
    // TODO(gesch): allow the character to move up a single block automatically
    /*if (!(m_stateFlags & IN_AIR)) {
        maxCoord.y += 0.3f;
    }*/
    std::vector<MapNode*> intersectingNodes;
    Geometry::AABB bigAABB(minCoord, maxCoord);
    m_worldMap->getIntersectingLeaves(bigAABB, intersectingNodes);
    if (intersectingNodes.size() > 0) {
        Vector3f minNormal;
        float32 minT;
        do {
            minT = 1.0f;
            for (auto* node : intersectingNodes) {
                float32 t;
                Vector3f normal;
                if (node->m_aabb.sweptIntersects(m_aabb, m_vel, t, normal)) {
                    if (t < minT) {
                        minT = t;
                        minNormal = normal;
                    }
                }
            }
            if (minT < 1.0f) {
                float32 dotProd = m_vel.Dot(minNormal) * (1.0f - minT);
                m_vel -= dotProd * minNormal;
            }
        } while (minT < 1.0f);
    } else {
        // Since there should always be some blocks in the list if the
        // character is on the ground, mark that the character is in the air.
        m_stateFlags |= IN_AIR;
    }

    // move the position of the player
    m_worldPos += m_vel * dt;
    m_aabb.fromCenterDims(m_worldPos, m_dims);
    if (m_attachedCamera) {
        Vector3f headPos(m_worldPos);
        headPos.y += 0.8f;
        m_attachedCamera->setPosition(headPos);
        m_attachedCamera->setRotation(m_heading, m_pitch);
    }

    intersectingNodes.clear();
    bigAABB.fromCenterDims(m_worldPos - Vector3f(0.0f, 0.025f, 0.0f), m_dims);
    m_worldMap->getIntersectingLeaves(bigAABB, intersectingNodes);
    // check to see if the player is on the ground.
    if (0 == intersectingNodes.size()) {
        m_stateFlags |= IN_AIR;
    } else {
        m_stateFlags &= ~JUMPING;
    }
}

void Character::draw() {
    return;
    Vector3f dims(0.5f, 1.75f, 0.2f);
    dims *= 0.5f;

    glPushMatrix();
    glTranslatef(m_worldPos.x, m_worldPos.y, m_worldPos.z);
    glRotatef(RadToDeg(m_heading), 0.0f, 1.0f, 0.0f);

    Vector3f lowBounds(-dims), highBounds(dims);
    
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINE_STRIP);
    glVertex3f(lowBounds.x, lowBounds.y, lowBounds.z);
    glVertex3f(highBounds.x, lowBounds.y, lowBounds.z);
    glVertex3f(highBounds.x, lowBounds.y, highBounds.z);
    glVertex3f(lowBounds.x, lowBounds.y, highBounds.z);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(highBounds.x, lowBounds.y, lowBounds.z);
    glVertex3f(highBounds.x, highBounds.y, lowBounds.z);
    glVertex3f(highBounds.x, highBounds.y, highBounds.z);
    glVertex3f(highBounds.x, lowBounds.y, highBounds.z);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(highBounds.x, highBounds.y, highBounds.z);
    glVertex3f(lowBounds.x, highBounds.y, highBounds.z);
    glVertex3f(lowBounds.x, highBounds.y, lowBounds.z);
    glVertex3f(highBounds.x, highBounds.y, lowBounds.z);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(lowBounds.x, highBounds.y, highBounds.z);
    glVertex3f(lowBounds.x, lowBounds.y, highBounds.z);
    glVertex3f(lowBounds.x, lowBounds.y, lowBounds.z);
    glVertex3f(lowBounds.x, highBounds.y, lowBounds.z);
    glEnd();

    glPopMatrix();
}

}
}