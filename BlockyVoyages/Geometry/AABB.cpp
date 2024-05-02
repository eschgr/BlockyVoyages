#include "stdafx.h"

#include "AABB.h"

namespace BlockyVoyages {
namespace Geometry {

AABB::AABB(const Vector3f& min, const Vector3f& max)
    : m_minCoord(min),
      m_maxCoord(max) 
{}

AABB::AABB(const Vector3f& point, float32 sideSize) {
    fromCenterSize(point, sideSize);
}

AABB::AABB(const AABB& other) {
    if (&other == this) {
        return;
    }
    m_minCoord = other.m_minCoord;
    m_maxCoord = other.m_maxCoord;
}

const AABB& AABB::operator= (const AABB& other) {
    if (&other == this) {
        return *this;
    }
    m_minCoord = other.m_minCoord;
    m_maxCoord = other.m_maxCoord;
    return *this;
}

const AABB& AABB::setCoords(const Vector3f& minCoord, const Vector3f& maxCoord) {
    m_minCoord = minCoord;
    m_maxCoord = maxCoord;

    return *this;
}

const AABB& AABB::fromCenterSize(const Vector3f& point, float32 size) {
    size *= 0.5f;
    m_minCoord.x = point.x - size;
    m_minCoord.y = point.y - size;
    m_minCoord.z = point.z - size;

    m_maxCoord.x = point.x + size;
    m_maxCoord.y = point.y + size;
    m_maxCoord.z = point.z + size;

    return *this;
}

const AABB& AABB::fromCenterDims(const Vector3f& point, const Vector3f& dims) {
    m_maxCoord.Mul(dims, 0.5f);
    m_minCoord.Sub(point, m_maxCoord);
    m_maxCoord.Add(point, m_maxCoord);
    return *this;
}

Vector3f AABB::getCenter() const {
    Vector3f point;
    point.Add(m_minCoord, m_maxCoord);
    point *= 0.5f;
    return point;
}

Vector3f AABB::getDimensions() const {
    Vector3f dims;
    dims.Sub(m_maxCoord, m_minCoord);
    return dims;
}

Vector3f AABB::getClosestVertex(const Vector3f& point) const {
    Vector3f diff;
    diff.Sub(point, getCenter());

    Vector3f farPoint;
    farPoint.x = diff.x > 0.0f ? m_maxCoord.x : m_minCoord.x;
    farPoint.y = diff.y > 0.0f ? m_maxCoord.y : m_minCoord.y;
    farPoint.z = diff.z > 0.0f ? m_maxCoord.z : m_minCoord.z;

    return farPoint;
}

bool AABB::intersects(const AABB& other) const {
    if (m_minCoord.x > other.m_maxCoord.x) {
        return false;
    }
    if (other.m_minCoord.x > m_maxCoord.x) {
        return false;
    }
    if (m_minCoord.y > other.m_maxCoord.y) {
        return false;
    }
    if (other.m_minCoord.y > m_maxCoord.y) {
        return false;
    }
    if (m_minCoord.z > other.m_maxCoord.z) {
        return false;
    }
    if (other.m_minCoord.z > m_maxCoord.z) {
        return false;
    }
    return true;
}

bool AABB::sweptIntersects(const AABB& other, const Vector3f& dir, float32& t, Vector3f& normal) const {
    Vector3f invEntry, invExit;

    // find the distance between the objects on the near and far sides for both x and y
    if (dir.x > 0.0f) {
        invEntry.x = m_minCoord.x - other.m_maxCoord.x;
        invExit.x = m_maxCoord.x - other.m_minCoord.x;
    } else {
        invEntry.x = m_maxCoord.x - other.m_minCoord.x;
        invExit.x = m_minCoord.x - other.m_maxCoord.x;
    }

    if (dir.y > 0.0f) {
        invEntry.y = m_minCoord.y - other.m_maxCoord.y;
        invExit.y = m_maxCoord.y - other.m_minCoord.y;
    } else {
        invEntry.y = m_maxCoord.y - other.m_minCoord.y;
        invExit.y = m_minCoord.y - other.m_maxCoord.y;
    }

    if (dir.z > 0.0f) {
        invEntry.z = m_minCoord.z - other.m_maxCoord.z;
        invExit.z = m_maxCoord.z - other.m_minCoord.z;
    } else {
        invEntry.z = m_maxCoord.z - other.m_minCoord.z;
        invExit.z = m_minCoord.z - other.m_maxCoord.z;
    }

    Vector3f entry, exit;

    if (dir.x == 0) {
        entry.x = -std::numeric_limits<float>::infinity();
        exit.x = std::numeric_limits<float>::infinity();
    } else {
        entry.x = invEntry.x / dir.x;
        exit.x = invExit.x / dir.x;
    }

    if (IsZero(dir.y)) {
        entry.y = -std::numeric_limits<float>::infinity();
        exit.y = std::numeric_limits<float>::infinity();
    } else {
        entry.y = invEntry.y / dir.y;
        exit.y = invExit.y / dir.y;
    }

    if (dir.z == 0) {
        entry.z = -std::numeric_limits<float>::infinity();
        exit.z = std::numeric_limits<float>::infinity();
    } else {
        entry.z = invEntry.z / dir.z;
        exit.z = invExit.z / dir.z;
    }

    // find the earliest/latest times of collision
    float entryTime, exitTime;
    if (entry.x > entry.y) {
        if (entry.x > entry.z) {
            entryTime = entry.x;
            normal.y = 0.0f;
            normal.z = 0.0f;
            if (invEntry.x < 0.0f) {
                normal.x = 1.0f;
            } else {
                normal.x = -1.0f;
            }
        } else {
            entryTime = entry.z;
            normal.x = 0.0f;
            normal.y = 0.0f;
            if (invEntry.z < 0.0f) {
                normal.z = 1.0f;
            } else {
                normal.z = -1.0f;
            }
        }
    } else if (entry.y > entry.z) {
        entryTime = entry.y;
        normal.x = 0.0f;
        normal.z = 0.0f;
        if (invEntry.y < 0.0f) {
            normal.y = 1.0f;
        } else {
            normal.y = -1.0f;
        }
    } else {
        entryTime = entry.z;
        normal.x = 0.0f;
        normal.y = 0.0f;
        if (invEntry.z < 0.0f) {
            normal.z = 1.0f;
        } else {
            normal.z = -1.0f;
        }
    }

    if (exit.x < exit.y) {
        if (exit.x < exit.z) {
            exitTime = exit.x;
        } else {
            exitTime = exit.z;
        }
    } else if (exit.y < exit.z) {
        exitTime = exit.y;
    } else {
        exitTime = exit.z;
    }

    // if there was no collision
    if ((entryTime > exitTime) || 
        (entry.x < 0.0f && entry.y < 0.0f && entry.z < 0.0f) || 
        (entry.x >= 1.0f || entry.y >= 1.0f || entry.z >= 1.0f)) {
        t = 1.0f;
        normal.Zero();
        return false;
    }

    // return the time of collision
    t = entryTime;
    return true;
}

bool AABB::isInside(const Vector3f& point) const {
    return (point.x > m_minCoord.x && point.x <= m_maxCoord.x &&
            point.y > m_minCoord.y && point.y <= m_maxCoord.y &&
            point.z > m_minCoord.z && point.z <= m_maxCoord.z);
}

}
}