#pragma once

#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Geometry {

class AABB {
public:
    AABB() {}
    AABB(const Vector3f& min, const Vector3f& max);
    AABB(const Vector3f& center, float32 sideSize);
    AABB(const AABB& other);

    const AABB& operator= (const AABB& other);

    const AABB& setCoords(const Vector3f& minCoord, const Vector3f& maxCoord);
    const AABB& fromCenterSize(const Vector3f& center, float32 size);
    const AABB& fromCenterDims(const Vector3f& center, const Vector3f& dims);

    const Vector3f& getMinimum() const { return m_minCoord; }
    const Vector3f& getMaximum() const { return m_maxCoord; }

    Vector3f getCenter() const;
    Vector3f getDimensions() const;

    Vector3f getClosestVertex(const Vector3f& point) const;

    inline bool onFrontHalf(Vector4f plane) const {
        Vector4f aabbVert;
        // get the vertex furthest from the plane
        aabbVert.x = (plane.x < 0.0f) ? m_minCoord.x : m_maxCoord.x;
        aabbVert.y = (plane.y < 0.0f) ? m_minCoord.y : m_maxCoord.y;
        aabbVert.z = (plane.z < 0.0f) ? m_minCoord.z : m_maxCoord.z;
        aabbVert.w = 1.0f;

        return plane.Dot(aabbVert) > 0.0f;
    }
    bool intersects(const AABB& other) const;
    bool sweptIntersects(const AABB& other, const Vector3f& dir, float32& t, Vector3f& normal) const;
    bool isInside(const Vector3f& point) const;
private:
    Vector3f m_minCoord;
    Vector3f m_maxCoord;
};

}
}