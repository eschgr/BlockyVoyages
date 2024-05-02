#if !defined(__CIRCLE_H__)
#define __CIRCLE_H__

#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Geometry {

class Circle {
public:
    // basic constructor.
    Circle(const Vector2f& center, float32 radius);

    // calculate the center and radius of a circle given three points on the circle
    Circle(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2);

    const Vector2f& GetCenter() const { return m_center; }
    float32 GetRadius() const { return m_radius; }

private:
    Vector2f m_center;
    float32 m_radius;

    bool IsPerpendicular(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2) const;
    bool CalculateCircle(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2);
};

}
}

#endif