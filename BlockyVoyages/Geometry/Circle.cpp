#include "stdafx.h"

#include "Circle.h"

namespace BlockyVoyages {
namespace Geometry {

Circle::Circle(const Vector2f& center, float32 radius)
    : m_center(center),
      m_radius(radius)
{}

// calculate the center and radius of a circle given three points on the circle
Circle::Circle(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2) {

    if (!IsPerpendicular(p0, p1, p2)) {
        CalculateCircle(p0, p1, p2);
    } else if (!IsPerpendicular(p0, p2, p1)) {
        CalculateCircle(p0, p2, p1);
    } else if (!IsPerpendicular(p1, p0, p2)) {
        CalculateCircle(p1, p0, p2);
    } else if (!IsPerpendicular(p1, p2, p0)) {
        CalculateCircle(p1, p2, p0);
    } else if (!IsPerpendicular(p2, p1, p0)) {
        CalculateCircle(p2, p1, p0);
    } else if (!IsPerpendicular(p2, p0, p1)) {
        CalculateCircle(p2, p0, p1);
    } else {
        //LOG(Debugability::DBG_LVL_LOW, "The three pts are perpendicular to axis");
        m_radius = -1.0f;
    }

    //LOG(Debugability::DBG_LVL_LOW, "Center: %f %f\n", m_center.x, m_center.y);
    //LOG(Debugability::DBG_LVL_LOW, "Radius: %f\n", m_radius);
}

// Check the given point are perpendicular to x or y axis 
bool Circle::IsPerpendicular(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2) const {
    Vector2f delta_a;
    Vector2f delta_b;

    delta_a.Sub(p1, p0);
    delta_b.Sub(p2, p1);
    
    //LOG(Debugability::DBG_LVL_LOW, " delta_a: (%f, %f)", delta_a.x, delta_a.y);
    //LOG(Debugability::DBG_LVL_LOW, " delta_b: (%f, %f)", delta_b.x, delta_b.y);

    // checking whether the line of the two pts are vertical
    if (IsCloseToZero(delta_a.x) && IsCloseToZero(delta_b.y)) {
        //LOG(Debugability::DBG_LVL_LOW, "The points are pependicular and parallel to x-y axis");
        return false;
    }

    if (IsCloseToZero(delta_a.y)) {
        //LOG(Debugability::DBG_LVL_LOW, "A line of two point are perpendicular to x-axis 1");
        return true;
    }
    else if (IsCloseToZero(delta_b.y)) {
        //LOG(Debugability::DBG_LVL_LOW, "A line of two point are perpendicular to x-axis 2");
        return true;
    }
    else if (IsCloseToZero(delta_a.x)) {
        //LOG(Debugability::DBG_LVL_LOW, "A line of two point are perpendicular to y-axis 1");
        return true;
    }
    else if (IsCloseToZero(delta_b.x)) {
        //LOG(Debugability::DBG_LVL_LOW, "A line of two point are perpendicular to y-axis 2");
        return true;
    }
    else return false ;
}

bool Circle::CalculateCircle(const Vector2f& p0, const Vector2f& p1, const Vector2f& p2) {
    Vector2f l10, l20;
    l10.Sub(p1, p0);
    l20.Sub(p2, p0);

    float32 area = (l10.x * l20.y - l10.y * l20.x);

    // ignore all points which makes a left turn or where the determinate is 0.
    if (IsCloseToZero(area)) {
      return false;
    }

    m_center.Set((l20.y * l10.LengthSquared() - l10.y * l20.LengthSquared()) * 0.5f / area,
                 (l10.x * l20.LengthSquared() - l20.x * l10.LengthSquared()) * 0.5f / area);
    double radius = m_center.Length();
    m_center += p0;

    return true;
}

}
}