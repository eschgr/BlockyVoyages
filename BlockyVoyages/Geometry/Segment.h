#if !defined(__SEGMENT_H__)
#define __SEGMENT_H__

#include "../Math/MathMain.h"

namespace BlockyVoyages {
namespace Geometry {

class Segment {
public:
    Segment(const Vector2f& start, const Vector2f& end);

    const Vector2f& GetStart() const { return m_start; }
    const Vector2f& GetEnd() const { return m_end; }
    const Vector2f& GetDirection() const { return m_dir; }
    const Vector2f& GetNormal() const { return m_normal; }

    void SetStart(const Vector2f& start) { m_start = start; }
    void SetEnd(const Vector2f& end) { m_end = end; }
    // sets this to the opposite of the other segment. (start = other.end, normal = -other.normal, etc...)
    const Segment& MirrorOf(const Segment& other);
    const Segment& SetEnds(const Vector2f& start, const Vector2f& end);

    bool Intersects(const Segment& other, Vector2f& intersection) const;
    float32 DistanceTo(const Vector2f& point) const { return sqrt(DistanceToSquared(point)); }
    float32 DistanceToSquared(const Vector2f& point) const;

    float32 Length() const { return m_dir.Length(); }
    float32 LengthSquared() const { return m_dir.LengthSquared(); }

private:
    Vector2f m_start;
    Vector2f m_end;
    Vector2f m_dir;
    Vector2f m_normal;
};

}
}

#endif