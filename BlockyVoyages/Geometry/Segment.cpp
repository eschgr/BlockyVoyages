#include "stdafx.h"

#include "Segment.h"

#include <assert.h>

namespace BlockyVoyages {
namespace Geometry {

Segment::Segment(const Vector2f& start, const Vector2f& end)
    : m_start(start),
      m_end(end)
{
    m_dir.Sub(end, start);
    m_normal.PerpendicularOf(m_dir);
}

// sets this to the opposite of the other segment
const Segment& Segment::MirrorOf(const Segment& other) {
    m_start = other.m_end;
    m_end = other.m_start;
    m_dir = -other.m_dir;
    m_normal = -other.m_normal;
    return *this;
}

const Segment& Segment::SetEnds(const Vector2f& start, const Vector2f& end) {
    m_start = start;
    m_end = end;
    m_dir.Sub(end, start);
    m_normal.PerpendicularOf(m_dir);
    return *this;
}

bool Segment::Intersects(const Segment& other, Vector2f& intersection) const {
    Vector2f w0 = m_start - other.m_start;
    float32 a = m_dir.Dot(m_dir);
    float32 b = m_dir.Dot(other.m_dir);
    float32 c = other.m_dir.Dot(other.m_dir);
    float32 d = m_dir.Dot(w0);
    float32 e = other.m_dir.Dot(w0);
    float32 sn = b * e - c * d;
    float32 tn = a * e - b * d;

    // check to ensure the intersection is on the segments
    if (sn < 0.0f) {
        return false;
    }
    if (tn < 0.0f) {
        return false;
    }

    float32 denom = a * c - b * b;
    if (sn > denom) {
        return false;
    }
    if (tn > denom) {
        return false;
    }

    intersection.Mul(m_dir, sn / denom);
    intersection += m_start;
    
    return true;
}

float32 Segment::DistanceToSquared(const Vector2f& point) const {
    Vector2f w = point - m_start;
    float32 proj = w.Dot(m_dir);
    if (proj <= 0.0f) {
        return w.Length();
    } else {
        float32 vsq = m_dir.LengthSquared();
        if (proj >= vsq) {
            return w.LengthSquared() - 2.0f * proj + vsq;
        } else {
            return w.LengthSquared() - proj * proj / vsq;
        }
    }
    return 0.0f;
}

}
}