#include "stdafx.h"
#include "BiomeNode.h"

#include "BiomeEdge.h"

#include <assert.h>

namespace BlockyVoyages {
namespace World {
namespace Biome {

const float32 kMaxMoisture = 10.0f;
const float32 kMinMoisture = 0.0f;
const float32 kLandFalloff = 0.02f; // the amount of moisture removed by land / km
const float32 kWaterAdd = 0.075f; // the amount of moisture added by water / km
const float32 kHeightFalloff = 0.125f; // the amount of moisture removed by height / km

bool BiomeNode::CalculateMoisture(const Vector2f& wind_dir) {
    float32 new_moisture = 0.0f;
    for (auto* edge : m_edges) {
        Geometry::Segment edge_seg = edge->GetSegment(this);
        float32 edge_in_dir = edge_seg.GetDirection().Dot(wind_dir);
        if (edge_in_dir < 0.0f) {
            BiomeEdge* bedge = dynamic_cast<BiomeEdge*>(edge);
            BiomeNode* other_node = dynamic_cast<BiomeNode*>(bedge->GetOtherNode(this));
            float32 moisture = other_node->m_moisture;
            if (bedge->IsLand()) {
                moisture -= kLandFalloff * edge_seg.Length();
                if (m_height > other_node->m_height) {
                    moisture -= kHeightFalloff * (m_height - other_node->m_height);
                }
            } else if (bedge->IsOcean()) {
                moisture += kWaterAdd * edge_seg.Length();
            }
            if (moisture > new_moisture) {
                new_moisture = moisture;
            }
        }
    }
    new_moisture = max(kMinMoisture, min(kMaxMoisture, new_moisture));
    if (new_moisture > m_moisture) {
        m_moisture = new_moisture;
        return true;
    }
    return false;
}

bool BiomeNode::IsShore() const {
    // if any of the edges are a shore edge, then this is a shore node
    for (auto* edge : m_edges) {
        BiomeEdge* bedge = dynamic_cast<BiomeEdge*>(edge);
        assert (nullptr != bedge);
        if (bedge->IsShore()) {
            return true;
        }
    }
    return false;
}

bool BiomeNode::IsLand() const {
    // all the edges must be a land edge for this to be land
    for (auto* edge : m_edges) {
        BiomeEdge* bedge = dynamic_cast<BiomeEdge*>(edge);
        assert (nullptr != bedge);
        if (!bedge->IsLand()) {
            return false;
        }
    }
    return true;
}

bool BiomeNode::IsOcean() const {
    // all the edges must be a ocean edge for this to be ocean
    for (auto* edge : m_edges) {
        BiomeEdge* bedge = dynamic_cast<BiomeEdge*>(edge);
        assert (nullptr != bedge);
        if (!bedge->IsOcean()) {
            return false;
        }
    }
    return true;
}

}
}
}