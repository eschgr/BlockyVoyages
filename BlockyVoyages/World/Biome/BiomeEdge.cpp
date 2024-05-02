#include "stdafx.h"

#include "BiomeNode.h"
#include "BiomeEdge.h"
#include "BiomeRegion.h"

#include <assert.h>

namespace BlockyVoyages {
namespace World {
namespace Biome {

BiomeEdge::BiomeEdge(const Graph::Edge* other)
    : Edge(other)
{}

bool BiomeEdge::IsShore() const {
    // if either side is a boundary edge, it's assumed this edge is in the deep
    // ocean and can't be a shore edge.
    if (nullptr == m_right_cell || nullptr == m_left_cell) {
        return false;
    }

    BiomeRegion* cell0 = dynamic_cast<BiomeRegion*>(m_right_cell);
    BiomeRegion* cell1 = dynamic_cast<BiomeRegion*>(m_left_cell);

    assert (nullptr != cell0);
    assert (nullptr != cell1);

    if (cell0->IsWater() == cell1->IsWater()) {
        return false;
    }

    // No subbiomes in the ocean yet
    return cell0->IsWater() || cell1->IsWater();
}

bool BiomeEdge::IsLand() const {
    // if either side is a boundary edge, it's assumed this edge is in the deep
    // ocean and can't be a land edge.
    if (nullptr == m_right_cell || nullptr == m_left_cell) {
        return false;
    }

    BiomeRegion* cell0 = dynamic_cast<BiomeRegion*>(m_right_cell);
    BiomeRegion* cell1 = dynamic_cast<BiomeRegion*>(m_left_cell);

    assert (nullptr != cell0);
    assert (nullptr != cell1);

    return !(cell0->IsWater() || cell1->IsWater());
}

bool BiomeEdge::IsOcean() const {
    // if either side is a boundary edge, it's assumed this edge is in the deep
    // ocean and has to be an ocean.
    if (nullptr == m_right_cell || nullptr == m_left_cell) {
        return true;
    }

    BiomeRegion* cell0 = dynamic_cast<BiomeRegion*>(m_right_cell);
    BiomeRegion* cell1 = dynamic_cast<BiomeRegion*>(m_left_cell);

    assert (nullptr != cell0);
    assert (nullptr != cell1);
    
    return cell0->IsWater() && cell1->IsWater();
}

}
}
}