#if !defined(__BIOME_EDGE_H__)
#define __BIOME_EDGE_H__

#include "../../Geometry/Segment.h"
#include "../../Graph/Edge.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class BiomeEdge : public Graph::Edge {
public:
    BiomeEdge(const Graph::Edge* other);

    bool IsShore() const;
    bool IsLand() const;
    bool IsOcean() const;

private:
};

}
}
}

#endif