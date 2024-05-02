#include "stdafx.h"

#include "WorldGenerator.h"
#include "Biome/Biome.h"
#include "Biome/BiomeFactory.h"
#include "Biome/BiomeMap.h"
#include "Biome/BiomeRegion.h"
#include "../Graph/Voronoi/Voronoi.h"
#include "../Math/SimplexNoise.h"

namespace BlockyVoyages {
namespace World {
using Math::SimplexNoise;

WorldGenerator::WorldGenerator(float32 feature_size)
    : m_seed(0),
      m_feature_size(feature_size),
      m_half_feature_size(feature_size * 0.5f),
      m_max_height(0.0f),
      m_area(0.0f),
      m_min_bounds(0.0f, 0.0f, 0.0f),
      m_max_bounds(0.0f, 0.0f, 0.0f),
      m_cur_region(nullptr)
{}

WorldGenerator::~WorldGenerator() {
}

Vector2f WorldGenerator::GetWorldDimensions() {
    if (m_biome_map.get()) {
        Vector2f bounds(m_max_bounds.x - m_min_bounds.x,
                        m_max_bounds.z - m_min_bounds.z);
        return bounds;
    }
    return Vector2f(0.0f, 0.0f);
}

void WorldGenerator::GenerateBiomeMap(int32 seed, float32 max_height, float32 area) {
    m_seed = seed;
    m_area = area;
    m_max_height = std::ceil(max_height / m_feature_size) * m_feature_size;

    World::Biome::BiomeFactory biome_factory;
    biome_factory.SetMoistureBounds(Vector2f(0.0f, 10.0f));
    biome_factory.SetTemperatureBounds(Vector2f(-20.0f, 60.0f));

    GenerateBounds();
    m_biome_map.reset(new World::Biome::BiomeMap(m_seed, max_height));
    m_biome_map->CopyFrom(GenerateBaseGraph().get());

    // create the map
    m_biome_map->AssignWater();
    m_biome_map->AssignHeight();
    m_biome_map->AssignTemperature();
    m_biome_map->AssignMoisture();
    m_biome_map->AssignBiomes(biome_factory);
}

// This starts the process of generating a large number of base blocks. To
// generate a region of blocks, start by calling InitializeRegion. Next, 
// call GetNextBlock in a loop until false is returned. Once done, call
// ReleaseRegion to clean up the working space.
bool WorldGenerator::InitializeRegion(const Vector2f& minCoord, const Vector2f& maxCoord) {
    m_min_bounds.x = std::floor(minCoord.x / m_feature_size) * m_feature_size + m_half_feature_size;
    m_min_bounds.z = std::floor(minCoord.y / m_feature_size) * m_feature_size + m_half_feature_size;

    m_max_bounds.x = std::ceil(maxCoord.x / m_feature_size) * m_feature_size - m_half_feature_size;
    m_max_bounds.z = std::ceil(maxCoord.y / m_feature_size) * m_feature_size - m_half_feature_size;

    m_int_dims.x = static_cast<int32>((m_max_bounds.x - m_min_bounds.x) / m_feature_size) + 1;
    m_int_dims.y = static_cast<int32>((m_max_bounds.z - m_min_bounds.z) / m_feature_size) + 1;
    m_buffer_dims.Set(m_int_dims.x + 2, m_int_dims.y + 2);

    m_cur_pos.x = 0;
    m_cur_pos.y = 0;

    m_height_cache.resize(m_buffer_dims.x * m_buffer_dims.y);
    int32 x, y, ind;
    Vector2f world_coord(m_min_bounds.x - m_feature_size, 0.0f);
    for (x = 0, ind = 0; x < m_buffer_dims.x; ++x, world_coord.x += m_feature_size) {
        world_coord.y = m_min_bounds.z - m_feature_size;
        for (y = 0; y < m_buffer_dims.y; ++y, ++ind, world_coord.y += m_feature_size) {
            m_height_cache[ind] = GetSurfaceHeight(world_coord);
        }
    }

    CulculateColumnStats();
    return true;
}

bool WorldGenerator::GetNextBlock(Vector3f& center, BlockType& type) {
    if (m_cur_pos.x >= m_int_dims.x) {
        return false;
    }

    center.Set(m_cur_pos.x * m_feature_size + m_min_bounds.x,
               m_cur_height,
               m_cur_pos.y * m_feature_size + m_min_bounds.z);

    {
        Vector2f coord(center.x, center.z);
        if (!m_cur_region || !m_cur_region->PointInside(coord)) {
            m_cur_region = dynamic_cast<Biome::BiomeRegion*>(m_biome_map->GetCellAt(coord));
        }
        if (nullptr == m_cur_region) {
            type = kDirt;
        } else {
            type = m_cur_region->GetBiome()->GetSurfaceType(coord);
        }
    }

    m_cur_height += m_feature_size;
    if (m_cur_height > m_max_bounds.y) {
        ++m_cur_pos.y;
        if (m_cur_pos.y >= m_int_dims.y) {
            ++m_cur_pos.x;
            m_cur_pos.y = 0;
        }
        CulculateColumnStats();
    }
    return true;
}

void WorldGenerator::ReleaseRegion() {
    std::vector<float32>().swap(m_height_cache);
}

// GetBlockType is meant to be used to determine the type of a leaf node.
// It is not meant to be used to generate a large area of terrain.
BlockType WorldGenerator::GetBlockType(const Vector3f& blockCenter) {
    return kAir;
}

// GetBlockLODType returns the type of a large block. It can find the type
// of a block that covers a large area. Not meant for leaf nodes or near
// leaf nodes. For near leaf nodes, it is more efficient and just as good
// to take the type that occurs the most of the child blocks.
BlockType WorldGenerator::GetLODBlockType(const Vector3f& blockCenter, float32 sideSize) {
    return kAir;
}

// Generates the bounds for the world.
void WorldGenerator::GenerateBounds() {
    std::mt19937 rng(m_seed);
    std::uniform_real_distribution<float32> aspect_dist(0.5, 1.5);
    float32 aspect = aspect_dist(rng);
    float32 half_sqrt_area = sqrt(m_area) * 0.5f / m_feature_size;
    float32 half_width = ceil(half_sqrt_area * aspect) * m_feature_size;
    float32 half_height = ceil(half_sqrt_area / aspect) * m_feature_size;

    m_max_bounds.Set(half_width, m_max_height, half_height);
    m_min_bounds.Set(-half_width, 0.0f, -half_height);
}

std::unique_ptr<Graph::Graph> WorldGenerator::GenerateBaseGraph() {
    static const int32 kNumSmoothingIters = 2;

    std::mt19937 rng(m_seed);
    // the distribution doesn't go to the edge to allow some space between sites and the boundary
    std::uniform_real_distribution<float32> width_dist(m_min_bounds.x, m_max_bounds.x);
    std::uniform_real_distribution<float32> height_dist(m_min_bounds.z, m_max_bounds.z);

    // Create the initial bunch of random points.
    std::vector<Vector2f> verts;
    verts.reserve(kNumRegions);
    while (verts.size() < static_cast<uint32>(kNumRegions)) {
        verts.push_back(Vector2f(width_dist(rng), height_dist(rng)));
    }

    // The Voronoi diagram calculation takes 2D bounds instead of 3D bounds.
    // Use x, z instead of x, y since the world lays in the x-z plane 
    Vector2f min_bounds(m_min_bounds.x, m_min_bounds.z);
    Vector2f max_bounds(m_max_bounds.x, m_max_bounds.z);
    
    // Create the initial Voronoi diagram.
    std::unique_ptr<Graph::Voronoi> out(new Graph::Voronoi(verts, min_bounds, max_bounds));
    for (int32 iter = 0; iter < kNumSmoothingIters; ++iter) {
        // get the vertices from the diagram
        // resize the vertices since it's possible some were removed
        const std::vector<Graph::Cell*>& cells = out->GetCells();
        verts.resize(cells.size());
        for (uint32 vert = 0; vert < verts.size(); ++vert) {
            verts[vert] = cells[vert]->GetCenter();
        }

        // recalculate the diagram
        out.reset(new Graph::Voronoi(verts, min_bounds, max_bounds));
    }
    return std::unique_ptr<Graph::Graph>(out.release());
}

float32 WorldGenerator::GetSurfaceHeight(const Vector2f& coord) {
    // Get the cell the coordinate is over
    if (!m_cur_region || !m_cur_region->PointInside(coord)) {
        m_cur_region = dynamic_cast<Biome::BiomeRegion*>(m_biome_map->GetCellAt(coord));
    }
    if (nullptr == m_cur_region) {
        return 1000.0f;
    }
    return ceil(m_cur_region->GetHeightAt(coord) / m_feature_size) * m_feature_size;
    //float32 offset = SimplexNoise::noise2DOctaves(coord.x, coord.y, 200.0f, 16);
    //return static_cast<int32>(((offset + 1.0f) * m_max_height * 0.5f / m_feature_size)) * m_feature_size;
}

float32 WorldGenerator::GetCachedHeight(int32 x, int32 y) {
    // -1 and the maximum dimension is acceptable because of the buffer
    if (m_height_cache.empty() ||
        x < -1 || x > m_int_dims.x ||
        y < -1 || y > m_int_dims.y) {
        return 0.0f;
    }
    // change the coordinates so that they take the buffer into account
    ++x;
    ++y;

    return m_height_cache[x * m_buffer_dims.y + y];
}

void WorldGenerator::CulculateColumnStats() {
    m_max_bounds.y = GetCachedHeight(m_cur_pos.x, m_cur_pos.y);
    m_min_bounds.y = m_max_bounds.y;

    float32 height = GetCachedHeight(m_cur_pos.x - 1, m_cur_pos.y);
    if (height < m_min_bounds.y) {
        m_min_bounds.y = height;
    }
    height = GetCachedHeight(m_cur_pos.x + 1, m_cur_pos.y);
    if (height < m_min_bounds.y) {
        m_min_bounds.y = height;
    }
    height = GetCachedHeight(m_cur_pos.x, m_cur_pos.y - 1);
    if (height < m_min_bounds.y) {
        m_min_bounds.y = height;
    }
    height = GetCachedHeight(m_cur_pos.x, m_cur_pos.y + 1);
    if (height < m_min_bounds.y) {
        m_min_bounds.y = height;
    }
    m_cur_height = m_min_bounds.y;
    if (m_min_bounds.y != m_max_bounds.y) {
        // there's no need to start at the same level as the neighbor as the
        // neighboring block will fill in the spot.
        m_cur_height += m_feature_size;
    }
}

} // World
} // BlockyVoyages