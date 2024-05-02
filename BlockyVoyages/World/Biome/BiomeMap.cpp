#include "stdafx.h"

#include "BiomeMap.h"

#include "BiomeNode.h"
#include "BiomeEdge.h"
#include "BiomeFactory.h"
#include "BiomeRegion.h"
#include "../../Graphics/OpenGL.h"
#include "../../Graphics/Shader.h"
#include "../../Graphics/VertexArrayObject.h"
#include "../../Math/SimplexNoise.h"
#include "../../Math/SparseMatrix.h"

#include <set>
#include <assert.h>

namespace BlockyVoyages {
namespace World {
namespace Biome {

BiomeMap::BiomeMap(int32 seed, float32 max_height)
    : m_rng(seed),
      m_max_height(max_height),
      m_is_south_hemisphere(false)
{
    std::uniform_real_distribution<float32> dist(0.0f, 100000.0f);
    m_simplex_offset.Set(dist(m_rng), dist(m_rng));
    m_simplex_scale = 125.0f;
}

void BiomeMap::CopyFrom(const Graph* other) {
    Graph::CopyFrom(other);
    BuildRegionGrid();
}

Cell* BiomeMap::GetCellAt(const Vector2f& pnt) const {
    int32 row = static_cast<int32>(ceil(pnt.y - m_min_bounds.y)) / kCellSize;
    int32 col = static_cast<int32>(ceil(pnt.x - m_min_bounds.x)) / kCellSize;

    if (row >= m_region_grid_height ||
        col >= m_region_grid_width) {
        return nullptr;
    }

    for (auto* cell : m_region_grid[row * m_region_grid_width + col]) {
        if (cell->PointInside(pnt)) {
            return cell;
        }
    }
    return nullptr;
}

void BiomeMap::AssignWater() {
    std::set<BiomeRegion*> water;
    std::queue<BiomeRegion*> ocean;
    for (auto* cell : m_cells) {
        BiomeRegion* region = dynamic_cast<BiomeRegion*>(cell);
        assert(region != nullptr);
        // if this cell is a boundary cell, just add it to the list
        if (region->AtBoundary()) {
            // add it to the water list and save it as a starting cell.
            ocean.push(region);
            region->SetIsWater();
        } else {
            Vector2f cell_center = cell->GetCenter();
            float32 noise = (Math::SimplexNoise::noise2DOctaves(cell_center.x + m_simplex_offset.x,
                                                                cell_center.y + m_simplex_offset.y,
                                                                m_simplex_scale, 4) + 1.0f) * 0.5f;
            Vector2f normalized_center;
            // set the min and max so the majority of the world will be land
            normalized_center.Sub(cell_center, m_min_bounds);
            normalized_center.x /= (m_max_bounds.x - m_min_bounds.x) * 0.5f;
            normalized_center.y /= (m_max_bounds.y - m_min_bounds.y) * 0.5f;
            normalized_center.x -= 1.0f;
            normalized_center.y -= 1.0f;

            float32 dist_to_center;
            if (fabs(normalized_center.x) > fabs(normalized_center.y)) {
                // x is the closer side
                dist_to_center = fabs(normalized_center.x);
            } else {
                dist_to_center = fabs(normalized_center.y);
            }
            float32 min_val = 0.2f;
            float32 max_val = 0.9f;
            if (noise > min_val + max_val * (1.0f - dist_to_center)) {
                water.insert(region);
            }
        }
    }

    // use flood fill to create the oceans
    while (!ocean.empty()) {
        BiomeRegion* region = ocean.front();
        ocean.pop();
        for (auto* edge : region->GetEdges()) {
            BiomeRegion* other = dynamic_cast<BiomeRegion*>(edge->GetOtherCell(region));
            // if other is null, the edge is a boundary edge and we can skip it.
            if (other == nullptr) {
                continue;
            }
            auto other_water = water.find(other);
            if (other_water != water.end()) {
                ocean.push(other);
                water.erase(other_water);
                other->SetIsWater();
            }
        }
    }
    
    CalculateDistanceToShore();
}

void BiomeMap::AssignHeight() {
    // first, get a random number of ridges
    std::uniform_int_distribution<int32> ridge_dist(2, 4);
    Vector2f border = (m_max_bounds - m_min_bounds) * 0.15f;
    float32 min_ridge_length = border.Length();
    std::uniform_real_distribution<float32> end_x_dist(m_min_bounds.x + border.x, m_max_bounds.x - border.x);
    std::uniform_real_distribution<float32> end_y_dist(m_min_bounds.y + border.y, m_max_bounds.y - border.y);
    int32 n_ridges = ridge_dist(m_rng);
    m_mountains.clear();
    m_mountains.reserve(n_ridges);
    for (int32 ridge = 0; ridge < n_ridges; ++ridge) {
        std::map<BiomeNode*, float32> set_heights;
        Vector2f start_pnt(end_x_dist(m_rng), end_y_dist(m_rng));
        Vector2f end_pnt(end_x_dist(m_rng), end_y_dist(m_rng));
        while (start_pnt.DistanceSquared(end_pnt) < border * border) {
            end_pnt.Set(end_x_dist(m_rng), end_y_dist(m_rng));
        }
        m_mountains.push_back(Geometry::Segment(start_pnt, end_pnt));

        // Set those nodes which are shore nodes
        for (Node* node : m_nodes) {
            BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
            if (bnode->GetShoreDistance() == 0.0f) {
                set_heights[bnode] = 0.0f;
            }
        }

        // find which cell the start point is in
        /*start_pnt.Add(m_min_bounds, m_max_bounds);
        start_pnt *= 0.5f;*/
        Cell* cur_cell = GetCellAt(start_pnt);
        // it's expected that this ridge be in the world.
        assert(cur_cell);
        // Set the heights for the ridge.
        do {
            for (Edge* edge : cur_cell->GetEdges()) {
                // Set only the start distance.
                BiomeNode* node = dynamic_cast<BiomeNode*>(edge->GetStart(cur_cell));
                float32 dist_to_edge = m_mountains.back().DistanceTo(node->GetPosition());
                float32 height = 1.0f;
                auto node_height = set_heights.find(node);
                if (node_height != set_heights.end()) {
                    node_height->second = height;
                } else {
                    set_heights[node] = height;
                }
            }
            cur_cell = cur_cell->GetNextCell(m_mountains.back());
        } while (nullptr != cur_cell);

        std::map<BiomeNode*, int32> matrix_node_map;
        int32 cur_ind = 0;
        // assign indices to all the nodes that need to be calculated.
        for (Node* node : m_nodes) {
            BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
            if (bnode->IsLand() && set_heights.find(bnode) == set_heights.end()) {
                matrix_node_map[bnode] = cur_ind++;
            }
        }

        // Construct the sparse matrix and the result vector
        VectorNf b(cur_ind, 0.0f);
        Math::SparseMatrix A(cur_ind);
        for (auto node : matrix_node_map) {
            Math::MatrixRow row;
            Math::SparseMatrixCell diag_cell;
            // set up the diagnal
            diag_cell.col = node.second;
            diag_cell.val = 0.0f;
            // calculate both the b vector and the cell
            for (auto edge : node.first->GetEdges()) {
                auto other = dynamic_cast<BiomeNode*>(edge->GetOtherNode(node.first));
                auto set_height = set_heights.find(other);
                Math::SparseMatrixCell cell;
                cell.val = -1.0f / edge->Length();
                diag_cell.val -= cell.val;
                if (set_height != set_heights.end()) {
                    b[node.second] -= set_height->second * cell.val;
                } else {
                    auto other_ind = matrix_node_map.find(other)->second;
                    cell.col = other_ind;
                    row.push_back(cell);
                }
            }
            row.push_back(diag_cell);
            A.SetRow(node.second, row);
        }

        VectorNf x;
        if (!A.ConjagateGradientSolve(b, x)) {
            assert(0 && "Failed to solve Laplacian.");
            LOG(Debugability::DBG_LVL_CRITICAL, "Failed to solve Laplacian.");
        }


        // copy over the result vector to the data
        for (auto height : set_heights) {
            if (height.second > height.first->GetHeight()) {
                height.first->SetHeight(height.second);
            }
        }

        for (auto node : matrix_node_map) {
            if (x[node.second] > node.first->GetHeight()) {
                node.first->SetHeight(x[node.second]);
            }
        }
    }

    float32 new_max_height = 0.0f;
    for (Node* node : m_nodes) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        float32 shore_dist = bnode->GetShoreDistance();
        if (bnode->IsOcean()) {
            bnode->SetHeight(max(-10, -shore_dist));
        } else {
            float32 height_scale = min(1.0f, sqrt(shore_dist / m_max_distance));
            float32 height = height_scale * bnode->GetHeight();
            bnode->SetHeight(height);
            if (height > new_max_height) {
                new_max_height = height;
            }
        }
    }
    for (Node* node : m_nodes) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        if (bnode->IsLand()) {
            float32 height = bnode->GetHeight();
            if (height == 0) {
                bnode->SetHeight(bnode->GetShoreDistance() / (m_max_distance) * m_max_height);
            } else {
                height *= height;
                bnode->SetHeight(height * m_max_height / new_max_height);
            }
        }
    }
}

void BiomeMap::AssignTemperature() {
    // Determine if this continent is in the North or South hemisphere.
    std::uniform_int_distribution<uint32> hemisphere_dist(0, 1);
    m_is_south_hemisphere = (1 == hemisphere_dist(m_rng));

    const float32 kMinTemp = -15.0f;
    const float32 kMaxTemp = 55.0f;
    float32 temp_spread = kMaxTemp - kMinTemp;

    for (auto* node : m_nodes) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        assert(nullptr != bnode);

        // the temperature is a simple mapping function for the height and the
        // latitude.
        float32 latitude = bnode->GetPosition().y;
        float32 height = bnode->GetHeight();

        // normalize the positions.
        height = height / m_max_height;
        latitude = (latitude - m_min_bounds.y) / (m_max_bounds.y - m_min_bounds.y);
        // inverse the latitude if in the south
        if (!m_is_south_hemisphere) {
            latitude = 1.0f - latitude;
        }
        latitude = latitude * 1.4f - 0.2f;
        height = height - 0.15f;

        bnode->SetTemperature((latitude - pow(height, 3.0f)) * temp_spread + kMinTemp);
    }
}

void BiomeMap::AssignMoisture() {
    {
        // start by finding the primary winds direction
        std::uniform_real_distribution<float32> wind_angle_dist(0.0f, PI2);
        float32 wind_angle = wind_angle_dist(m_rng);
        m_wind_dir.x = cosf(wind_angle);
        m_wind_dir.y = sinf(wind_angle);
    }

    const float32 kMaxMoisture = 10.0f;
    const float32 kMinMoisture = 0.0f;
    const float32 kLandFalloff = 0.02f; // the amount of moisture removed by land / km
    const float32 kWaterAdd = 0.075f; // the amount of moisture added by water / km
    const float32 kHeightFalloff = 0.125f; // the amount of moisture removed by height / km

    std::queue<BiomeNode*> process_nodes;
    for (auto* node : m_nodes) {
        const Vector2f& node_pos = node->GetPosition();
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        assert(nullptr != bnode);
        // this may add the corner nodes twice, but that isn't really a big deal.
        // The second time (and first really) the node should be ignored.
        if ((node_pos.x == m_min_bounds.x && m_wind_dir.x > 0.0f) ||
            (node_pos.x == m_max_bounds.x && m_wind_dir.x < 0.0f) ||
            (node_pos.y == m_min_bounds.y && m_wind_dir.y > 0.0f) ||
            (node_pos.y == m_max_bounds.y && m_wind_dir.y < 0.0f)) {
            // assign the edge moisture to max moisture value (10.0f) and add it to the processing queue
            bnode->SetMoisture(kMaxMoisture);
            process_nodes.push(bnode);
        }
    }

    // while the queue is still not empty
    while (!process_nodes.empty()) {
        BiomeNode* cur_node = process_nodes.front();
        process_nodes.pop();

        // for each edge in the forward direction of the wind
        for (auto* edge : cur_node->GetEdges()) {
            Geometry::Segment edge_seg = edge->GetSegment(cur_node);
            float32 edge_in_dir = edge_seg.GetDirection().Dot(m_wind_dir);
            if (edge_in_dir > 0.0f) {
                BiomeEdge* bedge = dynamic_cast<BiomeEdge*>(edge);
                BiomeNode* other_node = dynamic_cast<BiomeNode*>(bedge->GetOtherNode(cur_node));
                float32 moisture = cur_node->GetMoisture();
                if (bedge->IsLand()) {
                    moisture -= kLandFalloff * edge_seg.Length();
                    if (other_node->GetHeight() > cur_node->GetHeight()) {
                        moisture -= kHeightFalloff * (other_node->GetHeight() - cur_node->GetHeight());
                    }
                } else if (bedge->IsOcean()) {
                    moisture += kWaterAdd * edge_seg.Length();
                }
                moisture = max(min(moisture, kMaxMoisture), kMinMoisture);
                if (moisture > other_node->GetMoisture()) {
                    other_node->SetMoisture(moisture);
                    process_nodes.push(other_node);
                }
            }
        }
    }
}

void BiomeMap::AssignBiomes(const BiomeFactory& biome_factory) {
    for (auto* cell : m_cells) {
        BiomeRegion* region = dynamic_cast<BiomeRegion*>(cell);
        region->CalculateCenterStats();
        region->SetBiome(biome_factory);
    }
}

Node* BiomeMap::CopyNode(Node* other) {
    return new BiomeNode(other);
}

Edge* BiomeMap::CopyEdge(Edge* other) {
    return new BiomeEdge(other);
}

Cell* BiomeMap::CopyCell(Cell* other) {
    return new BiomeRegion(other);
}

struct MapVertex {
    Vector2f pos;
    float32 height;
    float32 temperature;
    float32 moisture;
};

void BiomeMap::Draw(Graphics::Shader* shader) {
    std::vector<Graphics::BufferComponent> comps(2);
    comps[0].attrib_name = "position";
    comps[0].element_count = 2;
    comps[1].attrib_name = "height";
    comps[1].element_count = 1;
    comps[1].offset = Graphics::BUFFER_OFFSET(offsetof(MapVertex, height));
    /*comps[2].attrib_name = "temperature";
    comps[2].element_count = 1;
    comps[2].offset = Graphics::BUFFER_OFFSET(offsetof(MapVertex, temperature));
    comps[3].attrib_name = "moisture";
    comps[3].element_count = 1;
    comps[3].offset = Graphics::BUFFER_OFFSET(offsetof(MapVertex, moisture));*/

    for (auto* cell : m_cells) {
        BiomeRegion* region = dynamic_cast<BiomeRegion*>(cell);
        region->Draw(shader);
    }

    return;

    Graphics::VertexArrayObject mountain_obj(shader);
    MapVertex build_vert;
    build_vert.height = m_max_height + 10.0f;
    build_vert.moisture = 20.0f;
    build_vert.temperature = 20.0f;
    std::vector<MapVertex> mountain_verts;
    mountain_verts.reserve(2 * m_mountains.size() + 2);
    for (auto& mountain : m_mountains) {
        build_vert.pos = mountain.GetStart();
        mountain_verts.push_back(build_vert);
        build_vert.pos = mountain.GetEnd();
        mountain_verts.push_back(build_vert);
    }
    build_vert.pos.Zero();
    mountain_verts.push_back(build_vert);
    build_vert.moisture = 0.0f;
    build_vert.pos.Mul(m_wind_dir, 100.0f);
    mountain_verts.push_back(build_vert);
    int32 ind = mountain_obj.AttachAttributeArray(
        &mountain_verts[0], sizeof(mountain_verts[0]), mountain_verts.size(), false, comps, false);
    mountain_obj.Draw(GL_LINES, mountain_verts.size());
    //Graph::Draw(shader);
}

class NodeShoreDistanceCompare {
public:
    // the processing queue should be ordered by distance. lowest first.
    bool operator()(const BiomeNode* n1, const BiomeNode* n2) {
        return n1->GetShoreDistance() > n2->GetShoreDistance();
    }
};

void BiomeMap::CalculateDistanceToShore() {
    std::queue<BiomeNode*> gray_nodes;
    for (auto* node : m_nodes) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        assert (nullptr != bnode);
        if(bnode->IsShore()) {
            // set all of the shore vertices to distance 0 and add them to the processing queue
            bnode->SetShoreDistance(0.0f);
            gray_nodes.push(bnode);
        } else {
            bnode->SetShoreDistance(-1.0f);
        }
    }
    while (!gray_nodes.empty()) {
        // take the top off of the queue.
        BiomeNode* node = gray_nodes.front();
        gray_nodes.pop();

        // find a white node on the edge
        for (auto* edge : node->GetEdges()) {
            BiomeNode* other = dynamic_cast<BiomeNode*>(edge->GetOtherNode(node));
            assert(nullptr != other);
            // get the length of the edge.
            float32 edge_len = edge->Length();
            float32 new_height = node->GetShoreDistance() + edge_len;
            if (other->GetShoreDistance() < 0.0f || new_height < other->GetShoreDistance()) {
                other->SetShoreDistance(new_height);
                gray_nodes.push(other);
            }
        }
    }
    m_max_distance = 0.0f;
    for (auto* node : m_nodes) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(node);
        if (bnode->IsLand() && bnode->GetShoreDistance() > m_max_distance) {
            m_max_distance = bnode->GetShoreDistance();
        }
    }
}

void BiomeMap::BuildRegionGrid() {
    m_region_grid_width = static_cast<int32>(ceil(m_max_bounds.x - m_min_bounds.x)) / kCellSize + 1;
    m_region_grid_height = static_cast<int32>(ceil(m_max_bounds.y - m_min_bounds.y)) / kCellSize + 1;

    float32 cell_x = m_min_bounds.x;
    float32 cell_y = m_min_bounds.y;

    m_region_grid.resize(m_region_grid_width * m_region_grid_height);
    for (auto* cell : m_cells) {
        Vector2f min_bounds;
        Vector2f max_bounds;
        cell->GetBounds(min_bounds, max_bounds);
        int32 start_row = static_cast<int32>(min_bounds.y - m_min_bounds.y) / kCellSize;
        int32 start_col = static_cast<int32>(min_bounds.x - m_min_bounds.x) / kCellSize;
        int32 end_row = static_cast<int32>(ceil(max_bounds.y - m_min_bounds.y)) / kCellSize;
        int32 end_col = static_cast<int32>(ceil(max_bounds.x - m_min_bounds.x)) / kCellSize;

        assert(start_row <= end_row);
        assert(start_col <= end_col);
        assert(end_row < m_region_grid_height);
        assert(end_col < m_region_grid_width);
        for (int32 row = start_row; row <= end_row; ++row) {
            for (int32 col = start_col; col <= end_col; ++col) {
                m_region_grid[row * m_region_grid_width + col].push_back(cell);
            }
        }
    }
}

}
}
}