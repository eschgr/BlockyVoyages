#include "stdafx.h"

#include "Biome.h"
#include "BiomeEdge.h"
#include "BiomeFactory.h"
#include "BiomeNode.h"
#include "BiomeRegion.h"

#include "../../Graphics/OpenGL.h"
#include "../../Graphics/Shader.h"
#include "../../Graphics/VertexArrayObject.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

BiomeRegion::BiomeRegion(const Cell* other)
    : Graph::Cell(other),
      m_is_water(false),
      m_vertex_array(nullptr),
      m_height(0.0f),
      m_temperature(0.0f),
      m_moisture(0.0f),
      m_biome(nullptr)
{}

BiomeRegion::~BiomeRegion() {
    delete m_vertex_array;
    delete m_biome;
}

struct MapVertex {
    Vector2f pos;
    float32 height;
    float32 temperature;
    float32 moisture;
};

void BiomeRegion::CalculateCenterStats() {
    Cell::CalculateCenterStats();

    m_height = 0.0f;
    m_temperature = 0.0f;
    m_moisture = 0.0f;

    for (auto* edge : m_edges) {
        BiomeNode* bnode = dynamic_cast<BiomeNode*>(edge->GetStart(this));
        m_height += bnode->GetHeight();
        m_temperature += bnode->GetTemperature();
        m_moisture += bnode->GetMoisture();
    }

    m_height /= m_edges.size();
    m_temperature /= m_edges.size();
    m_moisture /= m_edges.size();
}

void BiomeRegion::SetBiome(const BiomeFactory& biome_factory) {
    if (m_is_water) {
        m_biome = biome_factory.GetSeaBiome(this->AtBoundary());
    } else {
        m_biome = biome_factory.GetBiome(m_moisture, m_temperature);
    }
}

float32 BiomeRegion::GetHeightAt(const Vector2f& point) {
    for (auto* edge : m_edges) {
        BiomeNode* start = dynamic_cast<BiomeNode*>(edge->GetStart(this));
        BiomeNode* end = dynamic_cast<BiomeNode*>(edge->GetEnd(this));
        const Vector2f& pnt1 = start->GetPosition();
        const Vector2f& pnt2 = end->GetPosition();
        const Vector2f& pnt3 = m_center;
    
        float32 det_t = (pnt2.y - m_center.y) * (pnt1.x - m_center.x) + 
                        (m_center.x - pnt2.x) * (pnt1.y - m_center.y);
        float32 l1 = ((pnt2.y - m_center.y) * (point.x - m_center.x) +
                      (m_center.x - pnt2.x) * (point.y - m_center.y));
        float32 l2 = ((m_center.y - pnt1.y) * (point.x - m_center.x) +
                      (pnt1.x - m_center.x) * (point.y - m_center.y));
        if ((l1 >= 0.0f && l1 <= det_t) &&
            (l2 >= 0.0f && l2 <= det_t)) {
            float32 l3 = det_t - l1 - l2;
            return (l1 * start->GetHeight() + l2 * end->GetHeight() + l3 * m_height) / det_t;
        }
    }
    return -1000.0f;

    // find which edge the point is closest to.
    /*float32 min_dist = std::numeric_limits<float>::infinity();
    Graph::Edge* closest_edge = nullptr;
    for (auto* edge : m_edges) {
        float32 dist_to_edge = edge->GetSegment().DistanceTo(point);
        if (dist_to_edge < min_dist) {
            closest_edge = edge;
            min_dist = dist_to_edge;
        }
    }
    assert(nullptr != closest_edge);
    BiomeNode* start = dynamic_cast<BiomeNode*>(closest_edge->GetStart(this));
    BiomeNode* end = dynamic_cast<BiomeNode*>(closest_edge->GetEnd(this));
    const Vector2f& pnt1 = start->GetPosition();
    const Vector2f& pnt2 = end->GetPosition();
    const Vector2f& pnt3 = m_center;
    
    float32 det_t = (pnt2.y - m_center.y) * (pnt1.x - m_center.x) + 
                    (m_center.x - pnt2.x) * (m_center.y - pnt3.y);
    float32 l1 = ((pnt2.y - m_center.y) * (point.x - m_center.x) +
                  (m_center.x - pnt2.x) * (point.y - m_center.y));
    float32 l2 = ((m_center.y - pnt1.y) * (point.x - m_center.x) +
                  (pnt1.x - m_center.x) * (point.y - m_center.y));
    float32 l3 = det_t - l1 - l2;
    if (l1 < 0.0f || l1 > det_t ||
        l2 < 0.0f || l2 > det_t) {
        return m_height;
    }
    assert(l1 >= 0.0f && l1 <= det_t);
    assert(l2 >= 0.0f && l2 <= det_t);

    return (l1 * start->GetHeight() + l2 * end->GetHeight() + l3 * m_height) / det_t;*/
}

void BiomeRegion::Draw(Graphics::Shader* shader) {
    if (nullptr == m_biome) {
        return;
    }
    if (nullptr != m_vertex_array) {
        m_vertex_array->Bind();
    } else {
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

        m_vertex_array = new Graphics::VertexArrayObject(shader);
        std::vector<MapVertex> cell_verts;
        cell_verts.reserve(m_edges.size());
        for (auto* edge : m_edges) {
            BiomeNode* start_node = dynamic_cast<BiomeNode*>(edge->GetStart(this));
            float32 moved_dist = start_node->GetHeight();
            if (m_is_water) {
                moved_dist = min(-0.00001f, moved_dist);
            }
            MapVertex vert = {start_node->GetPosition(),
                              moved_dist,
                              m_temperature,
                              m_moisture};
            cell_verts.push_back(vert);
        }
        int32 ind = m_vertex_array->AttachAttributeArray(
            &cell_verts[0], sizeof(cell_verts[0]), cell_verts.size(), false, comps, false);
    }
    shader->setUniform("biome_color", m_biome->GetMapColor());
    m_vertex_array->Draw(GL_TRIANGLE_FAN, m_edges.size());
}

}
}
}