#pragma once

#include "OpenGL.h"

namespace BlockyVoyages {
namespace Graphics {

class VertexBuffer {
public:
    VertexBuffer(GLenum type,
                    void* buffer,
                    int32 count,
                    int32 stride,
                    bool dynamic);
    ~VertexBuffer();

    void ResetBuffer(void* buffer, int32 count);
    void* GetBufferPointer();
    void FinishBufferWrite();

    // set instancing for this buffer. This is a one way direction.
    void SetAttributeLocations(const std::vector<int32>& locs) { m_attrib_locs = locs; }
    void SetInstanced();

    bool IsInstanced() { return m_instanced; }

private:
    uint32 m_handle;
    bool m_dynamic;
    bool m_instanced;
    int32 m_stride;
    void* m_buffer_ptr;
    GLenum m_type;
    std::vector<int32> m_attrib_locs;
};


}
}