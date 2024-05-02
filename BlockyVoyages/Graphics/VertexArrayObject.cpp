#include "stdafx.h"

#include "VertexArrayObject.h"

#include <assert.h>

#include "OpenGL.h"
#include "Shader.h"
#include "VertexBuffer.h"

namespace BlockyVoyages {
namespace Graphics {
VertexArrayObject::VertexArrayObject(Shader* shader)
    : m_shader(shader),
      m_indices_ind(-1)
{
    glGenVertexArrays(1, &m_vao_handle);
}

VertexArrayObject::~VertexArrayObject() {
    for (auto* vbo : m_vbos) {
        delete vbo;
    }
    glDeleteVertexArrays(1, &m_vao_handle);
}

int32 VertexArrayObject::AttachAttributeArray(void* buffer,
                                              int32 stride,
                                              int32 count,
                                              bool dynamic,
                                              const std::vector<BufferComponent>& components,
                                              bool instanced) {
    assert (components.size() > 0);
    int32 glError;
    if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
    }
    Bind();
    int32 vbo_ind = m_vbos.size();
    if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
    }
    m_vbos.push_back(new VertexBuffer(GL_ARRAY_BUFFER, buffer, count, stride, dynamic));
    if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
    }
    // store the VBO parameters for later use
    VertexBuffer* new_vbo = m_vbos[vbo_ind];
    std::vector<int32> locs(components.size());
    for (uint32 i = 0; i < components.size(); ++i) {
        locs[i] = m_shader->getAttribLocation(components[i].attrib_name);
        if ((glError = glGetError()) != 0) {
            LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d on comp %d.", glError, i);
        }
        glEnableVertexAttribArray(locs[i]);
        if ((glError = glGetError()) != 0) {
            LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d on comp %d.", glError, i);
        }
        switch (components[i].data_type) {
        case GL_BYTE:
        case GL_UNSIGNED_BYTE:
        case GL_SHORT:
        case GL_UNSIGNED_SHORT:
        case GL_INT:
        case GL_UNSIGNED_INT:
            glVertexAttribIPointer(locs[i],
                                   components[i].element_count,
                                   components[i].data_type,
                                   stride,
                                   components[i].offset);
            break;
        case GL_DOUBLE:
            glVertexAttribLPointer(locs[i],
                                  components[i].element_count,
                                  components[i].data_type,
                                  stride,
                                  components[i].offset);
            break;
        default:
            glVertexAttribPointer(locs[i],
                                  components[i].element_count,
                                  components[i].data_type,
                                  components[i].normalized,
                                  stride,
                                  components[i].offset);
        }
        if ((glError = glGetError()) != 0) {
            LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d on comp %d.", glError, i);
        }
    }
    new_vbo->SetAttributeLocations(locs);
    if (instanced) {
        m_vbos[vbo_ind]->SetInstanced();
    }
    if ((glError = glGetError()) != 0) {
        LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
    }
    return vbo_ind;
}

int32 VertexArrayObject::AttachIndexArray(std::vector<uint32> inds) {
    // only one indexed array is allowed per vertex array object
    assert(m_indices_ind < 0);
    Bind();
    m_indices_ind = static_cast<int32>(m_vbos.size());
    m_vbos.push_back(new VertexBuffer(GL_ELEMENT_ARRAY_BUFFER, &inds[0], inds.size(), sizeof(uint32), false));
    return m_indices_ind;
}

void VertexArrayObject::ResetBufferData(int32 vbo, void* buffer, int32 count) {
    m_vbos[vbo]->ResetBuffer(buffer, count);
}

void* VertexArrayObject::GetBufferPointer(int32 vbo) {
    return m_vbos[vbo]->GetBufferPointer();
}

void VertexArrayObject::FinishBufferWrite(int32 vbo) {
    m_vbos[vbo]->FinishBufferWrite();
}

void VertexArrayObject::Bind() {
    glBindVertexArray(m_vao_handle);
}

void VertexArrayObject::Unbind() {
    glBindVertexArray(0);
}

void VertexArrayObject::Draw(GLenum primType, int32 primCount) {
    if (m_indices_ind >= 0) {
        glDrawElements(primType, primCount, GL_UNSIGNED_INT, nullptr);
    } else {
        glDrawArrays(primType, 0, primCount);
    }
}

void VertexArrayObject::DrawInstanced(GLenum primType, int32 primCount, int32 instances) {
    if (m_indices_ind >= 0) {
        glDrawElementsInstanced(primType, primCount, GL_UNSIGNED_INT, nullptr, instances);
    } else {
        glDrawArraysInstanced(primType, 0, primCount, instances);
    }
}

}
}