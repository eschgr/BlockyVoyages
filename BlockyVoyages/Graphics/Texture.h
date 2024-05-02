#pragma once

#include "../types.h"

#include "OpenGL.h"

#include <string>

namespace BlockyVoyages {
namespace Graphics {

// Manages the resources and state corresponding to a single texture.
// A Texture is meant to be used for putting an image on objects. For this
// reason, the internal format is always R8, RG8, RGB8, or RGBA8.
// If large buffers are needed, see DataBuffer (Not created yet.)
class Texture {
public:
    // load the texture from file. Only png files are supported.
    Texture(const std::string& filename, bool makeMips);
    Texture(const void* data, int32 width, int32 height, GLenum format, GLenum dataType, bool makeMips);
    ~Texture(void);

    void bind(void) const;
    void setParameter(GLenum param, uint32 value);

private:
    // OpenGL specific stuff
    uint32 m_texHandle;

    static size_t getElementCount(GLenum format);
    static size_t getDataSize(GLenum dataType);
    static GLenum GetSizedFormat(GLenum dataType, GLenum format);

    static void makeMipMaps(const void* data, int32 width, int32 height, GLenum format, GLenum dataType,
                            std::vector<void*>& out_data);
    static int32 GetMipMapLevelCount(int32 width, int32 height);
    void loadTexture(const void* data, int32 width, int32 height, GLenum format, GLenum dataType, bool makeMips);

    friend class TextureArray;
};

}
}