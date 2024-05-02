#pragma once

#include <string>
#include <vector>

#include "../types.h"

namespace BlockyVoyages {
namespace Graphics {

class Shader;

class TextureArray {
public:
    TextureArray(const std::vector<std::string>& filenames, bool make_mipsmaps);

    ~TextureArray();

    void Bind();
private:
    bool m_has_mipmaps;

    // used for verification of user provided parameters
    int32 m_size;
    int32 m_width;
    int32 m_height;
    int32 m_comp_count;

    uint32 m_tex_handle;

    void SetupTexture();
    void LoadTexture(const void* data, int32 index, GLenum format, GLenum dataType);
    bool SetTexture(const std::string& filename, int32 index);
};

}
}