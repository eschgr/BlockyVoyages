#pragma once

#include <string>
#include <vector>

#include "gl/glcorearb.h"

namespace BlockyVoyages {
namespace Graphics {

class PNGReader {
public:
    PNGReader(const std::string& filename);
    ~PNGReader();

    bool ReadImage(std::vector<ubyte>& buffer);

    // Get image information after reading
    int32 GetWidth() { return m_width; }
    int32 GetHeight() { return m_height; }
    GLenum GetDataType() { return GL_UNSIGNED_BYTE; }
    GLenum GetComponentType() { return m_component_type; }

private:
    bool CheckIfPNG();

    FILE* m_fp;
    std::string m_filename;
    
    GLenum m_component_type;
    int32 m_width;
    int32 m_height;
};

}
}