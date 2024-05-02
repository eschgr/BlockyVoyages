#pragma once

#include <string>
#include <vector>

#include <ft2build.h>
#include FT_FREETYPE_H

#include "../types.h"
#include "Font.h"

namespace BlockyVoyages {
namespace Graphics {

class Shader;

class FontManager {
public:
    FontManager();
    ~FontManager();

    // getting information from the font manager
    bool IsInitialize() { return m_initialized; }

    // font management
    Font* createFont(std::string fontFile, int32 pntSize);
    Font* getFont(std::string fontFile, int32 pntSize);
    Font* getDefaultFont() { return m_fonts[0]; }

    // will not release the default font
    void releaseFont(std::string fontFile, int32 pntSize);
        
    // will not release the default font
    void releaseFont(Font* font);

private:
    FT_Library m_lib;  // the FreeType library
    std::vector<Font*> m_fonts;
    bool m_initialized;
    Shader* m_fontShader;
};

extern FontManager* g_fontManager;

}
}