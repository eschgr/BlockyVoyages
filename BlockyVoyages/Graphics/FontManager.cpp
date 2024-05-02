#include "stdafx.h"

#include "FontManager.h"

#include <algorithm>
#include <cctype>

#include "Shader.h"
#include "../Debugability/Logger.h"

#if defined(_DEBUG)
#pragma comment(lib, "../freetype-2.5.0.1/objs/win32/vc2010/freetype250_D.lib")
#else
#pragma comment(lib, "../freetype-2.5.0.1/objs/win32/vc2010/freetype250.lib")
#endif

namespace BlockyVoyages {
namespace Graphics {

FontManager* g_fontManager = nullptr;

FontManager::FontManager() 
    : m_lib(nullptr),
      m_initialized(false)
{
    m_fontShader = new Shader();
    if (!m_fontShader->loadFromFile("..\\Data\\Shaders\\text.vs", "..\\Data\\Shaders\\text.fs")) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to load text rendering shader.");
        delete m_fontShader;
        m_fontShader = nullptr;
        return;
    }
    LOG(Debugability::DBG_LVL_NORMAL, "Successfully loaded text rendering shader.");

    if (FT_Init_FreeType(&m_lib)) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not initialize FreeType library.");
        return;
    }
    LOG(Debugability::DBG_LVL_NORMAL, "FreeType Library Initialized.");

    m_fonts.resize(1);
    m_fonts[0] = new Font(m_lib, "arial.ttf", 14, m_fontShader);
    if (!m_fonts[0]->isInitialized()) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not initialize default font 'arial.ttf'.");
        delete m_fonts[0];
        m_fonts[0] = nullptr;
        return;
    }

    LOG(Debugability::DBG_LVL_NORMAL, "Default font 'arial.ttf' initialized.");
    LOG(Debugability::DBG_LVL_NORMAL, "Font Manager Initialized.");
    m_initialized = true;
}

FontManager::~FontManager() {
    for(size_t font = 0; font < m_fonts.size(); ++font) {
        delete m_fonts[font];
    }

    FT_Done_FreeType(m_lib);

    delete m_fontShader;
}

Font* FontManager::createFont(std::string fontFile, int32 pntSize) {
    std::transform(fontFile.begin(), fontFile.end(), fontFile.begin(), std::tolower);

    // either find a place to create the font or find if the font already exists.
    size_t potentSpot = 0;
    for(size_t ind = 0; ind < m_fonts.size(); ++ind) {
        if (m_fonts[ind] == nullptr) {
            if (potentSpot == 0) {
                potentSpot = ind;
            }
        } else if (m_fonts[ind]->getName() == fontFile && 
                   m_fonts[ind]->getPntSize() == pntSize) {
            ++m_fonts[ind]->m_nRefs;
            return m_fonts[ind];
        }
    }

    if (potentSpot == 0) {
        potentSpot = m_fonts.size();
        m_fonts.resize(potentSpot + 1);
    }

    m_fonts[potentSpot] = new Font(m_lib, fontFile, pntSize, m_fontShader);
    if (!m_fonts[potentSpot]->isInitialized()) {
        LOG(Debugability::DBG_LVL_WARNING, "Could not initialize font \"%s\" with size %d.", fontFile.c_str(), pntSize);
        delete m_fonts[potentSpot];
        m_fonts[potentSpot] = nullptr;
        return nullptr;
    }
    LOG(Debugability::DBG_LVL_NORMAL, "Font \"%s\" with size %d initialized.", fontFile.c_str(), pntSize);
    ++m_fonts[potentSpot]->m_nRefs;
    return m_fonts[potentSpot];
}

Font* FontManager::getFont(std::string fontFile, int32 pntSize) {
    std::transform(fontFile.begin(), fontFile.end(), fontFile.begin(), std::tolower);

    // cycle through all the fonts and try to find the right one
    for(size_t ind = 0; ind < m_fonts.size(); ++ind) {
        if (m_fonts[ind] != nullptr && 
            m_fonts[ind]->getName() == fontFile && 
            m_fonts[ind]->getPntSize() == pntSize) 
        {
            ++m_fonts[ind]->m_nRefs;
            return m_fonts[ind];
        }
    }
    return nullptr;
}

void FontManager::releaseFont(std::string fontFile, int32 pntSize) {
    // do not release the default font
    std::transform(fontFile.begin(), fontFile.end(), fontFile.begin(), std::tolower);

    for(size_t ind = 1; ind < m_fonts.size(); ++ind) {
        if (m_fonts[ind]->getName() == fontFile && m_fonts[ind]->getPntSize() == pntSize) {
            --m_fonts[ind]->m_nRefs;
            if (m_fonts[ind]->m_nRefs == 0) {
                LOG(Debugability::DBG_LVL_NORMAL, "Font \"%s\" with size %d destroyed.", fontFile.c_str(), pntSize);
                delete m_fonts[ind];
                m_fonts[ind] = nullptr;
            }
            return;
        }
    }
}

void FontManager::releaseFont(Font* font) {
    if (font == nullptr || font == m_fonts[0]) {
        return;
    }

    // do not release the default font
    for(size_t ind = 1; ind < m_fonts.size(); ++ind) {
        if (m_fonts[ind] == font) {
            --m_fonts[ind]->m_nRefs;
            if (m_fonts[ind]->m_nRefs == 0) {
                LOG(Debugability::DBG_LVL_NORMAL, "Font \"%s\" with size %d destroyed.", font->m_name.c_str(), font->m_pntSize);
                delete m_fonts[ind];
                m_fonts[ind] = nullptr;
            }
            return;
        }
    }
}

}
}