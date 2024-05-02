#include "stdafx.h"

#include "OpenGL.h"

#include "Shader.h"
#include "../types.h"
#include "../Debugability/Logger.h"

#pragma comment(lib, "OpenGL32.lib")

// State and state requests
PFNGLGETSTRINGIPROC glGetStringi = nullptr;

// Texturing
PFNGLACTIVETEXTUREARBPROC glActiveTexture = nullptr;
PFNGLTEXIMAGE3DPROC glTexImage3D = nullptr;
PFNGLTEXSUBIMAGE3DPROC glTexSubImage3D = nullptr;
PFNGLTEXSTORAGE2DPROC glTexStorage2D = nullptr;
PFNGLTEXSTORAGE3DPROC glTexStorage3D = nullptr;

// vertex buffer support
PFNGLGENBUFFERSARBPROC glGenBuffers = nullptr;
PFNGLBINDBUFFERARBPROC glBindBuffer = nullptr;
PFNGLBUFFERDATAARBPROC glBufferData = nullptr;
PFNGLDELETEBUFFERSARBPROC glDeleteBuffers = nullptr;
PFNGLGENVERTEXARRAYSPROC glGenVertexArrays = nullptr;
PFNGLBINDVERTEXARRAYPROC glBindVertexArray = nullptr;
PFNGLDELETEVERTEXARRAYSPROC glDeleteVertexArrays = nullptr;
PFNGLDRAWARRAYSINSTANCEDPROC glDrawArraysInstanced = nullptr;
PFNGLDRAWELEMENTSINSTANCEDPROC glDrawElementsInstanced = nullptr;
PFNGLVERTEXATTRIBDIVISORPROC glVertexAttribDivisor = nullptr;
PFNGLMAPBUFFERPROC glMapBuffer = nullptr;
PFNGLUNMAPBUFFERPROC glUnmapBuffer = nullptr;

// Vsync control
PFNWGLSWAPINTERVALEXTPROC wglSwapIntervalEXT = nullptr;

// shaders
PFNGLCREATESHADERPROC glCreateShader = nullptr;
PFNGLSHADERSOURCEPROC glShaderSource = nullptr;
PFNGLCOMPILESHADERPROC glCompileShader = nullptr;
PFNGLCREATEPROGRAMPROC glCreateProgram = nullptr;
PFNGLATTACHSHADERPROC glAttachShader = nullptr;
PFNGLLINKPROGRAMPROC glLinkProgram = nullptr;
PFNGLUSEPROGRAMPROC glUseProgram = nullptr;
PFNGLGETSHADERIVPROC glGetShaderiv = nullptr;
PFNGLGETPROGRAMIVPROC glGetProgramiv = nullptr;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog = nullptr;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog = nullptr;
PFNGLDETACHSHADERPROC glDetachShader = nullptr;
PFNGLDELETESHADERPROC glDeleteShader = nullptr;
PFNGLDELETEPROGRAMPROC glDeleteProgram = nullptr;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation = nullptr;
PFNGLUNIFORM1FPROC glUniform1f = nullptr;
PFNGLUNIFORM2FPROC glUniform2f = nullptr;
PFNGLUNIFORM3FPROC glUniform3f = nullptr;
PFNGLUNIFORM4FPROC glUniform4f = nullptr;
PFNGLUNIFORM1IPROC glUniform1i = nullptr;
PFNGLUNIFORM2IPROC glUniform2i = nullptr;
PFNGLUNIFORM3IPROC glUniform3i = nullptr;
PFNGLUNIFORM4IPROC glUniform4i = nullptr;
PFNGLUNIFORM1UIPROC glUniform1ui = nullptr;
PFNGLUNIFORM2UIPROC glUniform2ui = nullptr;
PFNGLUNIFORM3UIPROC glUniform3ui = nullptr;
PFNGLUNIFORM4UIPROC glUniform4ui = nullptr;
PFNGLUNIFORM1FVPROC glUniform1fv = nullptr;
PFNGLUNIFORM2FVPROC glUniform2fv = nullptr;
PFNGLUNIFORM3FVPROC glUniform3fv = nullptr;
PFNGLUNIFORM4FVPROC glUniform4fv = nullptr;
PFNGLUNIFORM1IVPROC glUniform1iv = nullptr;
PFNGLUNIFORM2IVPROC glUniform2iv = nullptr;
PFNGLUNIFORM3IVPROC glUniform3iv = nullptr;
PFNGLUNIFORM4IVPROC glUniform4iv = nullptr;
PFNGLUNIFORM1UIVPROC glUniform1uiv = nullptr;
PFNGLUNIFORM2UIVPROC glUniform2uiv = nullptr;
PFNGLUNIFORM3UIVPROC glUniform3uiv = nullptr;
PFNGLUNIFORM4UIVPROC glUniform4uiv = nullptr;
PFNGLUNIFORMMATRIX2FVPROC glUniformMatrix2fv = nullptr;
PFNGLUNIFORMMATRIX3FVPROC glUniformMatrix3fv = nullptr;
PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv = nullptr;
PFNGLGETATTRIBLOCATIONPROC glGetAttribLocation = nullptr;
PFNGLVERTEXATTRIB1FPROC glVertexAttrib1f = nullptr;
PFNGLVERTEXATTRIB2FPROC glVertexAttrib2f = nullptr;
PFNGLVERTEXATTRIB3FPROC glVertexAttrib3f = nullptr;
PFNGLVERTEXATTRIB4FPROC glVertexAttrib4f = nullptr;
PFNGLVERTEXATTRIB1FVPROC glVertexAttrib1fv = nullptr;
PFNGLVERTEXATTRIB2FVPROC glVertexAttrib2fv = nullptr;
PFNGLVERTEXATTRIB3FVPROC glVertexAttrib3fv = nullptr;
PFNGLVERTEXATTRIB4FVPROC glVertexAttrib4fv = nullptr;
PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray = nullptr;
PFNGLDISABLEVERTEXATTRIBARRAYPROC glDisableVertexAttribArray = nullptr;
PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer = nullptr;
PFNGLVERTEXATTRIBIPOINTERPROC glVertexAttribIPointer = nullptr;
PFNGLVERTEXATTRIBLPOINTERPROC glVertexAttribLPointer = nullptr;
PFNGLGETACTIVEUNIFORMPROC glGetActiveUniform = nullptr;
PFNGLBINDATTRIBLOCATIONPROC glBindAttribLocation = nullptr;

namespace BlockyVoyages {
namespace Graphics {

OpenGL* g_opengl = nullptr;

OpenGL::OpenGL(HWND hWnd)
            : m_hWnd(hWnd),
              m_hGLRC(0),
              m_initialized(false),
              m_activeShader(nullptr) {
    if (!SetupDC()) {
        LOG(Debugability::DBG_LVL_ALWAYS_TERMINAL, "Failed to setup device context.");
        return;
    }

    if (!CreateContext()) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to create OpenGL context.");
        return;
    }

    // clear the front buffer so we don't have the desktop color
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glDrawBuffer(GL_FRONT);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawBuffer(GL_BACK);

    // must get this procedure address now or else we can't continue
    if (!LoadProcAddress(glGetStringi, "glGetStringi")) {
        return;
    }

    // check for support for the extensions that we care about
    LoadExtensionSupport();
    if (!LoadTexturingProcs()) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not load procedures for texturing.");
        return;
    }
    if (!LoadBufferObjectProcs()) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not load procedures for buffer objects.");
        return;
    }
    if (!LoadShaderProcs()) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not load procedures for shaders.");
        return;
    }

    if (!LoadProcAddress(wglSwapIntervalEXT, "wglSwapIntervalEXT")) {
        return;
    }
    
    LOG(Debugability::DBG_LVL_MESSAGE, "OpenGL initialized.");
    m_initialized = true;
}

OpenGL::~OpenGL() {
    // shut down the texture manager
    if (m_hGLRC) {
        if (!wglMakeCurrent(nullptr, nullptr)) {
            LOG(Debugability::DBG_LVL_CRITICAL, "Release Of DC And RC Failed.");
        }

        if (!wglDeleteContext(m_hGLRC)) {
            LOG(Debugability::DBG_LVL_CRITICAL, "Unable to delete rendering context.");
        }
    }

    if (m_hDC) {
        if (!ReleaseDC(m_hWnd, m_hDC)) {
            LOG(Debugability::DBG_LVL_CRITICAL, "Unable to release DC.");
        }
    }
}

void OpenGL::present() {
    SwapBuffers(m_hDC);
}

void OpenGL::SetActiveShader(Shader* shader) {
    // if the shader being set is the same as the previous one, just return
    if (m_activeShader == shader) {
        return;
    }
    if (nullptr == shader) {
        glUseProgram(0);
    } else {
        shader->bind();
    }
    m_activeShader = shader;
}

bool OpenGL::SetupDC() {
    int color_depth = 32;

    m_hDC = GetDC(m_hWnd);
    if (0 == m_hDC) {
        LOG(Debugability::DBG_LVL_ALWAYS_TERMINAL, "Could not get window device context.");
        return false;
    }

    PIXELFORMATDESCRIPTOR pfd = {
        sizeof(PIXELFORMATDESCRIPTOR),    // Size Of This Pixel Format Descriptor
            1,                            // Version Number
            PFD_DRAW_TO_WINDOW |          // Format Must Support Window
            PFD_SUPPORT_OPENGL |          // Format Must Support OpenGL
            PFD_DOUBLEBUFFER,             // Must Support Double Buffering
            PFD_TYPE_RGBA,                // Request An RGBA Format
            color_depth,                   // Select Our Color Depth
            0, 0, 0, 0, 0, 0,             // Color Bits Ignored
            0,                            // No Alpha Buffer
            0,                            // Shift Bit Ignored
            0,                            // No Accumulation Buffer
            0, 0, 0, 0,                   // Accumulation Bits Ignored
            16,                           // 16Bit Z-Buffer (Depth Buffer)  
            0,                            // No Stencil Buffer
            0,                            // No Auxiliary Buffer
            PFD_MAIN_PLANE,               // Main Drawing Layer
            0,                            // Reserved
            0, 0, 0                       // Layer Masks Ignored
    };
    unsigned int pixel_format = ChoosePixelFormat(m_hDC, &pfd);
    if (0 == pixel_format) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not find requested pixel format.");
        return false;
    }

    if (!SetPixelFormat(m_hDC, pixel_format, &pfd)) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Could not set pixel format.");
        return false;
    }
    return true;
}

bool OpenGL::CreateContext() {
    // Create a temporary old style context.
    HGLRC old_context = wglCreateContext(m_hDC);
    if (0 == old_context) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to create old OpenGL context.");
        return false;
    }
    if (!wglMakeCurrent(m_hDC, old_context)) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to make old context current.");
        return false;
    }

    // check for support of OpenGL 3.3+.
    int major, minor;
    glGetIntegerv(GL_MAJOR_VERSION, &major);
    glGetIntegerv(GL_MINOR_VERSION, &minor);
    if (major < 3 || (major == 3 && minor < 3)) {
        LOG(Debugability::DBG_LVL_CRITICAL, "OpenGL Version 3.3 or higher not supported.");
        LOG(Debugability::DBG_LVL_CRITICAL, "Found version %d.%d instead.", major, minor);
        return false;
    }

    // Create the OpenGL 3.3 core profile context
    int attribs[] =
    {
        WGL_CONTEXT_MAJOR_VERSION_ARB, major,
        WGL_CONTEXT_MINOR_VERSION_ARB, minor, 
        WGL_CONTEXT_FLAGS_ARB, WGL_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB,
        WGL_CONTEXT_PROFILE_MASK_ARB, WGL_CONTEXT_CORE_PROFILE_BIT_ARB,
        0
    };
    PFNWGLCREATECONTEXTATTRIBSARBPROC wglCreateContextAttribsARB = nullptr;
    if (!LoadProcAddress(wglCreateContextAttribsARB, "wglCreateContextAttribsARB")) {
        return false;
    }
    m_hGLRC = wglCreateContextAttribsARB(m_hDC, 0, attribs);
    if (!m_hGLRC) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to create OpenGL 3.3 core context.");
        return false;
    }

    // clean up the old context
    wglMakeCurrent(0, 0);
    wglDeleteContext(old_context);

    if (!wglMakeCurrent(m_hDC, m_hGLRC)) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to make OpenGL 3.3 context current.");
        return false;
    }
    return true;
}

void OpenGL::LoadExtensionSupport() {
    int32 num_exts, i;
    glGetIntegerv(GL_NUM_EXTENSIONS, &num_exts);
    LOG(Debugability::DBG_LVL_NORMAL, "Supported Extensions:");
    for (i = 0; i < num_exts; i++) {
        const char* ext_str = reinterpret_cast<const char*>(glGetStringi(GL_EXTENSIONS, i));
        LOG(Debugability::DBG_LVL_NORMAL, "   %s", ext_str);
        m_extsSupported[ext_str] = true;
    }
}

bool OpenGL::isExtensionSupported(const char* extension) const {
    auto supported = m_extsSupported.find(extension);
    if (supported != m_extsSupported.end()) {
        return false;
    }
    return supported->second;
}

template <class Type>
bool OpenGL::LoadProcAddress(Type &func_ptr, const std::string& name) {
    func_ptr = reinterpret_cast<Type>(wglGetProcAddress(name.c_str()));
    if (nullptr == func_ptr) {
        LOG(Debugability::DBG_LVL_HIGH, "Could not load address for method '%s'.", name.c_str());
        return false;
    }
    return true;
}

bool OpenGL::LoadTexturingProcs() {
    return LoadProcAddress(glActiveTexture, "glActiveTexture") &&
           LoadProcAddress(glTexImage3D, "glTexImage3D") &&
           LoadProcAddress(glTexSubImage3D, "glTexSubImage3D") &&
           LoadProcAddress(glTexStorage2D, "glTexStorage2D") &&
           LoadProcAddress(glTexStorage3D, "glTexStorage3D");
}

bool OpenGL::LoadBufferObjectProcs() {
    return LoadProcAddress(glGenBuffers, "glGenBuffers") &&
           LoadProcAddress(glBindBuffer, "glBindBuffer") &&
           LoadProcAddress(glBufferData, "glBufferData") &&
           LoadProcAddress(glDeleteBuffers, "glDeleteBuffers") &&
           LoadProcAddress(glBindAttribLocation, "glBindAttribLocation") &&
           LoadProcAddress(glGenVertexArrays, "glGenVertexArrays") &&
           LoadProcAddress(glBindVertexArray, "glBindVertexArray") &&
           LoadProcAddress(glDeleteVertexArrays, "glDeleteVertexArrays") &&
           LoadProcAddress(glDrawArraysInstanced, "glDrawArraysInstanced") &&
           LoadProcAddress(glDrawElementsInstanced, "glDrawElementsInstanced") &&
           LoadProcAddress(glVertexAttribDivisor, "glVertexAttribDivisor") &&
           LoadProcAddress(glMapBuffer, "glMapBuffer") &&
           LoadProcAddress(glUnmapBuffer, "glUnmapBuffer");
}

bool OpenGL::LoadShaderProcs() {
    return (LoadProcAddress(glCreateShader, "glCreateShader") &&
            LoadProcAddress(glShaderSource, "glShaderSource") &&
            LoadProcAddress(glCompileShader, "glCompileShader") &&
            LoadProcAddress(glCreateProgram, "glCreateProgram") &&
            LoadProcAddress(glAttachShader, "glAttachShader") &&
            LoadProcAddress(glLinkProgram, "glLinkProgram") &&
            LoadProcAddress(glUseProgram, "glUseProgram") &&
            LoadProcAddress(glGetShaderiv, "glGetShaderiv") &&
            LoadProcAddress(glGetProgramiv, "glGetProgramiv") &&
            LoadProcAddress(glGetShaderInfoLog, "glGetShaderInfoLog") &&
            LoadProcAddress(glGetProgramInfoLog, "glGetProgramInfoLog") &&
            LoadProcAddress(glDetachShader, "glDetachShader") &&
            LoadProcAddress(glDeleteShader, "glDeleteShader") &&
            LoadProcAddress(glDeleteProgram, "glDeleteProgram") &&
            LoadProcAddress(glGetUniformLocation, "glGetUniformLocation") &&
            LoadProcAddress(glUniform1f, "glUniform1f") &&
            LoadProcAddress(glUniform2f, "glUniform2f") &&
            LoadProcAddress(glUniform3f, "glUniform3f") &&
            LoadProcAddress(glUniform4f, "glUniform4f") &&
            LoadProcAddress(glUniform1i, "glUniform1i") &&
            LoadProcAddress(glUniform2i, "glUniform2i") &&
            LoadProcAddress(glUniform3i, "glUniform3i") &&
            LoadProcAddress(glUniform4i, "glUniform4i") &&
            LoadProcAddress(glUniform1ui, "glUniform1ui") &&
            LoadProcAddress(glUniform2ui, "glUniform2ui") &&
            LoadProcAddress(glUniform3ui, "glUniform3ui") &&
            LoadProcAddress(glUniform4ui, "glUniform4ui") &&
            LoadProcAddress(glUniform1fv, "glUniform1fv") &&
            LoadProcAddress(glUniform2fv, "glUniform2fv") &&
            LoadProcAddress(glUniform3fv, "glUniform3fv") &&
            LoadProcAddress(glUniform4fv, "glUniform4fv") &&
            LoadProcAddress(glUniform1iv, "glUniform1iv") &&
            LoadProcAddress(glUniform2iv, "glUniform2iv") &&
            LoadProcAddress(glUniform3iv, "glUniform3iv") &&
            LoadProcAddress(glUniform4iv, "glUniform4iv") &&
            LoadProcAddress(glUniform1uiv, "glUniform1uiv") &&
            LoadProcAddress(glUniform2uiv, "glUniform2uiv") &&
            LoadProcAddress(glUniform3uiv, "glUniform3uiv") &&
            LoadProcAddress(glUniform4uiv, "glUniform4uiv") &&
            LoadProcAddress(glUniformMatrix2fv, "glUniformMatrix2fv") &&
            LoadProcAddress(glUniformMatrix3fv, "glUniformMatrix3fv") &&
            LoadProcAddress(glUniformMatrix4fv, "glUniformMatrix4fv") &&
            LoadProcAddress(glGetAttribLocation, "glGetAttribLocation") &&
            LoadProcAddress(glVertexAttrib1f, "glVertexAttrib1f") &&
            LoadProcAddress(glVertexAttrib2f, "glVertexAttrib2f") &&
            LoadProcAddress(glVertexAttrib3f, "glVertexAttrib3f") &&
            LoadProcAddress(glVertexAttrib4f, "glVertexAttrib4f") &&
            LoadProcAddress(glVertexAttrib1fv, "glVertexAttrib1fv") &&
            LoadProcAddress(glVertexAttrib2fv, "glVertexAttrib2fv") &&
            LoadProcAddress(glVertexAttrib3fv, "glVertexAttrib3fv") &&
            LoadProcAddress(glVertexAttrib4fv, "glVertexAttrib4fv") &&
            LoadProcAddress(glEnableVertexAttribArray, "glEnableVertexAttribArray") &&
            LoadProcAddress(glDisableVertexAttribArray, "glDisableVertexAttribArray") &&
            LoadProcAddress(glVertexAttribPointer, "glVertexAttribPointer") &&
            LoadProcAddress(glVertexAttribIPointer, "glVertexAttribIPointer") &&
            LoadProcAddress(glVertexAttribLPointer, "glVertexAttribLPointer") &&
            LoadProcAddress(glGetActiveUniform, "glGetActiveUniform"));
}

}
}