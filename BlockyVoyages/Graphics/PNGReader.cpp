#include "stdafx.h"
#include "PNGReader.h"

#include "../types.h"
#include "../Debugability/Logger.h"
#include "../../lpng163/png.h"
#include "../../lpng163/pngstruct.h"
#include "../../lpng163/pnginfo.h"

#if defined(_DEBUG)
#pragma comment(lib, "../Debug/libpng16.lib")
#pragma comment(lib, "../Debug/zlib.lib")
#else
#pragma comment(lib, "../Release/libpng16.lib")
#pragma comment(lib, "../Release/zlib.lib")
#endif

const static int32 PNG_BYTES_TO_CHECK = 8;

const static png_infopp png_infopp_NULL = (png_infopp)nullptr;
const int32* int_p_NULL = (int*)nullptr;

namespace BlockyVoyages {
namespace Graphics {

static void pngLoadError(png_structp png_ptr, png_const_charp message) {
    printf("Error Loading PNG: %s.", message);
        
    /* Return control to the setjmp point */
    longjmp(png_ptr->jmp_buf_local, 1);
}

static void pngLoadWarning(png_structp png_ptr, png_const_charp message) {
    printf("Warning During PNG Load: %s.", message);
}

void* PNGMalloc(png_structp png_ptr, png_size_t size) {
    return ::operator new(size);
}

void PNGFree(png_structp png_ptr, png_voidp mem) {
    ::operator delete(mem);
}

PNGReader::PNGReader(const std::string& filename)
    : m_fp(fopen(filename.c_str(), "rb")),
      m_filename(filename),
      m_component_type(0),
      m_width(0),
      m_height(0)
{}

PNGReader::~PNGReader() {
    if (nullptr != m_fp) {
        fclose(m_fp);
    }
}

bool PNGReader::ReadImage(std::vector<ubyte>& buffer) {
    png_structp png_ptr;
    png_infop info_ptr;
    unsigned int sig_read = 0;

    if (CheckIfPNG() == false) {
        LOG(Debugability::DBG_LVL_HIGH, "Image '%s' failed PNG check.", m_filename.c_str());
        return false;
    }

    /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply nullptr for the last three parameters.  We also supply the
    * the compiler header file version, so that we know if the application
    * was compiled with a compatible version of the library.  REQUIRED
    */
    png_ptr = png_create_read_struct_2(PNG_LIBPNG_VER_STRING, 
                                        nullptr, 
                                        pngLoadError, 
                                        pngLoadWarning,
                                        nullptr,
                                        PNGMalloc,
                                        PNGFree
                                    );

    if (png_ptr == nullptr) {
        LOG(Debugability::DBG_LVL_HIGH, "Image '%s' is not a well formed PNG image.", m_filename.c_str());
        return false;
    }

    /* Allocate/initialize the memory for image information.  REQUIRED. */
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        LOG(Debugability::DBG_LVL_HIGH, "Image '%s' is not a well formed PNG image.", m_filename.c_str());
        png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
        return false;
    }

    /* Set error handling if you are using the setjmp/longjmp method (this is
    * the normal method of doing things with libpng).  REQUIRED unless you
    * set up your own error handlers in the png_create_read_struct() earlier.
    */

    if (setjmp(png_jmpbuf(png_ptr))) {
        // Free all of the memory associated with the png_ptr and info_ptr
        LOG(Debugability::DBG_LVL_HIGH, "Image '%s' is not a well formed PNG image.", m_filename.c_str());
        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
        // If we get here, we had a problem reading the file
        return false;
    }

    // Set up the input control if you are using standard C streams
    png_init_io(png_ptr, m_fp);

    // If we have already read some of the signature
    png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);

    /*
    * If you have enough memory to read in the entire image at once,
    * and you need to specify only transforms that can be controlled
    * with one of the PNG_TRANSFORM_* bits (this presently excludes
    * dithering, filling, setting background, and doing gamma
    * adjustment), then you can read the entire image (including
    * pixels) into the info structure with this call:
    */
    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, nullptr);

    /* At this point you have read the entire image */
    m_width = static_cast<int32>(png_ptr->width);
    m_height = static_cast<int32>(png_ptr->height);
    png_bytep* rowPointers = png_get_rows(png_ptr, info_ptr);
        
    switch(info_ptr->color_type) {
        case PNG_COLOR_TYPE_PALETTE:
            buffer.resize(m_width * m_height * 3);
            for(int32 row = 0; row < m_height; ++row) {
                for(int32 col = 0; col < m_width; ++col) {
                    buffer[3 * (row * m_width + col) + 0] =
                            info_ptr->palette[rowPointers[m_height - row - 1][col]].red;
                    buffer[3 * (row * m_width + col) + 1] =
                            info_ptr->palette[rowPointers[m_height - row - 1][col]].green;
                    buffer[3 * (row * m_width + col) + 2] =
                            info_ptr->palette[rowPointers[m_height - row - 1][col]].blue;
                }
            }
            m_component_type = GL_RGB;
            break;
        case PNG_COLOR_TYPE_RGB:
            buffer.resize(m_width * m_height * 3);
            for(int32 row = 0; row < m_height; ++row) {
                memcpy(&buffer[3 * row * m_width], rowPointers[m_height - row - 1], m_width * 3);
            }
            m_component_type = GL_RGB;
            break;
        case PNG_COLOR_TYPE_RGB_ALPHA:
            buffer.resize(m_width * m_height * 4);
            for(int32 row = 0; row < m_height; ++row) {
                memcpy(&buffer[4 * row * m_width], rowPointers[m_height - row - 1], m_width * 4);
            }
            m_component_type = GL_RGBA;
            break;

        default:
            // the image format isn't supported.
            LOG(Debugability::DBG_LVL_HIGH, "PNG image format %d not supported for image '%s'.",
                    info_ptr->color_type, m_filename.c_str());
            return false;
    }

    // clean up after the read, and free any memory allocated - REQUIRED
    png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

    return true;
}

bool PNGReader::CheckIfPNG() {
    char buf[PNG_BYTES_TO_CHECK];
    // Check to ensure the PNG file was opened.
    if (nullptr == m_fp) {
        return false;
    }

    // Read in some of the signature bytes
    if (fread(buf, 1, PNG_BYTES_TO_CHECK, m_fp) != PNG_BYTES_TO_CHECK) {
        return false;
    }

    // Compare the first PNG_BYTES_TO_CHECK bytes of the signature.
    //    Return nonzero (true) if they match
    if (png_sig_cmp((png_bytep)buf, (png_size_t)0, PNG_BYTES_TO_CHECK)) {
        return false;
    }
    return true;
}

}
}