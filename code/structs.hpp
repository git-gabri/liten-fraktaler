#ifndef STRUCTS_HPP_INCLUDED
#define STRUCTS_HPP_INCLUDED

#include <string>
#include <vector>
#include <complex>

//Fracyal type enum
enum class ftype{
    mandelbrot, tippets, burnship, mandelbar, magnet1, magnet2, cactus, cactus0, zubieta, zubitheta,
    logmap, unpolsquare, moth, cave, wankel, seaangel, smith, spade,
    test0, test1, test2, test3, test4, test5, test6, test7, test8, test9,
    unknown
};

//Renderer type enum
enum class rtype{
    basic, mibc,
    test0, test1, test2, test3, test4, test5, test6, test7, test8, test9,
    unknown
};

//Coloring mode enum
enum class coloring_mode{
    binary, linear, ln, sqrt, cbrt, scurve, lastangle,
    unknown
};

//Struct containing all the settings for the fractal
struct fractalsettings_t {
    ftype fractal_type;

    long double offset_re;
    long double offset_im;

    bool julia_mode;
    bool nebula_mode;

    long double julia_re;
    long double julia_im;

    long double bailout;
    long double scaling_factor;
    size_t max_iter;

    size_t max_history_length;

    fractalsettings_t(
        ftype _fractal_type = ftype::mandelbrot,
        long double _offset_re = 0, long double _offset_im = 0,
        bool _juliaMode = false, bool _nebula_mode = false,
        long double _julia_re  = 0, long double _julia_im  = 0,
        long double _bailout = -1,
        long double _scalingFactor = 1,
        int _max_iter = 2000
    ) :
    fractal_type(_fractal_type),
    offset_re(_offset_re), offset_im(_offset_im),
    julia_mode(_juliaMode), nebula_mode(_nebula_mode),
    julia_re(_julia_re), julia_im(_julia_im),
    bailout(_bailout),
    scaling_factor(_scalingFactor),
    max_iter(_max_iter) {}
};

//Struct containing all the settings for the output image
struct imagesettings_t {
    size_t image_width;
    size_t image_height;
    std::string image_name;

    imagesettings_t(
        size_t _imageWidth = 1920,
        size_t _imageHeight = 1080,
        std::string _imageName = {"fractal"}
    ) :
    image_width(_imageWidth), image_height(_imageHeight),
    image_name(_imageName) {}
};

//Struct containing all the settings for the coloring of the image
struct colorsettings_t {
    coloring_mode cmode;
    std::string palette_filename;

    bool draw_crosshair;

    colorsettings_t(
        coloring_mode _cmode = coloring_mode::cbrt,
        std::string _palette_filename = {"palette"},
        bool _draw_crosshair = false
    ) :
    cmode(_cmode),
    palette_filename(_palette_filename),
    draw_crosshair(_draw_crosshair) {}
};

//Struct containing all the settings for the rendering of the fractal
struct rendersettings_t {
    rtype renderer_type;

    size_t max_sector_size;
    size_t max_threads;

    size_t mibc_max_iter_between_checks;

    rendersettings_t(
        const rtype& _renderer_type = rtype::basic,
        const size_t& _max_sector_size = 64,
        const size_t& _max_threads = 1,
        const size_t& _mibc_max_iter_between_checks = (size_t(1) << size_t(DEFAULT_MIBC_MAX_ITER_POW2EXP))
    ) :
    renderer_type(_renderer_type),
    max_sector_size(_max_sector_size),
    max_threads(_max_threads)
    {
        //Round up _mibc_max_iter_between_checks to next power of 2 >= than it
        size_t init_mibc_max = 1;
        while(init_mibc_max < _mibc_max_iter_between_checks && init_mibc_max != 0)
            init_mibc_max *= 2;

        //Checks if init_mibc_max got an overflow
        if(init_mibc_max == 0)
            mibc_max_iter_between_checks = (size_t(1) << size_t(DEFAULT_MIBC_MAX_ITER_POW2EXP));
        
        mibc_max_iter_between_checks = init_mibc_max;
    }
};

//Struct containing information of a single rendered pixel
//struct pixel_t{
//    unsigned char red, green, blue, alpha;
//
//    //T -> to compute
//    //C -> directly computed
//    //A -> approximated (from rectangle guessing, ecc...)
//    unsigned char status;
//
//    std::vector<std::complex<long double>> history;
//    
//    long double pickover_min_distance;
//};

//Struct containing all the settings for printing text on the console
struct consolesettings_t {
    int verbose_output;
    int colored_output;
    int suppress_warnings;
    int suppress_errors;

    consolesettings_t(
        const int& _verbose_output = 0,
        const int& _colored_output = 0,
        const int& _suppress_warnings = 0,
        const int& _suppress_errors = 0
    ) :
    verbose_output(_verbose_output),
    colored_output(_colored_output),
    suppress_warnings(_suppress_warnings),
    suppress_errors(_suppress_errors) {}
};

#endif