#ifndef STRUCTS_H_INCLUDED
#define STRUCTS_H_INCLUDED

#include <string>

struct fractalsettings_t {
    int fractal_type;

    long double offset_re;
    long double offset_im;

    bool juliaMode;
    bool buddhaMode;

    long double julia_re;
    long double julia_im;

    long double bailout;
    long double scalingFactor;
    int max_iter;

    fractalsettings_t(
        int _fractal_type = 0,
        long double _offset_re = 0, long double _offset_im = 0,
        bool _juliaMode = false, bool _buddhaMode = false,
        long double _julia_re  = 0, long double _julia_im  = 0,
        long double _bailout = -1,
        long double _scalingFactor = 1,
        int _max_iter = 2000
    ) :
    fractal_type(_fractal_type),
    offset_re(_offset_re), offset_im(_offset_im),
    juliaMode(_juliaMode), buddhaMode(_buddhaMode),
    julia_re(_julia_re), julia_im(_julia_im),
    bailout(_bailout),
    scalingFactor(_scalingFactor),
    max_iter(_max_iter) {}
};

struct imagesettings_t {
    size_t imageWidth;
    size_t imageHeight;
    size_t maxSectorSide;
    std::string imageName;

    imagesettings_t(
        size_t _imageWidth = 1920,
        size_t _imageHeight = 1080,
        size_t _maxSectorSide = 64,
        std::string _imageName = {"fractal"}
    ) :
    imageWidth(_imageWidth), imageHeight(_imageHeight),
    maxSectorSide(_maxSectorSide),
    imageName(_imageName) {}
};

struct colorsettings_t {
    int colorMode;
    std::string paletteConfig;

    colorsettings_t(
        int _colorMode = 4,
        std::string _paletteConfig = {"palette"}
    ) :
    colorMode(_colorMode),
    paletteConfig(_paletteConfig) {}
};

#endif // STRUCTS_H_INCLUDED
