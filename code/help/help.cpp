#include "help.hpp"

using namespace std;

void print_help(ostream& os){
    os << "Fractal generator by git-gabri, v" << VERSION << "\n";
    os << "Usage: liten-fraktaler [OPTIONS]\n";
    os << R"foo(
The name of this program was inspired by a much more popular fractal rendering
program called "Kalles Fraktaler".

HELP
    -H
    --help              show this help
    -v                  
    --verbose           print verbose output during rendering, more user friendly.

IMAGE RELATED FLAGS
    -w SIZE_T
    --width SIZE_T
                        sets width of the image.  The default is 1920.
    -h SIZE_T
    --height SIZE_T
                        sets height of the image. The default is 1080.
    -o STRING   
    --output-image-filename STRING
                        sets output file name to STRING. Do not include the extension, the
                        output image is always a PNG file.
                        The default name is "fractal.png"

RENDERING RELATED FLAGS
    -S SIZE_T
    -sector-size SIZE_T 
                        sets maximum sector side.
                        To utilize multiple processors to speed up the rendering,
                        this program subdivides the image in sectors, whose dimensions are
                        at maximum 64x64 pixels, by default.
                        This flag allows the user to change the maximum side of the sectors.
    -T SIZE_T
    --max-threads SIZE_T
                        sets maximum number of threads working on the rendering of the image.
                        By default, the program tries to automatically detect the number of
                        available threads on the host machine and uses all of them. If this
                        detection fails, the maximum number of threads is set to 1.
    --renderer NAME
                        sets the renderer to use for the current image.
                        The default is "basic". The recognized names are
                        "basic"         -> classic renderer, iterates the fractal function for
                                           each pixel, checks for exit condition after calculating
                                           every iteration
                        "mibc"          -> Multiple Iteration Between Checks, works similarly to
                                           the "basic" renderer but instead of checking for the
                                           exit condition after every iteration, it performs 2^x
                                           iterations then checks. If the exit condition is not satisfied,
                                           it tries to perform 2^x more iterations. If it's satisfied
                                           it rolls back befor the initial 2^x iterations were performed
                                           and tries to iterate for 2^(x-1). And so on util only 1
                                           iteration is performed.
                                           Works well when the maximum number of iterations is high and the
                                           points in the image take a lot of time to diverge.

    --mibc-max-iter SIZE_T
                        sets the number of iterations initially performed by the "mibc" renderer.
                        If the renderer is not set to "mibc" this has no effect.
                        The number is internally rounded to the closest next power of 2 that is >= than
                        the given argument.

FRACTAL RELATED FLAGS
    -f NAME
    --fractal NAME
                        sets fractal type. The default is "mandelbrot". The recognized names are
                        "mandelbrot"    -> the classic Mandelbrot set
                        "tippets"       -> Tippets M. set (M. set implementation w bug)
                        "burnship"      -> Burning Ship
                        "mandelbar"     -> Mandelbar / Tricorn
                        "magnet1"       -> Magnet type I
                        "magnet2"       -> Magnet type II
                        "cactus"        -> Cactus fractal (z^3 + (c-1)*z - c)
                        "cactus0"       -> Cactus with z = 0 starting condition
                        "zubieta"       -> Zubieta (z^2 + c/z, fancy Julia sets)
                        "zubitheta"     -> Zubitheta (z^2 + z/c)
                        "logmap"        -> Logistic map fractal (c*z*(1-z))
                        "unpolsquare"   -> Unpolished square
                        "moth"          -> Moth fractal
                        "cave"          -> Cave fractal
                        "wankel"        -> Wankel engine fractal
                        "seaangel"      -> Sea angel fractal
                        "smith"         -> Smith Chart fractal
                        "spade"         -> Spade fractal
                        "__testX"       -> X in [0; 9]. These are used for testing and development.
                        NOTE: not all the listed fractals will render well on default settings.
    -r LNG_DBL
    --real LNG_DBL
                        sets real part of the offset from the origin.
    -i LNG_DBL  
    --imag LNG_DBL
                        sets imag part of the offset from the origin.
                        The rendering of the fractal is done around a point whose coordinates
                        are specified by these two flags. The default is the origin (0, 0);
    -s LNG_DBL 
    --scaling-factor LNG_DBL
    --zoom LNG_DBL      
                        sets scaling factor.
                        Changes the amount of zoom done around the centeral rendering point.
                        The default is 1.
    -b LNG_DBL  
    --bailout LNG_DBL
                        sets bailout radius for the fractal.
                        The default varies from fractal to fractal.
    -j
    --julia-mode        enables Julia mode
                        In Julia mode the constant c in the iterative formulas to render the
                        different fractal is fixed, and doesn't depend on the starting point
                        in the complex plane. These three flags allow the user to run the
                        rendering in Julia mode and to set the real and imaginary part of this
                        constant. The default c is the origin (0, 0).
    -A LNG_DBL
    --julia-real LNG_DBL
                        sets real part of the constant c in Julia mode.
    -B LNG_DBL          
    --julia-imag LNG_DBL
                        sets imag part of the constant c in Julia mode.
    -t SIZE_T
    --max-iter SIZE_T
                        sets maximum number of iterations. The default is 2000.

COLOR RELATED FLAGS
    -p STRING
    --palette-filename STRING
                        sets filename of input config file for color palette to STRING.
                        This means that when this flag is used, the program looks in the folder
                        into which the executable is placed, for a file whose name is the
                        content of STRING. It must be a simple .txt file.
                        Each line in this file represents a color that is then loaded in the
                        color palette. The colors must be specified by writing their
                        individual R(ed), G(reen) and B(lue) channel values separated by a
                        space.
                        For example, a line containing "255 255 255" represents white,
                        "0 0 0" represents black and "255 255 0" is yellow.
                        In general, the format for a single line is "uint8_t uint8_t uint8_t".
                        By default, the program looks for a file named "palette" (note that the
                        name has no extension. If you want to load "palette.txt" you will have
                        to specify it) and if the palette file is not found or an error occurs
                        while trying to load the colors, the palette is loaded with only black
                        and white.
    -c NAME   
    --coloring-mode NAME   
                        sets coloring mode. The default is "cbrt".
                        The supported coloring modes are:
                        "binary"        -> binary coloring (black or white)
                        "linear"        -> linear coloring
                        "ln"            -> a failed attempt at smooth coloring
                        "sqrt"          -> square root of iteration
                        "cbrt"          -> cubic root of iteraion
                        "scurve"        -> s shaped curve iteration mapping
                        "lastangle"     -> angle of the last iterated point
    -C          
    --crosshair         enables crosshair.
                        Draws a crosshair in the middle of the fractal image by inverting the colors
                        along the equator and the prime meridian.
)foo";
    return;
}