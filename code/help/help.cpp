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
                        Has no effect if the specified fractal type is "script".

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
                        "script"        -> Fractal defined through a simple scripting language
                        "__testX"       -> X in [0; 9]. These are used for testing and development.
                        NOTE: not all the listed fractals will render well on default settings.
    --fsf STRING
    --fractal-script-filename STRING
                        specifies the name of the file containing the script defining the custom
                        fractal.
                        More on the syntax of the scripting language in the dedicated section below.
                        Has no effect if the fractal type is not "script", see option above.
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

FRACTAL SCRIPTING LANGUAGE
    This scripting language is interpreted in a small virtual register machine, with 16 registers (0 to 15)
    and a dynamically sized vector of constants.
    All of the constants and the registers are complex numbers of long doubles. Before writing code, it's
    recommended to read all of the text below, especially the notes on the machine operation.

    The machine can perform the following operations:
        *) loading of a constant
           Syntax: "c<dest_const> <- (<real>,<imag>)"
           Regex: ^c\d+ <- \S+$
           Example: "c3 <- (15,1)" loads in constant 3 the value 15+i
        *) addition
           Syntax: "<dest_reg> <- <src1_reg> + <src2_reg>"
           Regex: ^\d{1,2} <- \d{1,2} \+ \d{1,2}$
           Example: "0 <- 1 + 2" writes in register 0 the result of (register 1 + register 2)
        *) subtraction
           Syntax: "<dest_reg> <- <src1_reg> - <src2_reg>"
           Regex: ^\d{1,2} <- \d{1,2} \- \d{1,2}$
           Example: "1 <- 5 - 6" writes in register 1 the result of (register 5 - register 6)
        *) multiplication
           Syntax: "<dest_reg> <- <src1_reg> * <src2_reg>"
           Regex: ^\d{1,2} <- \d{1,2} \* \d{1,2}$
           Example: "2 <- 4 * 8" writes in register 1 the result of (register 4 * register 8)
        *) division
           Syntax: "<dest_reg> <- <src1_reg> / <src2_reg>"
           Regex: ^\d{1,2} <- \d{1,2} \/ \d{1,2}$
           Example: "3 <- 10 / 12" writes in register 3 the result of (register 10 / register 12)
        *) copy of a constant to a register
           Syntax: "<dest_reg> <- c<src_const>"
           Regex: ^\d{1,2} <- c\d+$
           Example: "9 <- c0" writes in register 9 the constant 0
        *) copy of a register to a register
           Syntax: "<dest_reg> <- <src_reg>"
           Regex: ^\d{1,2} <- \d{1,2}+$
           Example: "14 <- 15" copies the content of register 15 into register 14
    
    Important notes on the register machine operation:
    - if in the code the constant cN is loaded with a value, all of the other constants from 0 to N-1
      will be created if not already present;
    - (almost) NO BOUNDARY CHECKING IS PRESENT, if you want to access register 16 or 42, you can, but
      it's unallocated space, the program will probably crash.
      Also, if in the code loads a constant cN but then accesses cM with M > N, it's again unallocated
      space and again, the code will probably crash.
    - when rendering a fractal, the value of z is loaded in register 0, the value of c is loaded in register 1.
      The new value of z is taken as what's present in register 0 after the script terminates.
)foo";
    return;
}