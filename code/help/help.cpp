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
    This scripting language is used to specify the mathematical operations to perform in a single iteration of the fractal.
    It's interpreted in a small virtual register machine (RM), with 16 registers (0 to 15) and a read-only vector of constants.
    All of the constants and the registers are complex numbers of long doubles.

    A script contains one "instruction" per line and each line can be of 2 types:
    - constant initialization
    - instruction for the RM
    When loading the script, the constants are stored in a dedicated vector and the instructions in another one.
    When executing the script, only the instruction vector is iterated until the end of execution.

    NOTE:
    Before writing code, it's HIGHLY RECOMMENDED to READ all of the text below, especially the notes on the
    machine operation.

SCRIPTING SYNTAX
    CONSTANT INITIALIZATION
    The syntax of a line for loading a constant is the following:

    c<dest_const> = (<real>,<imag>)

    and it has to match the following regex: ^c\d+ = \S+$
    Example: "c3 = (15,1)" loads in constant 3 the value 15+i
    If in the code the constant cN is loaded with a value, all of the other constants from 0 to N-1
    will be created and default-initialized if not already present.
    If a constant is initialized more than once, only the last line initializing it will have effect.
    Constants are READ-ONLY during the script execution, the machine has no instructions to write to them.

    RM INSTRUCTIONS
    The machine can perform a certain set of operations, some are binary, in the sense that they
    take 2 arguments, and some are unary. All of them always output one result.
    The syntax of all the instructions is the following:
    
    Unary:      "<dest_reg> = <op> <src_reg>"
    Binary:     "<dest_reg> = <src1_reg> <op> <src2_reg>"

    Below is a complete list of all the operations that the RM can perform:
        *) addition         ^\d{1,2} = \d{1,2} \+ \d{1,2}$
        *) subtraction      ^\d{1,2} = \d{1,2} \- \d{1,2}$
        *) multiplication   ^\d{1,2} = \d{1,2} \* \d{1,2}$
        *) division         ^\d{1,2} = \d{1,2} \/ \d{1,2}$
            Basic arithmetic operations.
            Examples: "0 = 1 + 2"   writes in reg. 0 the result of (reg. 1 + reg. 2)
                      "3 = 4 - 5"   writes in reg. 3 the result of (reg. 4 - reg. 5)
                      "6 = 7 * 8"   writes in reg. 6 the result of (reg. 7 * reg. 8)
                      "9 = 10 / 11" writes in reg. 9 the result of (reg. 10 / reg. 11)
        *) swap             ^\d{1,2} = swap \d{1,2}$
            Swaps the real and imaginary part of the specified number.
            Example:  "0 = swap 3" copies reg. 3 in reg. 0 swapping the real and imaginary parts.
        *)inverse           ^\d{1,2} = inv \d{1,2}$
        *)inverse real      ^\d{1,2} = invre \d{1,2}$
        *)inverse imag      ^\d{1,2} = invim \d{1,2}$
            Inverse of the number (1/x) or of the real/imaginary part only.
            Examples: "0 = inv 1"   writes in reg. 0 the result of 1/(reg. 1)
                      "2 = invre 3" copies reg. 3 in reg. 2 but with the real part inverted (1/real)
                      "4 = invim 5" copies reg. 5 in reg. 4 but with the imaginary part inverted (1/imag)
        *)flip sign         ^\d{1,2} = flip \d{1,2}$
        *)flip sign real    ^\d{1,2} = flipre \d{1,2}$
        *)flip sing imag    ^\d{1,2} = flipim \d{1,2}$
            Flips the sign of the number or of the real/imaginary part only.
            Examples: "0 = flip 1"   writes in reg. 0 the result of -1*(reg. 1)
                      "2 = flipre 3" copies reg. 3 in reg. 2 but with the real part multiplied by -1
                      "4 = flipim 5" copies reg. 5 in reg. 4 but with the imaginary part multiplied by -1
        *)abs. value        ^\d{1,2} = abs \d{1,2}$
        *)abs. value real   ^\d{1,2} = absre \d{1,2}$
        *)abs. value imag   ^\d{1,2} = absim \d{1,2}$
            Takes the absolute value of the number or of the real/imaginary part only.
            Examples: "0 = abs 1"   writes in reg. 0 the result of |reg. 1|
                          Note: this operation is very different from the other 2
                      "2 = absre 3" copies reg. 3 in reg. 2 but with the absolute value of the real part
                      "4 = absim 5" copies reg. 5 in reg. 4 but with the absolute value of the imaginary part
        *)nullify real      ^\d{1,2} = nullre \d{1,2}$
        *)nullify imag      ^\d{1,2} = nullim \d{1,2}$
            Nullifies the real/imaginary part of a number.
            Examples: "0 = nullre 1" copies the imaginary part of reg. 1 in reg. 0 and sets the real part to 0
                      "2 = nullim 3" copies the real part of reg. 3 in reg. 2 and sets the imaginary part to 0
        *)sine              ^\d{1,2} = sin \d{1,2}$
        *)cosine            ^\d{1,2} = cos \d{1,2}$
        *)tangent           ^\d{1,2} = tan \d{1,2}$
            Trigonometric functions.
            Examples: "0 = sin 1" writes in reg. 0 the result of sin(reg. 1)
                      "2 = cos 3" writes in reg. 2 the result of cos(reg. 3)
                      "4 = tan 5" writes in reg. 4 the result of tan(reg. 5)
        *)hyper. sine       ^\d{1,2} = sinh \d{1,2}$
        *)hyper. cosine     ^\d{1,2} = cosh \d{1,2}$
        *)hyper. tangent    ^\d{1,2} = tanh \d{1,2}$
            Hyperbolic functions.
            Examples: "0 = sinh 1" writes in reg. 0 the result of sinh(reg. 1)
                      "2 = cosh 3" writes in reg. 2 the result of cosh(reg. 3)
                      "4 = tanh 5" writes in reg. 4 the result of tanh(reg. 5)
        *)generic power     ^\d{1,2} = \d{1,2} pow \d{1,2}$
        *)square            ^\d{1,2} = pow2 \d{1,2}$
        *)cube              ^\d{1,2} = pow3 \d{1,2}$
            Power functions.
            Examples: "0 = 1 pow 2" writes in reg. 0 the result of (reg. 1)^(reg. 2)
                      "3 = pow2 4" writes in reg. 3 the result of (reg. 4)*(reg. 4)
                      "5 = pow3 6" writes in reg. 5 the result of (reg. 6)*(reg. 6)*(reg. 6)
        *)natural log       ^\d{1,2} = log \d{1,2}$
        *)base 2 log        ^\d{1,2} = log2 \d{1,2}$
            Logarithm functions.
            Examples: "0 = log 1" writes in reg. 0 the result of log_e(reg. 1)
                      "2 = log2 3" writes in reg. 2 the result of log_2(reg. 3)
        *)copy constant to register                 ^\d{1,2} = c\d+$
        *)copy constant real to register            ^\d{1,2} = c\d+re$
        *)copy constant imag to register            ^\d{1,2} = c\d+im$
            Copying constants to registers.
            Examples: "0 = c0" copies the constant 0 in reg. 0
                      "1 = c1re" copies the real part of constant 1 in reg. 1
                      "2 = c2im" copies the imaginary part of constant 2 in reg. 2
        *)copying register to register              ^\d{1,2} = \d{1,2}$
        *)copying register real to register         ^\d{1,2} = \d{1,2}re$
        *)copying register imag to register         ^\d{1,2} = \d{1,2}im$
            Examples: "0 = 1" copies reg. 1 into reg. 0
                      "2 = 3re" copies the real part of reg. 3 into reg. 2
                      "4 = 5im" copies the imaginary part of reg. 5 into reg. 4
                          Note: the last 2 operation leave the other part of the complex number as is.
    
IMPORTANT NOTES ON RM OPERATION
    - (almost) NO BOUNDARY CHECKING IS PRESENT, if you want to access register 16 or 42, you can, but
      it's unallocated space, the program will probably crash.
      Also, if in the code loads a constant cN but then accesses cM with M > N, it's again unallocated
      space and again, the code will probably crash.
    - when rendering a fractal, the value of z is loaded in register 0, the value of c is loaded in register 1.
      The script is then executed and the new value of z is taken as what's present in register 0 at the end.
    - registers are not reset between one iteration and another one, even of different pixels.
)foo";
    return;
}