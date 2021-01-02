#ifndef MISC_FUNCTIONS_HPP_INCLUDED
#define MISC_FUNCTIONS_HPP_INCLUDED

#include <iostream>                 //For debugging purposes
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <png++/png.hpp>

#define CONST_PI 3.141592653589793238l

using namespace std;

///Function to compute the color of a pixel in the image
//Coloring mode
/*Supported coloring modes
* 0  -> binary (black or white)
* 1  -> ln
* 2  -> square root of iteration
* 3  -> cubic root of iteration
* 4  -> angle of the last iterated point (calculated through arctan)

TODO
* FIX LN
* 5  -> bailout / norm of last iter
* 6  -> imaginary part of the last iterated point (normalized to idk what)
* 7  -> real + imag (normalized to idk what)
* 8  -> real * imag (normalized to idk what)
* 9  -> linear
* Y  -> ln(ln()) but make it option number 2 and shift all the others
*/
png::rgb_pixel compute_color(
    const int iter, const int max_iter, const long double last_re, const long double last_im, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette){

    if(iter > max_iter || iter < 0 || max_iter < 0) return png::rgb_pixel(0, 0, 0);

    //Variables utilized for all the different coloring techniques
    //Normalized iteration
    const long double normalized_iter = (long double)iter/(long double)max_iter;
    /*These three indices are used to indicate the two closest (by index) colors in the palette to the calculated correct (probably fractional) index
    * - calculated_index represents what the calculated color should be in the palette. It's usually computed by mapping something which ranges from
    *   [0,1] to [0, palette.size()-1]. Since this variable will probably never be an integer, we need the two other indices
    * - color_index is floor(calculated_index). Represents the starting color
    * - interpolation_index is ceil(calculated_index). Represents the other color for interpolation
    *   We add to calculated color simply round_error * deltaColor between the two
    */
    long double calculated_index = 0;
    size_t color_index = 0;
    size_t interpolation_index = 0;
    //This variable is used to measure the distance between the calculated_index and color_index
    long double round_error = 0;
    //Returned variable
    png::rgb_pixel ret_color(0, 0, 0);

    //Switching color modes
    switch(color_mode){
        case 1:         //SMOOTH / LOG
        {
            if(iter == max_iter) return png::rgb_pixel(0, 0, 0);
            const long double smooth_iter = (long double)iter - log(log(last_re * last_re + last_im * last_im)/(2.0l * log(bailout))) / log(2.0l);
            if(smooth_iter < 0.0l) return palette.front();
            calculated_index = fmod(smooth_iter, (palette.size() - 1));
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
        }
            break;

        case 2:         //SQRT
            calculated_index = sqrt(normalized_iter) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case 3:         //CBRT
            calculated_index = cbrt(normalized_iter) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case 4:         //ANGLE
            calculated_index = atan2(last_im, last_re);     //Return range [-pi, +pi]
            calculated_index = calculated_index/(2.0l * CONST_PI) + 0.5l;
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case 0:         //BINARY / DEFAULT
        default:
            if(iter == max_iter)
                return png::rgb_pixel(255, 255, 255);
            else
                return png::rgb_pixel(0, 0, 0);
            break;
    }

    //Check if the interpolation index is out of bounds
    if(interpolation_index == palette.size()) --interpolation_index;

    //Assign rounded value to ret_color
    ret_color = palette[color_index];
    //Apply correction with linear interpolation
    ret_color.red   += round_error * (palette[interpolation_index].red   - palette[color_index].red);       //round_error * deltaRed
    ret_color.green += round_error * (palette[interpolation_index].green - palette[color_index].green);     //round_error * deltaGreen
    ret_color.blue  += round_error * (palette[interpolation_index].blue  - palette[color_index].blue);      //round_error * deltaBlue

    return ret_color;
}


///Function to load the color palette from file
//Returns 0 if the loading was succesful, 1 if some errors occurred
/*The config file of the color palette should have the following format
* uint8 uint8 uint8
* uint8 uint8 uint8
* uint8 uint8 uint8
* [...]
*
* Each line represents a color to be loaded in the color palette, the first row corresponds
* to the RED channel, the middle row to the GREEN channel and the last row to the BLUE channel
*/
int load_palette(vector<png::rgb_pixel>& palette, string file_name){
    //Clear color palette
    palette.clear();

    //Open input text file
    ifstream input_file(file_name, ifstream::in);
    if(!input_file.is_open()) return 1;

    //Buffer to read the file line by line
    string line;
    while(getline(input_file, line)){
        //String stream to parse the line
        istringstream iss(line);
        //Temporary integers to store the values of the color, because the character 0 in the file gets interpreted as color 48,
        //so we have to store the reading from the file in an integer, so that it gets interpreted as a number
        int tmpred = 0, tmpgreen = 0, tmpblue = 0;

        //Read from iss and check that we read 3 color values. If not, return 1
        if(!(iss >> tmpred >> tmpgreen >> tmpblue)) return 1;

        palette.push_back(png::rgb_pixel(uint8_t(tmpred), uint8_t(tmpgreen), uint8_t(tmpblue)));
    }

    return 0;
}

///Function that prints info about the program and its usage
//Print help of the program, usage and flags
void printHelp(char *argv[]){
    cout << "Fractal generator by git-gabri, v0.7" << endl;
    cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
    cout << R"foo(
The name of this program was inspired by a much more popular fractal rendering
program called "Kalles Fraktaler".

HELP
    -H          show this help

IMAGE RELATED FLAGS
    -w SIZE_T   sets width of the image.  The default is 1920
    -h SIZE_T   sets height of the image. The default is 1080

    -S SIZE_T   sets maximum sector side.
                To utilize multiple processors to speed up
                the rendering this program subdivides the image in sectors, whose
                dimensions are at maximum 256x256 pixels, by default. This flag
                allows the user to change the maximum side of the sectors.

    -o STRING   sets output file name to STRING. Do not include the extension, the
                output image is always a PNG file.
                The default name is "fractal.png"

FRACTAL RELATED FLAGS
    -f INT      sets fractal type. Must be a number in [0, 3].
                0 -> Mandelbrot set
                1 -> Burning ship
                2 -> Mandelbar
                3 -> Spade fractal (z^z + c)
                The default is 0

    -r LNG_DBL  sets real part of the offset from the origin
    -i LNG_DBL  sets imag part of the offset from the origin
                The rendering of the fractal is done around a point whose coordinates
                are specified by these two flags. The default is the origin (0, 0);

    -j          enable Julia mode
    -A LNG_DBL  sets real part of the constant c in Julia mode
    -B LNG_DBL  sets imag part of the constant c in Julia mode
                In Julia mode the constant c in the iterative formulas to render the
                different fractal is fixed, and doesn't depend on the starting point
                in the complex plane. These three flags allow the user to run the
                rendering in Julia mode and to set the real and imaginary part of this
                constant. The default c is the origin (0, 0);

    -s LNG_DBL  sets scaling factor.
                Changes the amount of zoom done around the centeral rendering point.
                The default is 1

    -t INT      sets maximum number of iterations. The default is 2000

    -b LNG_DBL  sets bailout radius for the fractal. The default is 2

    -T INT      sets maximum number of threads working on the rendering of the image.
                The default is 8

COLOR RELATED FLAGS
    -p STRING   sets filename of input config file for color palette to STRING.
                This means that when this flag is used, the program looks in the folder
                into which the executable is placed, for a file whose name is the
                content of STRING. It must be a simple .txt file.
                Each line in this file represents a color that is then loaded in the
                color palette. The colors must be specified by writing their
                individual R(ed), G(reen) and B(lue) channel values separated by a
                space.
                For example, a line containing "255 255 255" represents white,
                "0 0 0" represents black and "255 255 0" is yellow.
                In general,
                the format for a single line is "uint8_t uint8_t uint8_t".
                By default, or if an error is found in the config file, the palette
                is loaded with just black and white.

    -c INT      sets coloring mode. Must be a number in [0, 4].
                The supported coloring modes are:
                0  -> binary (black or white)
                1  -> ln (which is the most popular option but I'm really struggling in
                          implementing it properly. It doesn't look that great)
                2  -> square root of iteration
                3  -> cubic root of iteration
                4  -> angle of the last iterated point (calculated through arctan)
                The default is 3, because it's the best looking one.

OTHER FLAGS
    -v          print verbose output, more user friendly
)foo";
    return;
}

#endif // MISC_FUNCTIONS_HPP_INCLUDED
