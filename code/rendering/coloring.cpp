#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "coloring.hpp"

#define CONST_PI 3.141592653589793238l

using namespace std;

png::rgb_pixel compute_color(
    const int& iter, const int& max_iter, const long double& last_re, const long double& last_im,
    const long double& bailout, const coloring_mode& cmode, const vector<png::rgb_pixel>& palette){

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
    switch(cmode){
        case coloring_mode::linear:         //LINEAR
            calculated_index = (normalized_iter) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case coloring_mode::ln:             //SMOOTH / LOG
        {
            if(iter == max_iter) return png::rgb_pixel(0, 0, 0);
            const long double smooth_iter = (long double)iter - log(log(last_re * last_re + last_im * last_im)/(2.0l * log(bailout))) / log(2.0l);
            if(smooth_iter < 0.0l) return palette.front();
            calculated_index = fmod(smooth_iter, palette.size());
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = (size_t)ceil(calculated_index) % palette.size();
        }
            break;

        case coloring_mode::sqrt:           //SQRT
            calculated_index = sqrt(normalized_iter) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case coloring_mode::cbrt:           //CBRT
            calculated_index = cbrt(normalized_iter) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case coloring_mode::scurve:         //S CURVE
            calculated_index = ((normalized_iter >= 0.5 ? powl((normalized_iter*2 - 2), 3) + 2 : powl(normalized_iter*2, 3)) / 2) * (palette.size() - 1);
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case coloring_mode::lastangle:      //ANGLE
            calculated_index = atan2(last_im, last_re);     //Return range [-pi, +pi]
            calculated_index = calculated_index/(2.0l * CONST_PI) + 0.5l;
            color_index = floor(calculated_index);
            round_error = calculated_index - color_index;
            interpolation_index = ceil(calculated_index);
            break;

        case coloring_mode::binary:         //BINARY / DEFAULT
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

png::rgb_pixel invert_color(const png::rgb_pixel& c){
    return png::rgb_pixel(255 - c.red, 255 - c.green, 255 - c.blue);
}