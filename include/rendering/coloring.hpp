#include "structs.hpp"

#include <string>
#include <vector>
#include <png++/png.hpp>

//Function to compute the color of a pixel in the image
//Coloring mode
/*Supported coloring modes
* 0  -> binary (black or white)
* 1  -> linear
* 2  -> ln
* 3  -> square root of iteration
* 4  -> cubic root of iteration
* 5  -> s curve, made with x^3 from 0 to 1 and (x-2)^3 + 2 from 1 to 2
* 6  -> angle of the last iterated point (calculated through arctan)

TODO
* FIX LN
* N  -> bailout / norm of last iter
* N  -> imaginary part of the last iterated point (normalized to idk what)
* N  -> real + imag (normalized to idk what)
* N  -> real * imag (normalized to idk what)
* N  -> ln(ln()) but make it option number 2 and shift all the others
*/
png::rgb_pixel compute_color(
    const int& iter, const int& max_iter, const long double& last_re, const long double& last_im,
    const long double& bailout, const coloring_mode& cmode, const std::vector<png::rgb_pixel>& palette);

//Function to invert a single color
png::rgb_pixel invert_color(const png::rgb_pixel& c);

//Function to load the color palette from file
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
int load_palette(std::vector<png::rgb_pixel>& palette, std::string file_name);