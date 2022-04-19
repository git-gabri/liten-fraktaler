#ifndef CALC_AND_COLOR_HPP_INCLUDED
#define CALC_AND_COLOR_HPP_INCLUDED

#include "structs.hpp"

#include <vector>
#include <complex>
#include <png++/png.hpp>

//Type alias for the pointer of the function implementing the fractal iteration
using fractal_fn_ptr_t = 
    std::complex<long double> (*)(const std::complex<long double>& last_z,
                                  const std::complex<long double>& c,
                                  const std::vector<std::complex<long double>>& history,
                                  const std::vector<std::complex<long double>>& extra_params,
                                  const fractalsettings_t& fsettings);

void block_renderer(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t& fsettings, const colorsettings_t& csettings,
    const std::vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write);

long double init_bailout_radius(const ftype& f);
fractal_fn_ptr_t init_fractal_fn_pointer(const fractalsettings_t& fsettings);

//Functions
//Function that makes multiple threads to render the image
png::image<png::rgb_pixel> launch_render(   imagesettings_t& isettings,
                                            fractalsettings_t& fsettings,
                                            colorsettings_t& csettings,
                                            rendersettings_t& rsettings,
                                            consolesettings_t& consettings,
                                            const std::vector<png::rgb_pixel> palette);

#endif // CALC_AND_COLOR_HPP_INCLUDED