#ifndef RENDER_DEFS_HPP_INCLUDED
#define RENDER_DEFS_HPP_INCLUDED

#include "structs.hpp"

#include <complex>
#include <vector>
#include <png++/png.hpp>

//Type alias for the pointer to a function implementing the fractal iteration
using fractal_fn_ptr_t = 
    std::complex<long double> (*)(const std::complex<long double>& last_z,
                                  const std::complex<long double>& c,
                                  const std::vector<std::complex<long double>>& history,
                                  const std::vector<std::complex<long double>>& extra_params,
                                  const fractalsettings_t& fsettings);

//Type alias for the pointer to a block renderer
using block_renderer_fn_ptr_t =
    void (*)(const size_t startX, const size_t startY,
             const size_t endX, const size_t endY,
             png::image<png::rgb_pixel>& image_to_write);

#endif