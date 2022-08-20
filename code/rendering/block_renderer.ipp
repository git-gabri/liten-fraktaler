#ifndef BLOCK_RENDERER_IPP_INCLUDED
#define BLOCK_RENDERER_IPP_INCLUDED

#include "liten.hpp"
#include "fractals.hpp"
#include "coloring.hpp"

#include <algorithm>

//Main function that runs in a thread whose job is to render the colors in a specific block in the image
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto square_scale = std::min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout = fsettings.bailout;
    const long double bailout_2 = bailout * bailout;

    std::complex<long double> z;
    std::complex<long double> c;
    bool stop_iterating = false;
    std::vector<std::complex<long double>> history;
    std::vector<std::complex<long double>> extra_params;

    //Calculate the values of the pixels for each point in the assigned sector of the image
    for(size_t y = startY; y < endY; y++) {
        for(size_t x = startX; x < endX; x++) {
            //Initialize fractal parameters z, c, etc...
            init_fractal(x, y, width, height, square_scale, z, c, history, extra_params, fsettings);

            //Auxiliary variables to count the number of iterations and when to stop iterating
            size_t iter_count = (FIRST_ITERATION_IS_0 ? -1 : 0);
            stop_iterating = false;
            do {
                //compute new value for fractal
                z = (*fractal_func)(z, c, history, extra_params, fsettings);

                if(std::norm(z) > bailout_2)
                    stop_iterating = true;

                ++iter_count;
            } while((iter_count < fsettings.max_iter) && !stop_iterating);

            image_to_write[y][x] = compute_color(iter_count, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.cmode, palette);
        }
    }

    return;
}

#endif