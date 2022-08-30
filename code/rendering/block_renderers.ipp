#ifndef BLOCK_RENDERER_IPP_INCLUDED
#define BLOCK_RENDERER_IPP_INCLUDED

#include "liten.hpp"
#include "fractals.hpp"
#include "coloring.hpp"

#include <algorithm>
#include <cmath>

//--------------------------------------------------------------------------------------------------------------------------------------------------
// REGULAR BLOCK RENDERERS
//--------------------------------------------------------------------------------------------------------------------------------------------------
//Main functions that run in threads whose job is to render the colors in a specific block in the image

//Basic version:
//iterate over every pixel in the block, iterate the fractal function, check after every iteration
template<fractal_fn_ptr_t fractal_func>
void lf::internals::basic_block_renderer(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
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

//MIBC (Multiple Iterations Between Checks) version:
//iterate over every pixel in the block, iterate the fractal function,
//checks every 2^x iterations, if we overshoot rollback and make 2^(x-1) iterations and continue until 1 iteration is done. After that, exit
template<fractal_fn_ptr_t fractal_func>
void lf::internals::mibc_block_renderer(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
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
            size_t total_iter_count = 0;    //Total iteration count
            size_t iter_count = 0;          //Iteration count between a check and another. Used to not overshoot the maximum number of iterations

            //Number of iterations of the function to do between a check of the norm and another one
            size_t num_iter_between_checks = rsettings.mibc_max_iter_between_checks;

            stop_iterating = false;

            do {
                //compute new value for fractal skipping checks
                const std::complex<long double> prev_z = z;
                for(iter_count = 0; iter_count < num_iter_between_checks && iter_count < (fsettings.max_iter - total_iter_count); ++iter_count)
                    z = (*fractal_func)(z, c, history, extra_params, fsettings);

                //If we didn't overshoot we increment the value of total_iter_count by the amount of iterations performed
                if(std::isfinite(z.real()) && std::isfinite(z.imag()) && std::norm(z) < bailout_2){
                    total_iter_count += iter_count;
                }
                //If we overshot and performed too many iterations between a check and another
                else{
                    //If we performed just an iteration, we're already at the limit, so we must exit the do-while
                    if(num_iter_between_checks == 1){
                        stop_iterating = true;
                    }
                    else {
                        z = prev_z;
                        num_iter_between_checks /= 2;
                    }
                }
            } while(total_iter_count < fsettings.max_iter && (stop_iterating == false));

            image_to_write[y][x] = compute_color(total_iter_count, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.cmode, palette);
        }
    }

    return;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
// BLOCK RENDERERS FOR TESTING
//--------------------------------------------------------------------------------------------------------------------------------------------------
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test0(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test1(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test2(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test3(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test4(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test5(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test6(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test7(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test8(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
template<fractal_fn_ptr_t fractal_func>
void lf::internals::block_renderer_test9(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write){
    lf::internals::basic_block_renderer<fractal_func>(startX, startY, endX, endY, image_to_write);
}
#endif