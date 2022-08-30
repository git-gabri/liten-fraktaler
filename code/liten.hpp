#ifndef LITEN_HPP_INCLUDED
#define LITEN_HPP_INCLUDED

#include "structs.hpp"
#include "render_defs.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <png++/png.hpp>

//This namespace is more like a singleton
namespace lf{
    //Initialization of the namespace
    void init();

    //Parse command line options
    //Implemented in:   parse_options.cpp
    int parse_options(std::vector<std::string> options);

    //Load palette from file
    //Implemented in:   load_palette.cpp
    int load_palette();

    //Render the image
    png::image<png::rgb_pixel> launch_render();

    //Anonymous namespace containing internal functions and data
    namespace internals{
        //---------------------------------------------------------------------------------------
        //Private members
        extern imagesettings_t isettings;               //struct containing various image related settings
        extern fractalsettings_t fsettings;             //struct containing various fractal related settings
        extern colorsettings_t csettings;               //struct containing various color related settings
        extern rendersettings_t rsettings;              //struct containing various information on how to render the image
        extern consolesettings_t consettings;           //struct containing various information on what the program should output on the console

        extern std::vector<png::rgb_pixel> palette;     //color palette

        //---------------------------------------------------------------------------------------
        //Private methods

        //Print various information about the current render
        void print_render_info();

        //Get function pointer to block renderer
        //Switch wrt the type of block renderer used
        block_renderer_fn_ptr_t get_block_renderer_ptr();
        //Switch wrt the type of fractal to render
        block_renderer_fn_ptr_t gbr_select_fractal();       //gbr = get_block_renderer, function called by get_block_renderer_ptr
        template<fractal_fn_ptr_t fractal_func>
        block_renderer_fn_ptr_t gbr_select_renderer();      //gbr = get_block_renderer, function called by gbr_select_fractal

        //Different block renderers used in the program
        //Basic block renderer which iterates over all pixels in a block and colors them
        template<fractal_fn_ptr_t fractal_func>
        void basic_block_renderer(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void mibc_block_renderer(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);

        //Block renderers for the different test fractals. Each of them has its own
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test0(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test1(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test2(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test3(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test4(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test5(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test6(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test7(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test8(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);
        template<fractal_fn_ptr_t fractal_func>
        void block_renderer_test9(const size_t startX, const size_t startY, const size_t endX, const size_t endY, png::image<png::rgb_pixel>& image_to_write);

        //Get default bailout radius for different fractals
        long double default_bailout_radius(const ftype& f);
    }
}

#include "block_renderers.ipp"

#endif