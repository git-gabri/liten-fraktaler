#ifndef CALC_AND_COLOR_HPP_INCLUDED
#define CALC_AND_COLOR_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <functional>
#include <png++/png.hpp>

#include "threadpool.hpp"
#include <future>

#include <structs.hpp>
#include "fractals.hpp"

#define vcout if(consettings.verbose_output) cout        //print only if verbose output enabled (verbose-cout)

using namespace std;

//printinfo prototype
void print_info( const imagesettings_t& isettings,
                const fractalsettings_t& fsettings,
                const colorsettings_t& csettings,
                const rendersettings_t& rsettings,
                const consolesettings_t& consettings);

//Functions
//main function to compute the entire fractal image
png::image<png::rgb_pixel> launch_render(  imagesettings_t& isettings,
                                            fractalsettings_t& fsettings,
                                            colorsettings_t& csettings,
                                            rendersettings_t& rsettings,
                                            consolesettings_t& consettings,
                                            const vector<png::rgb_pixel> palette){

    //To better distribute the workload between all the threads, the image gets divided into sectors, then
    //when one thread working on a sector finishes, we launch another one to work on another sector, so that we
    //(almost) always have maxThreads working on a sector
    vector<array<size_t, 4>> sectors = {};
    //Image of the fractal
    png::image<png::rgb_pixel> fractal_image(isettings.image_width, isettings.image_height);

    //Creating sectors onto which thread can work in parallel
    for(size_t i = 0; i < isettings.image_height; i += rsettings.max_sector_size) {
        for(size_t j = 0; j < isettings.image_width; j += rsettings.max_sector_size) {
            size_t start_x = j;
            size_t start_y = i;
            size_t end_x   = j + min((size_t)rsettings.max_sector_size, isettings.image_width - j);
            size_t end_y   = i + min((size_t)rsettings.max_sector_size, isettings.image_height - i);

            sectors.push_back({start_x, start_y, end_x, end_y});
        }
    }
    size_t totalSectors = sectors.size();

    /*
    //DEBUG STUFF
    size_t totalAreaSectors = 0;
    cout << "Image size: " << image_width << " x " << image_height << endl;
    cout << "# of sectors: " << totalSectors << endl;
    for(int i = 0; i < totalSectors; ++i){
        cout << "Sector #" << i << ":" << endl;
        cout << "   Start : (" << sectors[i][0] << ", " << sectors[i][1] << ")" << endl;
        cout << "   Finish: (" << sectors[i][2] << ", " << sectors[i][3] << ")" << endl;
        size_t area = (sectors[i][3] - sectors[i][1]) * (sectors[i][2] - sectors[i][0]);
        cout << "   Area  :  " << area << " pixels" << endl;
        totalAreaSectors += area;
    }
    cout << "Total area: " << totalAreaSectors << "    Image area: " << image_width*image_height << endl;
    return 0;
    */

    //Function pointer to the fractal function to call
    void (*fractal_fn_pointer)(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write);
    //Initialize function pointer based on selected fractal
    switch(fsettings.fractal_type){
        case ftype::mandelbrot:     //MANDELBROT
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &mandelbrot;
            break;

        case ftype::tippets:        //TIPPETS MANDELBROT
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &tippets_mandelbrot;
            break;

        case ftype::burnship:       //BURNING SHIP
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &burning_ship;
            break;

        case ftype::mandelbar:      //MANDELBAR
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &mandelbar;
            break;

        case ftype::spade:          //SPADE FRACTAL
            if(fsettings.bailout < 0) fsettings.bailout = 4;
            fractal_fn_pointer = &spadefract;
            break;

        case ftype::magnet1:        //MAGNET 1
            if(fsettings.bailout < 0) fsettings.bailout = 128;
            fractal_fn_pointer = &magnet_type1;
            break;

        case ftype::magnet2:        //MAGNET 2
            if(fsettings.bailout < 0) fsettings.bailout = 1024;
            fractal_fn_pointer = &magnet_type2;
            break;

        case ftype::cactus:         //CACTUS
            if(fsettings.bailout < 0) fsettings.bailout = 8;
            fractal_fn_pointer = &cactus;
            break;

        case ftype::zubieta:        //ZUBIETA
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &zubieta;
            break;

        case ftype::zubitheta:      //ZUBITHETA
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &zubitheta;
            break;

        case ftype::logmap:         //LOGISTIC MAP
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &logistic_map;
            break;

        case ftype::unpolsquare:    //UNPOLISHED SQUARE
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &unpol_square;
            break;

        case ftype::moth:           //MOTH
            if(fsettings.bailout < 0) fsettings.bailout = 4;
            fractal_fn_pointer = &moth;
            break;

        case ftype::cave:           //CAVE
            if(fsettings.bailout < 0) fsettings.bailout = 16;
            fractal_fn_pointer = &cave;
            break;

        case ftype::wankel:         //WANKEL
            if(fsettings.bailout < 0) fsettings.bailout = 4;
            fractal_fn_pointer = &wankel;
            break;

        case ftype::seaangel:            //SEA ANGEL
            if(fsettings.bailout < 0) fsettings.bailout = 64;
            fractal_fn_pointer = &sea_angel;
            break;

        case ftype::test1:           //TESTING FORMULA
            if(fsettings.bailout < 0) fsettings.bailout = 4;
            fractal_fn_pointer = &test;
            break;

        case ftype::test2:
            if(fsettings.bailout < 0) fsettings.bailout = 16;
            fractal_fn_pointer = &test2;
            break;

        default:            //DEFAULT = MANDELBROT
            fractal_fn_pointer = &mandelbrot;
            break;
    }

    //Print info if required
    if(consettings.verbose_output) print_info(isettings, fsettings, csettings, rsettings, consettings);

    //Create threadpool for rendering
    threadpool renderpool(rsettings.max_threads);
    //vector of void to ensure the main threads proceeds only after all the sectors have rendered
    vector<future<void>> completed_sectors;

    //For every sector
    for(auto s : sectors){
        size_t start_x = s[0];
        size_t start_y = s[1];
        size_t end_x   = s[2];
        size_t end_y   = s[3];

        //Enqueue a job to the renderpool
        completed_sectors.emplace_back(
            renderpool.enqueue(
                fractal_fn_pointer,          //Fractal to launch
                start_x, start_y,            //(x,y) starting position
                end_x, end_y,                //(x,y) ending position
                fsettings,                   //fractal settings
                csettings,                   //image settings
                ref(palette),                //Reference to color palette
                ref(fractal_image)           //Reference to image to update pixels
            )
        );
    }

    //Once all the jobs are enqueued, wait for all of them to finish
    for(size_t i = 0; i < completed_sectors.size(); ++i){
        completed_sectors[i].get();
        vcout << "Completed sectors: " << i << "/" << totalSectors << "\r" << flush;
    }

    //If specified, draw the crosshair on the output image
    if(csettings.draw_crosshair){
        const size_t halfWidth = isettings.image_width / 2;
        const bool evenWidth = (isettings.image_width % 2 == 0);
        const size_t halfHeight = isettings.image_height / 2;
        const bool evenHeight = (isettings.image_height % 2 == 0);

        //Draw vertical center line
        for(size_t i = 0; i < isettings.image_height; ++i){
            fractal_image[i][halfWidth] = invert_color(fractal_image[i][halfWidth]);
            //If the image has an even number of pixels in width, make the line 2 pixels thick
            if(evenWidth)
                fractal_image[i][halfWidth - 1] = invert_color(fractal_image[i][halfWidth - 1]);
        }

        //Draw horizontal center line
        for(size_t i = 0; i < isettings.image_width; ++i){
            fractal_image[halfHeight][i] = invert_color(fractal_image[halfHeight][i]);
            //If the image has an even number of pixels in height, make the line 2 pixels thick
            if(evenHeight)
                fractal_image[halfHeight - 1][i] = invert_color(fractal_image[halfHeight - 1][i]);
        }
    }

    return fractal_image;
}

//----------------------------------------------------------------------------------------------------------------------------
void print_info( const imagesettings_t& isettings,
                const fractalsettings_t& fsettings,
                const colorsettings_t& csettings,
                const rendersettings_t& rsettings,
                const consolesettings_t& consettings
              ){

    string fractal_name = ftype_to_string(fsettings.fractal_type);

    const string coloring_mode_name = coloring_mode_to_strign(csettings.cmode);

    cout << "Rendering      : " << fractal_name << (fsettings.julia_mode ? " (J)" : "") << " (" << rsettings.max_threads << " thread(s))" << endl;
    cout << "Image size     : " << isettings.image_width << " x " << isettings.image_height << endl;
    cout << "Sectors up to  : " << rsettings.max_sector_size << "x" << rsettings.max_sector_size << " = " << rsettings.max_sector_size * rsettings.max_sector_size << " pixels^2" << endl;

    cout << "Iterations     : " << fsettings.max_iter << endl;
    cout << "Bailout radius : " << fsettings.bailout << endl;
    cout << "Coloring mode  : " << coloring_mode_name << endl;
    cout << "Centered at    : (" << fsettings.offset_re << ", " << fsettings.offset_im << " i)" << endl;
    if(fsettings.julia_mode) cout << "Julia c const. : (" << fsettings.julia_re << ", " << fsettings.julia_im << " i)" << endl;

    //Auxiliary variable, used when calculating scaling from pixel space to complex plane, to prevent stretching
    const auto squareScale = min(isettings.image_width, isettings.image_height);
    //Top left corner
    const long double start_re = 4 * (- (long double)isettings.image_width / 2)  / ((long double)squareScale * fsettings.scaling_factor) + fsettings.offset_re;
    const long double start_im = 4 * (- (long double)isettings.image_height / 2) / ((long double)squareScale * fsettings.scaling_factor) + fsettings.offset_im;
    //Bottom right
    const long double stop_re  = 4 * ((long double)isettings.image_width / 2)  / ((long double)squareScale * fsettings.scaling_factor) + fsettings.offset_re;
    const long double stop_im  = 4 * ((long double)isettings.image_height / 2) / ((long double)squareScale * fsettings.scaling_factor) + fsettings.offset_im;
    cout << "Span  : (" << stop_re - start_re << " x " << stop_im - start_im << ")" << endl;
    cout << "Scale : " << fsettings.scaling_factor << endl;
}

#endif // CALC_AND_COLOR_HPP_INCLUDED
