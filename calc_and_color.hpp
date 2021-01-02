#ifndef CALC_AND_COLOR_HPP_INCLUDED
#define CALC_AND_COLOR_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <thread>
#include <atomic>
#include <functional>
#include <png++/png.hpp>
#include <unistd.h>                         //for usleep
#include "fractals.hpp"

#define vcout if(verboseOutput) cout        //print only if verbose output enabled (verbose-cout)

using namespace std;

atomic<int> runningThreads(0);
extern int maxThreads;

///Prototypes
void printInfo( const size_t imageWidth, const size_t imageHeight, const size_t maxSectorSide,
                const int fractal_type, const long double offset_re, const long double offset_im,
                const bool juliaMode, const long double julia_re, const long double julia_im,
                const long double scalingFactor, const int max_iter, const long double bailout,
                const int colorMode);

///Functions
//main function to compute the entire fractal image
png::image<png::rgb_pixel> calc_and_color(
    const size_t imageWidth, const size_t imageHeight, const size_t maxSectorSide,                                          //image data
    const int fractal_type, const long double offset_re, const long double offset_im,                                       //fractal type, offset position
    const bool juliaMode, const long double julia_re, const long double julia_im,                                           //Julia mode, julia offset position
    const long double scalingFactor, const int max_iter, const long double bailout,                                         //zoom, iterations, bailout radius
    const int colorMode, const vector<png::rgb_pixel> palette, const bool verboseOutput){                                   //coloring mode, palette, user friendly console output

    //To better distribute the workload between all the threads, the image gets divided into sectors, then
    //when one thread working on a sector finishes, we launch another one to work on another sector, so that we
    //(almost) always have maxThreads working on a sector
    vector<array<size_t, 4>> sectors = {};
    //Image of the fractal
    png::image<png::rgb_pixel> fractal_image(imageWidth, imageHeight);

    //Creating sectors onto which thread can work in parallel
    for(size_t i = 0; i < imageHeight; i += maxSectorSide) {
        for(size_t j = 0; j < imageWidth; j += maxSectorSide) {
            size_t start_x = j;
            size_t start_y = i;
            size_t end_x   = j + min((size_t)maxSectorSide, imageWidth - j);
            size_t end_y   = i + min((size_t)maxSectorSide, imageHeight - i);

            sectors.push_back({start_x, start_y, end_x, end_y});
        }
    }
    size_t totalSectors = sectors.size();
    size_t finishedSectors = 0;

    /*
    //DEBUG STUFF
    size_t totalAreaSectors = 0;
    cout << "Image size: " << imageWidth << " x " << imageHeight << endl;
    cout << "# of sectors: " << totalSectors << endl;
    for(int i = 0; i < totalSectors; ++i){
        cout << "Sector #" << i << ":" << endl;
        cout << "   Start : (" << sectors[i][0] << ", " << sectors[i][1] << ")" << endl;
        cout << "   Finish: (" << sectors[i][2] << ", " << sectors[i][3] << ")" << endl;
        size_t area = (sectors[i][3] - sectors[i][1]) * (sectors[i][2] - sectors[i][0]);
        cout << "   Area  :  " << area << " pixels" << endl;
        totalAreaSectors += area;
    }
    cout << "Total area: " << totalAreaSectors << "    Image area: " << imageWidth*imageHeight << endl;
    return 0;
    */

    //Function pointer to the fractal function to call
    void (*fractal_fn_pointer)(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const long double offset_re, const long double offset_im, const bool juliaMode, const long double julia_re, const long double julia_im,
    const long double scalingFactor, const int max_iter, const long double bailout,
    const int color_mode, const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write);
    //Initialize function pointer based on selected fractal
    switch(fractal_type){
        case 0:             //MANDELBROT
            fractal_fn_pointer = &mandelbrot_color;
            break;

        case 1:             //BURNING SHIP
            fractal_fn_pointer = &burningShip_color;
            break;

        case 2:             //MANDELBAR
            fractal_fn_pointer = &mandelbar_color;
            break;

        case 3:             //SPADE FRACTAL
            fractal_fn_pointer = &spadefract_color;
            break;

        default:            //DEFAULT = MANDELBROT
            fractal_fn_pointer = &mandelbrot_color;
            break;
    }


    //Heavy lifting: spawning all the threads to compute the fractal
    if(verboseOutput) printInfo(imageWidth, imageHeight, maxSectorSide, fractal_type, offset_re, offset_im, juliaMode, julia_re, julia_im, scalingFactor, max_iter, bailout, colorMode);
    //While I still have sectors to perform computations on
    while(sectors.size() > 0) {
        //If I can spawn threads
        if(atomic_load(&runningThreads) < maxThreads) {
            size_t start_x = sectors.back()[0];
            size_t start_y = sectors.back()[1];
            size_t end_x   = sectors.back()[2];
            size_t end_y   = sectors.back()[3];
            sectors.pop_back();
            thread(fractal_fn_pointer,          //Fractal to launch
                   start_x, start_y,            //(x,y) starting position
                   end_x, end_y,                //(x,y) ending position
                   offset_re, offset_im,        //Offset in the complex plane
                   juliaMode,                   //Whether it should run in Julia mode or not
                   julia_re, julia_im,          //Offsets of the constant c in Julia mode
                   scalingFactor,               //Scaling factor
                   max_iter,                    //Maximum number of iterations
                   bailout,                     //Bailout radius
                   colorMode,                   //Coloring mode of the image
                   ref(palette),                //Reference to color palette
                   ref(fractal_image)           //Reference to image to update pixels
                   ).detach();

            runningThreads += 1;
            ++finishedSectors;
        }

        vcout << finishedSectors << "/" << totalSectors << " sectors finished\r";
        usleep(100);
    }

    //Wait for all the threads to finish
    while(atomic_load(&runningThreads) > 0){
        vcout << "Waiting " << runningThreads.load(memory_order_seq_cst) << " threads to finish\r";
        usleep(1000);
    }

    return fractal_image;
}

/*
size_t calc_all_iter(){
    return 0;
}

png::rgb_pixel calc_all_colors(){
    return png::rgb_pixel(0, 0, 0);
}

void color_image(){

}
*/

void printInfo( const size_t imageWidth, const size_t imageHeight, const size_t maxSectorSide,                          //Image data and sector side
                const int fractal_type, const long double offset_re, const long double offset_im,                       //Fractal type and offset from origin
                const bool juliaMode, const long double julia_re, const long double julia_im,                           //Julia mode and value of c
                const long double scalingFactor, const int max_iter, const long double bailout,                         //zoom, maximum iterations, bailout radius
                const int colorMode){                                                                                   //coloring mode

    string fractal_name;
    switch (fractal_type){
        case 0: fractal_name = "Mandelbrot set"; break;
        case 1: fractal_name = "Burning ship";   break;
        case 2: fractal_name = "Mandelbar set";  break;
        case 3: fractal_name = "Spade fractal";  break;
        default:
            fractal_name = "Unknown, defaulting to Mandelbrot";
            break;
    }

    string coloring_mode_name;
    switch (colorMode){
        case 0: coloring_mode_name = "binary (b/w)";  break;
        case 1: coloring_mode_name = "ln";            break;
        case 2: coloring_mode_name = "sqrt(iter)";    break;
        case 3: coloring_mode_name = "cbrt(iter)";    break;
        case 4: coloring_mode_name = "atan2(im, re)"; break;
        default:
            coloring_mode_name = "unknown, defaulting to binary";
            break;
    }

    cout << "Rendering: " << fractal_name << (juliaMode ? " (J)\n" : "\n");
    cout << "Image properties: " << imageWidth << " x " << imageHeight << endl;
    cout << "Sectors up to " << maxSectorSide << "x" << maxSectorSide << " = " << maxSectorSide * maxSectorSide << " pixels^2" << endl;
    cout << "Iterations: " << max_iter << endl;
    cout << "Bailout r.: " << bailout << endl;
    cout << "Coloring mode: " << coloring_mode_name << " (" << colorMode << ")" << endl;
    cout << "Centered at: (" << offset_re << ", " << offset_im << " i)" << endl;
    if(juliaMode) cout << "Julia c: (" << julia_re << ", " << julia_im << " i)" << endl;
    //Auxiliary variable, used when calculating scaling from pixel space to complex plane, to prevent stretching
    const auto squareScale = min(imageWidth, imageHeight);
    //Top left corner
    const long double start_re = 4 * (- (long double)imageWidth / 2)  / ((long double)squareScale * scalingFactor) + offset_re;
    const long double start_im = 4 * (- (long double)imageHeight / 2) / ((long double)squareScale * scalingFactor) + offset_im;
    //Bottom right
    const long double stop_re  = 4 * ((long double)imageWidth / 2)  / ((long double)squareScale * scalingFactor) + offset_re;
    const long double stop_im  = 4 * ((long double)imageHeight / 2) / ((long double)squareScale * scalingFactor) + offset_im;
    cout << "Span: (" << stop_re - start_re << " x " << stop_im - start_im << ")" << endl;
    cout << "Scale: " << scalingFactor << endl;
}

#endif // CALC_AND_COLOR_HPP_INCLUDED
