#ifndef CALC_AND_COLOR_HPP_INCLUDED
#define CALC_AND_COLOR_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <functional>
#include <png++/png.hpp>

#include <ThreadPool.h>
#include <future>

#include <structs.h>
#include "fractals.hpp"

#define vcout if(verboseOutput) cout        //print only if verbose output enabled (verbose-cout)

using namespace std;

extern int maxThreads;

//printinfo prototype
void printInfo(const imagesettings_t isettings, const fractalsettings_t fsettings, const colorsettings_t csettings);

///Functions
//main function to compute the entire fractal image
png::image<png::rgb_pixel> calc_and_color(imagesettings_t isettings, fractalsettings_t fsettings, colorsettings_t csettings,
                                          const vector<png::rgb_pixel> palette, const bool verboseOutput){

    //To better distribute the workload between all the threads, the image gets divided into sectors, then
    //when one thread working on a sector finishes, we launch another one to work on another sector, so that we
    //(almost) always have maxThreads working on a sector
    vector<array<size_t, 4>> sectors = {};
    //Image of the fractal
    png::image<png::rgb_pixel> fractal_image(isettings.imageWidth, isettings.imageHeight);

    //Creating sectors onto which thread can work in parallel
    for(size_t i = 0; i < isettings.imageHeight; i += isettings.maxSectorSide) {
        for(size_t j = 0; j < isettings.imageWidth; j += isettings.maxSectorSide) {
            size_t start_x = j;
            size_t start_y = i;
            size_t end_x   = j + min((size_t)isettings.maxSectorSide, isettings.imageWidth - j);
            size_t end_y   = i + min((size_t)isettings.maxSectorSide, isettings.imageHeight - i);

            sectors.push_back({start_x, start_y, end_x, end_y});
        }
    }
    size_t totalSectors = sectors.size();

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
    const fractalsettings_t fsettings, const colorsettings_t csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write);
    //Initialize function pointer based on selected fractal
    switch(fsettings.fractal_type){
        case 0:             //MANDELBROT
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &mandelbrot;
            break;

        case 1:             //TIPPETS MANDELBROT
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &tippets_mandelbrot;
            break;

        case 2:             //BURNING SHIP
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &burningShip;
            break;

        case 3:             //MANDELBAR
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &mandelbar;
            break;

        case 4:             //SPADE FRACTAL
            if(fsettings.bailout < 0) fsettings.bailout = 4;
            fractal_fn_pointer = &spadefract;
            break;

        case 5:             //MAGNET 1
            if(fsettings.bailout < 0) fsettings.bailout = 128;
            fractal_fn_pointer = &magnet_type1;
            break;

        case 6:             //MAGNET 2
            if(fsettings.bailout < 0) fsettings.bailout = 1024;
            fractal_fn_pointer = &magnet_type2;
            break;

        case 7:             //CACTUS
            if(fsettings.bailout < 0) fsettings.bailout = 8;
            fractal_fn_pointer = &cactus;
            break;

        case 8:             //ZUBIETA
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &zubieta;
            break;

        case 9:             //ZUBITHETA
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &zubitheta;
            break;

        case 10:            //LOGISTIC MAP
            if(fsettings.bailout < 0) fsettings.bailout = 2;
            fractal_fn_pointer = &logistic_map;
            break;

        case 99:             //BINOCULAR m = 4
            if(fsettings.bailout < 0) fsettings.bailout = 1024;
            fractal_fn_pointer = &binocular_m4;
            break;

        default:            //DEFAULT = MANDELBROT
            fractal_fn_pointer = &mandelbrot;
            break;
    }

    //Print info if required
    if(verboseOutput) printInfo(isettings, fsettings, csettings);

    //Create threadpool for rendering
    ThreadPool renderpool(maxThreads);
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
        if(verboseOutput) cout << "Completed sectors: " << i << "/" << totalSectors << "\r";
    }

    return fractal_image;
}

//----------------------------------------------------------------------------------------------------------------------------
void printInfo(const imagesettings_t isettings, const fractalsettings_t fsettings, const colorsettings_t csettings){

    string fractal_name;
    switch (fsettings.fractal_type){
        case  0: fractal_name = "Mandelbrot set";   break;
        case  1: fractal_name = "Tippets Mandel";   break;
        case  2: fractal_name = "Burning ship";     break;
        case  3: fractal_name = "Mandelbar set";    break;
        case  4: fractal_name = "Spade fractal";    break;
        case  5: fractal_name = "Magnet type I";    break;
        case  6: fractal_name = "Magnet type II";   break;
        case  7: fractal_name = "Cactus";           break;
        case  8: fractal_name = "Zubieta";          break;
        case  9: fractal_name = "Zubitheta";        break;
        case 10: fractal_name = "Logistic map";     break;
        default:
            fractal_name = "Unknown/undocumented";
            break;
    }

    string coloring_mode_name;
    switch (csettings.colorMode){
        case 0: coloring_mode_name = "binary (b/w)";  break;
        case 1: coloring_mode_name = "linear";        break;
        case 2: coloring_mode_name = "ln";            break;
        case 3: coloring_mode_name = "sqrt(iter)";    break;
        case 4: coloring_mode_name = "cbrt(iter)";    break;
        case 5: coloring_mode_name = "S curve";       break;
        case 6: coloring_mode_name = "atan2(im, re)"; break;
        default:
            coloring_mode_name = "Unknown/undocumented";
            break;
    }

    cout << "Rendering      : " << fractal_name << (fsettings.juliaMode ? " (J)" : "") << " (" << maxThreads << " thread(s))" << endl;
    cout << "Image size     : " << isettings.imageWidth << " x " << isettings.imageHeight << endl;
    cout << "Sectors up to  : " << isettings.maxSectorSide << "x" << isettings.maxSectorSide << " = " << isettings.maxSectorSide * isettings.maxSectorSide << " pixels^2" << endl;

    cout << "Iterations     : " << fsettings.max_iter << endl;
    cout << "Bailout radius : " << fsettings.bailout << endl;
    cout << "Coloring mode  : " << coloring_mode_name << " (" << csettings.colorMode << ")" << endl;
    cout << "Centered at    : (" << fsettings.offset_re << ", " << fsettings.offset_im << " i)" << endl;
    if(fsettings.juliaMode) cout << "Julia c const. : (" << fsettings.julia_re << ", " << fsettings.julia_im << " i)" << endl;

    //Auxiliary variable, used when calculating scaling from pixel space to complex plane, to prevent stretching
    const auto squareScale = min(isettings.imageWidth, isettings.imageHeight);
    //Top left corner
    const long double start_re = 4 * (- (long double)isettings.imageWidth / 2)  / ((long double)squareScale * fsettings.scalingFactor) + fsettings.offset_re;
    const long double start_im = 4 * (- (long double)isettings.imageHeight / 2) / ((long double)squareScale * fsettings.scalingFactor) + fsettings.offset_im;
    //Bottom right
    const long double stop_re  = 4 * ((long double)isettings.imageWidth / 2)  / ((long double)squareScale * fsettings.scalingFactor) + fsettings.offset_re;
    const long double stop_im  = 4 * ((long double)isettings.imageHeight / 2) / ((long double)squareScale * fsettings.scalingFactor) + fsettings.offset_im;
    cout << "Span  : (" << stop_re - start_re << " x " << stop_im - start_im << ")" << endl;
    cout << "Scale : " << fsettings.scalingFactor << endl;
}

#endif // CALC_AND_COLOR_HPP_INCLUDED
