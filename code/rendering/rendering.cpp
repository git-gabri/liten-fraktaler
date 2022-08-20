#include "liten.hpp"
#include "threadpool.hpp"
#include "misc_functions.hpp"
#include "coloring.hpp"
#include "fractals.hpp"

#include <complex>
#include <functional>
#include <map>
#include <future>

using namespace std;
using namespace lf::internals;

#define vcout if(consettings.verbose_output) cout        //print only if verbose output enabled (verbose-cout)

//Function that makes multiple threads to render the image
png::image<png::rgb_pixel> lf::launch_render(){
    //To better distribute the workload between all the threads, the image gets divided into sectors, then
    //when one thread working on a sector finishes, we launch another one to work on another sector, so that we
    //(almost) always have all the threads working on a sector
    vector<array<size_t, 4>> sectors = {};
    //Image of the fractal
    vcout << "Allocating image in RAM... " << flush;
    png::image<png::rgb_pixel> fractal_image(isettings.image_width, isettings.image_height);
    vcout << "Done!" << endl;

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
    size_t total_sectors = sectors.size();

    /*
    //DEBUG STUFF
    size_t totalAreaSectors = 0;
    cout << "Image size: " << image_width << " x " << image_height << endl;
    cout << "# of sectors: " << total_sectors << endl;
    for(int i = 0; i < total_sectors; ++i){
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
    block_renderer_fn_ptr_t block_renderer_pointer = get_block_renderer_ptr(fsettings.fractal_type);

    //Init bailout radius to default if necessary
    if(fsettings.bailout < 0)
        fsettings.bailout = default_bailout_radius(fsettings.fractal_type);

    //Print info if required
    if(consettings.verbose_output) print_render_info();

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
                block_renderer_pointer,      //Block renderer to launch to launch
                start_x, start_y,            //(x,y) starting position
                end_x, end_y,                //(x,y) ending position
                ref(fractal_image)           //Reference to image to update pixels
            )
        );
    }

    //Once all the jobs are enqueued, wait for all of them to finish
    for(size_t i = 0; i < completed_sectors.size(); ++i){
        completed_sectors[i].get();
        vcout << "Completed sectors: " << i << "/" << total_sectors << "\r" << flush;
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

block_renderer_fn_ptr_t lf::internals::get_block_renderer_ptr(const ftype& fractal_type){
    block_renderer_fn_ptr_t ret_ptr;

    switch(fractal_type){
        default:
        case ftype::mandelbrot:     ret_ptr = &block_renderer<mandelbrot<long double>>;             break;
        case ftype::tippets:        ret_ptr = &block_renderer<tippets_mandelbrot<long double>>;     break;
        case ftype::burnship:       ret_ptr = &block_renderer<burning_ship<long double>>;           break;
        case ftype::mandelbar:      ret_ptr = &block_renderer<mandelbar<long double>>;              break;
        case ftype::magnet1:        ret_ptr = &block_renderer<magnet_type1<long double>>;           break;
        case ftype::magnet2:        ret_ptr = &block_renderer<magnet_type2<long double>>;           break;
        case ftype::cactus:         ret_ptr = &block_renderer<cactus<long double>>;                 break;
        case ftype::cactus0:        ret_ptr = &block_renderer<cactus<long double>>;                 break;
        case ftype::zubieta:        ret_ptr = &block_renderer<zubieta<long double>>;                break;
        case ftype::zubitheta:      ret_ptr = &block_renderer<zubitheta<long double>>;              break;
        case ftype::logmap:         ret_ptr = &block_renderer<logistic_map<long double>>;           break;
        case ftype::unpolsquare:    ret_ptr = &block_renderer<unpol_square<long double>>;           break;
        case ftype::moth:           ret_ptr = &block_renderer<moth<long double>>;                   break;
        case ftype::cave:           ret_ptr = &block_renderer<cave<long double>>;                   break;
        case ftype::wankel:         ret_ptr = &block_renderer<wankel<long double>>;                 break;
        case ftype::seaangel:       ret_ptr = &block_renderer<sea_angel<long double>>;              break;
        case ftype::smith:          ret_ptr = &block_renderer<smith<long double>>;                  break;

        case ftype::spade:          ret_ptr = &block_renderer<spadefract<long double>>;             break;

        case ftype::test0:          ret_ptr = &block_renderer<test0<long double>>;                  break;
        case ftype::test1:          ret_ptr = &block_renderer<test1<long double>>;                  break;
        case ftype::test2:          ret_ptr = &block_renderer<test2<long double>>;                  break;
        case ftype::test3:          ret_ptr = &block_renderer<test3<long double>>;                  break;
        case ftype::test4:          ret_ptr = &block_renderer<test4<long double>>;                  break;
        case ftype::test5:          ret_ptr = &block_renderer<test5<long double>>;                  break;
        case ftype::test6:          ret_ptr = &block_renderer<test6<long double>>;                  break;
        case ftype::test7:          ret_ptr = &block_renderer<test7<long double>>;                  break;
        case ftype::test8:          ret_ptr = &block_renderer<test8<long double>>;                  break;
        case ftype::test9:          ret_ptr = &block_renderer<test9<long double>>;                  break;
    }

    return ret_ptr;
}

long double lf::internals::default_bailout_radius(const ftype& f){
    const map<ftype, long double> default_bailout_radiuses{
        {ftype::mandelbrot,     2},
        {ftype::tippets,        2},
        {ftype::burnship,       2},
        {ftype::mandelbar,      2},
        {ftype::magnet1,      128},
        {ftype::magnet2,     1024},
        {ftype::cactus,         8},
        {ftype::cactus0,        8},
        {ftype::zubieta,        2},
        {ftype::zubitheta,      2},
        {ftype::logmap,         2},
        {ftype::unpolsquare,    2},
        {ftype::moth,           4},
        {ftype::cave,          16},
        {ftype::wankel,         4},
        {ftype::seaangel,      64},
        {ftype::smith,         64},
        {ftype::spade,        128},
        
        {ftype::test0,   128},
        {ftype::test1,   128},
        {ftype::test2,   128},
        {ftype::test3,   128},
        {ftype::test4,   128},
        {ftype::test5,   128},
        {ftype::test6,   128},
        {ftype::test7,   128},
        {ftype::test8,   128},
        {ftype::test9,   128}
    };

    long double ret_rad;
    if(default_bailout_radiuses.contains(f))
        ret_rad = default_bailout_radiuses.at(f);
    else
        ret_rad = -1;

    return ret_rad;
}