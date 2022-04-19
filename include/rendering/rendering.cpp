#include "threadpool.hpp"
#include "rendering.hpp"
#include "coloring.hpp"
#include "fractals.hpp"
#include "misc_functions.hpp"

#include <string>
#include <functional>
#include <map>
#include <future>

using namespace std;

#define vcout if(consettings.verbose_output) cout        //print only if verbose output enabled (verbose-cout)

//Function that makes multiple threads to render the image
png::image<png::rgb_pixel> launch_render(   imagesettings_t& isettings,
                                            fractalsettings_t& fsettings,
                                            colorsettings_t& csettings,
                                            rendersettings_t& rsettings,
                                            consolesettings_t& consettings,
                                            const vector<png::rgb_pixel> palette){

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
    void (*block_renderer_pointer)(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t& fsettings, const colorsettings_t& csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write);
    block_renderer_pointer = &block_renderer;

    //Init bailout radius to default if necessary
    if(fsettings.bailout < 0)
        fsettings.bailout = init_bailout_radius(fsettings.fractal_type);

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
                block_renderer_pointer,      //Block renderer to launch to launch
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

//Main function that runs in a thread whose job is to render the colors in a specific block in the image
void block_renderer(
    const size_t startX, const size_t startY, const size_t endX, const size_t endY,
    const fractalsettings_t& fsettings, const colorsettings_t& csettings,
    const vector<png::rgb_pixel>& palette, png::image<png::rgb_pixel>& image_to_write){

    const auto width = image_to_write.get_width();
    const auto height = image_to_write.get_height();
    //In order not to have an image which is stretched either along the x or the y axis,
    //we introduce an auxiliary variable which will be useful for scaling, which represents
    //the side of a square whose length is the minimum of the width and the height
    const auto square_scale = min(width, height);

    //Square of the bailout radius to chech the norm against it in the main iteration loop
    const long double bailout = fsettings.bailout;
    const long double bailout_2 = bailout * bailout;

    complex<long double> z;
    complex<long double> c;
    bool stop_iterating = false;
    vector<complex<long double>> history;
    vector<complex<long double>> extra_params;

    //Assign fractal_fn_pointer wrt selected fractal
    fractal_fn_ptr_t fractal_fn_pointer = init_fractal_fn_pointer(fsettings);

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
                z = (*fractal_fn_pointer)(z, c, history, extra_params, fsettings);

                if(norm(z) > bailout_2)
                    stop_iterating = true;

                ++iter_count;
            } while((iter_count < fsettings.max_iter) && !stop_iterating);

            image_to_write[y][x] = compute_color(iter_count, fsettings.max_iter, z.real(), z.imag(), fsettings.bailout, csettings.cmode, palette);
        }
    }

    return;
}

long double init_bailout_radius(const ftype& f){
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

fractal_fn_ptr_t init_fractal_fn_pointer(const fractalsettings_t& fsettings){
    fractal_fn_ptr_t ret_ptr;

    switch(fsettings.fractal_type){
        default:
        case ftype::mandelbrot:     ret_ptr = &mandelbrot<long double>;         break;
        case ftype::tippets:        ret_ptr = &tippets_mandelbrot<long double>; break;
        case ftype::burnship:       ret_ptr = &burning_ship<long double>;       break;
        case ftype::mandelbar:      ret_ptr = &mandelbar<long double>;          break;
        case ftype::magnet1:        ret_ptr = &magnet_type1<long double>;       break;
        case ftype::magnet2:        ret_ptr = &magnet_type2<long double>;       break;
        case ftype::cactus:         ret_ptr = &cactus<long double>;             break;
        case ftype::cactus0:        ret_ptr = &cactus<long double>;             break;
        case ftype::zubieta:        ret_ptr = &zubieta<long double>;            break;
        case ftype::zubitheta:      ret_ptr = &zubitheta<long double>;          break;
        case ftype::logmap:         ret_ptr = &logistic_map<long double>;       break;
        case ftype::unpolsquare:    ret_ptr = &unpol_square<long double>;       break;
        case ftype::moth:           ret_ptr = &moth<long double>;               break;
        case ftype::cave:           ret_ptr = &cave<long double>;               break;
        case ftype::wankel:         ret_ptr = &wankel<long double>;             break;
        case ftype::seaangel:       ret_ptr = &sea_angel<long double>;          break;
        case ftype::smith:          ret_ptr = &smith<long double>;              break;

        case ftype::spade:          ret_ptr = &spadefract;                      break;

        case ftype::test0:          ret_ptr = &test0<long double>;              break;
        case ftype::test1:          ret_ptr = &test1<long double>;              break;
        case ftype::test2:          ret_ptr = &test2<long double>;              break;
        case ftype::test3:          ret_ptr = &test3<long double>;              break;
        case ftype::test4:          ret_ptr = &test4<long double>;              break;
        case ftype::test5:          ret_ptr = &test5<long double>;              break;
        case ftype::test6:          ret_ptr = &test6<long double>;              break;
        case ftype::test7:          ret_ptr = &test7<long double>;              break;
        case ftype::test8:          ret_ptr = &test8<long double>;              break;
        case ftype::test9:          ret_ptr = &test9<long double>;              break;
    }

    return ret_ptr;
}