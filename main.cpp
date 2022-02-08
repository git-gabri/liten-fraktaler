#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <png++/png.hpp>
#include "structs.hpp"
#include "option_parsing.hpp"
#include "rendering.hpp"
#include "coloring.hpp"
#include "misc_functions.hpp"

#define MAX_IMAGE_WIDTH  1000000        //An image 1000000*1000000 would require >3TB of RAM
#define MAX_IMAGE_HEIGHT 1000000
#define MAX_SECTOR_SIDE  1000000

#define FALLBACK_NUM_THREADS 1
#define vcout if(consettings.verbose_output) cout        //print only if verbose output enabled (verbose-cout)

using namespace std;
using namespace png;

int main(int argc, char *argv[]) {
    fractalsettings_t fsettings{};
    imagesettings_t isettings{};
    colorsettings_t csettings{};
    rendersettings_t rsettings{};
    consolesettings_t consettings{};

    rsettings.max_threads = thread::hardware_concurrency();
    if(rsettings.max_threads == 0) rsettings.max_threads = FALLBACK_NUM_THREADS;

    vector<png::rgb_pixel> palette = {};

    {
        const vector<string> arv_vec(argv + 1, argv + argc);
        if(parse_options(arv_vec, isettings, fsettings, csettings, rsettings, consettings)) return 2;

        //If the loading of the palette returns 1, it failed
        if(load_palette(palette, csettings.palette_filename)){
            print_warning("palette could not be loaded correctly, defaulting to b/w");
            palette.clear();
            palette.push_back(png::rgb_pixel{0, 0, 0});         //Load black
            palette.push_back(png::rgb_pixel{255, 255, 255});   //Load white
        }
    }
    if(palette.size() == 0){cerr << "THIS SHOULDN'T HAVE HAPPENED" << endl; return 1;}

    //Declaring and initializing image after parsing user input, by assigning it the value of the computed fractal image
    chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();
    image<rgb_pixel> fractalImage = launch_render(isettings, fsettings, csettings, rsettings, consettings, palette);
    chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();

    vcout << "                                        \r";
    vcout << "Rendering took " << chrono::duration_cast<chrono::duration<double>>(end_time - start_time).count() << " seconds" << endl;
    vcout << "Writing image..." << endl;
    isettings.image_name = isettings.image_name + ".png";
    fractalImage.write(isettings.image_name);
    vcout << "Done!" << endl;

    return EXIT_SUCCESS;
}