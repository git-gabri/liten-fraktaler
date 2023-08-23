#include <iostream>
#include <vector>
#include <string>
#include <png++/png.hpp>

#include "liten.hpp"
#include "stopwatch.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    lf::init();
    if(lf::parse_options(vector<string>(argv + 1, argc + argv)))
        return EXIT_FAILURE;

    const auto verbose = lf::internals::consettings.verbose_output;

    lf::load_palette();

    if(lf::internals::palette.size() == 0){cerr << "THIS SHOULDN'T HAVE HAPPENED" << endl; return 1;}

    //Declaring and initializing image after parsing user input, by assigning it the value of the computed fractal image
    stopwatch s;
    s.tic();
    png::image<png::rgb_pixel> fractal_image;
    const auto render_status = lf::launch_render(fractal_image);
    s.toc();

    if(render_status != 0){
        cerr << "[ERROR]: couldn't properly render the image with the specified options" << endl;
        return EXIT_FAILURE;
    }

    if(verbose){
        cout << "Rendering took " << s << " seconds" << endl;
        cout << "Writing image..." << endl;
    }
    fractal_image.write(lf::internals::isettings.image_name + ".png");
    if(verbose)
        cout << "Done!" << endl;

    return EXIT_SUCCESS;
}