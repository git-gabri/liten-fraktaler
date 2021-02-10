#include <iostream>
#include <string>
#include <thread>
#include <png++/png.hpp>
#include <getopt.h>         //Needed for parsing command line arguments
#include <structs.h>
#include "calc_and_color.hpp"
#include "misc_functions.hpp"

//TODO: allow image segmentation, i.e. subdivide big images that wouldn't fit in RAM, into multiple
//smaller images, compute each of the individually, and save them on disk
#define MAX_IMAGE_WIDTH  1000000        //An image 1000000*1000000 would require >3TB of RAM
#define MAX_IMAGE_HEIGHT 1000000
#define MAX_SECTOR_SIDE  1000000

#define FALLBACK_NUM_THREADS 4


using namespace std;
using namespace png;

int maxThreads = FALLBACK_NUM_THREADS;

/*Supported flags
* -H            show help
* -w SIZE_T     sets width of the image
* -h SIZE_T     sets height of the image
* -S SIZE_T     sets maximum sector side
* -o STRING     output file name
* -f INT        sets fractal type
* -r LNG_DBL    sets real part of the offset from the origin
* -i LNG_DBL    sets imag part of the offset from the origin
* -j            enable Julia mode
* -A LNG_DBL    sets real part of the constant c in Julia mode              //This will be replaced with somehting more readable when I'll figure out how multi letter flags work
* -B LNG_DBL    sets imag part of the constant c in Julia mode              //This will be replaced with somehting more readable when I'll figure out how multi letter flags work
* -s LNG_DBL    sets scaling factor
* -t INT        sets maximum iteration
* -b LNG_DBL    sets bailout
* -T INT        sets maximum number of threads
* -p STRING     sets filename of config file for color palette
* -c INT        sets coloring mode
* -v            verbose output, more user friendly

TODO
* ADD MULTI LETTER AND MULTI PARAMETER FLAGS
* -R LNG_DBL    sets fractal rotation
* plane warping option
* -d            automatically set iterations
* pickover stalks
*/
//Parse flags and user input, sanitize it and if something is wrong print a corresponfing error
int parseArgv(
    int argc, char *argv[],             //argc and argv from main
    imagesettings_t& isettings,         //struct containing various image related settings
    fractalsettings_t& fsettings,       //struct containing various fractal related settings
    colorsettings_t& csettings,         //struct containing various color related settings
    bool& verbose                       //verbose output flag
    ){

    ///Parsing command line flags
    for(;;) {
        switch(getopt(argc, argv, "w:h:S:o:f:r:i:jA:B:s:t:b:T:p:c:vH")) {
        //------------------------------------------
        //Image related flags
        case 'w':       //Width
            isettings.imageWidth = stoul(optarg);
            continue;
        case 'h':       //Height
            isettings.imageHeight = stoul(optarg);
            continue;
        case 'S':       //Sector side
            isettings.maxSectorSide = stoul(optarg);
            continue;
        case 'o':       //Output filename
            isettings.imageName = optarg;
            continue;
        //------------------------------------------
        //Fractal related flags
        case 'f':       //Fractal type
            fsettings.fractal_type = stoi(optarg);
            continue;
        case 'r':       //Offset real
            fsettings.offset_re = stold(optarg);
            continue;
        case 'i':       //Offset imaginary
            fsettings.offset_im = stold(optarg);
            continue;
        case 'j':       //Run in Julia mode
            fsettings.juliaMode = true;
            continue;
        case 'A':       //Real part of c in Julia mode
            fsettings.julia_re = stold(optarg);
            continue;
        case 'B':       //Imag part of c in Julia mode
            fsettings.julia_im = stold(optarg);
            continue;
        case 's':       //Scaling
            fsettings.scalingFactor = stold(optarg);
            continue;
        case 't':       //Maximum iterations
            fsettings.max_iter = stoi(optarg);
            continue;
        case 'b':
            fsettings.bailout = stold(optarg);
            continue;
        case 'T':       //Maximum threads
            maxThreads = stoi(optarg);
            continue;
        //------------------------------------------
        //Colors related flags
        case 'p':       //Filename of config file for color palette
            csettings.paletteConfig = optarg;
            continue;
        case 'c':       //Coloring mode
            csettings.colorMode = stoi(optarg);
            continue;
        //------------------------------------------
        //Other flags
        case 'v':       //Verbose output
            verbose = true;
            continue;
        //------------------------------------------
        //Help
        case 'H':
        case '?':
        default :
            printHelp(argv);
            return 1;
            break;

        case -1:
            break;
        }

        break;
    }

    ///Checking user input
    //This flag is used so that if there are multiple errors in the user input, all the errors will be displayed at once
    //and the user doesn't have to restart the program over and over to check that all of the flags are correct
    bool retval = 0;
    //-----------------------------------------------------------------------------------------------------
    //Checking image related variables
    if(isettings.imageWidth > MAX_IMAGE_WIDTH){
        cerr << "ERROR: image width too large. Maximum allowed: " << MAX_IMAGE_WIDTH << endl;
        retval = 1;
    }
    if(isettings.imageHeight > MAX_IMAGE_HEIGHT){
        cerr << "ERROR: image height too high. Maximum allowed: " << MAX_IMAGE_HEIGHT << endl;
        retval = 1;
    }
    if(isettings.maxSectorSide > MAX_SECTOR_SIDE){
        cerr << "ERROR: sector side too large. Maximum allowed: " << MAX_SECTOR_SIDE << endl;
        retval = 1;
    }
    if(isettings.imageName == ""){
        cerr << "ERROR: empty image name" << endl;
        retval = 1;
    }
    //-----------------------------------------------------------------------------------------------------
    //Checking fractal related variables
    if(!isfinite(fsettings.offset_re)){
        cerr << "ERROR: real offset must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(fsettings.offset_im)){
        cerr << "ERROR: imaginary offset must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(fsettings.julia_re)){
        cerr << "ERROR: real offset of c in Julia mode must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(fsettings.offset_im)){
        cerr << "ERROR: imaginary offset of c in Julia mode must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(fsettings.scalingFactor)){
        cerr << "ERROR: scaling must be a finite number" << endl;
        retval = 1;
    }
    if(fsettings.scalingFactor <= 0){
        cerr << "ERROR: scaling must be a positive number" << endl;
        retval = 1;
    }
    if(fsettings.max_iter < 0){
        cerr << "ERROR: maximum number of iterations must be a positive number" << endl;
        retval = 1;
    }
    if(fsettings.bailout > 1.0e+2400l){
        cerr << "ERROR: bailout must be less than 1.0e+2400" << endl;
        retval = 1;
    }
    if(maxThreads < 0){
        cerr << "ERROR: maximum number of threads must be a positive number" << endl;
        retval = 1;
    }
    //-----------------------------------------------------------------------------------------------------
    //Checking colors related variables
    if(csettings.paletteConfig == ""){
        cerr << "ERROR: empty palette config file name" << endl;
        retval = 1;
    }

    return retval;
}

//Minibrot at
// -1.870344820894
// 0
// zoom = 1.0e10
//FOUND IT :)

//Minibrot at
// -1.992577750638
// 0
// zoom = 1.0e12
//FOUND IT! :)

int main(int argc, char *argv[]) {
    imagesettings_t isettings{};
    fractalsettings_t fsettings{};
    colorsettings_t csettings{};

    maxThreads = thread::hardware_concurrency();
    if(maxThreads == 0) maxThreads = FALLBACK_NUM_THREADS;

    vector<png::rgb_pixel> palette = {};

    bool verboseOutput = false;
    {
        bool errors = parseArgv(argc, argv, isettings, fsettings, csettings, verboseOutput);

        if(errors) return 1;

        //If the loading of the palette returns 1, it failed
        if(load_palette(palette, csettings.paletteConfig)){
            cerr << "WARN: palette could not be loaded correctly, defaulting to b/w" << endl;
            palette.clear();
            palette.push_back(png::rgb_pixel{0, 0, 0});         //Load black
            palette.push_back(png::rgb_pixel{255, 255, 255});   //Load white
        }
    }
    if(palette.size() == 0){cerr << "THIS SHOULDN'T HAVE HAPPENED" << endl; return 1;}

    //Declaring and initializing image after parsing user input, by assigning it the value of the computed fractal image
    chrono::high_resolution_clock::time_point start_time = chrono::high_resolution_clock::now();
    image<rgb_pixel> fractalImage = calc_and_color(isettings, fsettings, csettings, palette, verboseOutput);
    chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();

    vcout << "                                        \r";
    vcout << "Rendering took " << chrono::duration_cast<chrono::duration<double>>(end_time - start_time).count() << " seconds" << endl;
    vcout << "Writing image..." << endl;
    isettings.imageName = isettings.imageName + ".png";
    fractalImage.write(isettings.imageName);
    vcout << "Done!" << endl;

    return EXIT_SUCCESS;
}
