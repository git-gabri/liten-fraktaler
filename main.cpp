#include <iostream>
#include <string>
#include <png++/png.hpp>
#include <getopt.h>         //Needed for parsing command line arguments
#include "calc_and_color.hpp"
#include "misc_functions.hpp"

//TODO: allow image segmentation, i.e. subdivide big images that wouldn't fit in RAM, into multiple
//smaller images, compute each of the individually, and save them on disk
#define MAX_IMAGE_WIDTH  1000000        //An image 1000000*1000000 would require >3TB of RAM
#define MAX_IMAGE_HEIGHT 1000000
#define MAX_SECTOR_SIDE  1000000
#define MAX_FRACTAL_TYPE 3              //Because only 4 fractal types [0,3] are implemented
#define MAX_COLOR_MODE   4              //Because only 5 color modes [0,4] are implemented


using namespace std;
using namespace png;

int maxThreads = 8;

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
    int argc, char *argv[],                                                                             //argc and argv from main
    size_t &imageWidth, size_t &imageHeight, size_t &maxSectorSide,                                     //image width, height and maximum sector side
    int &fractal_type, long double &offset_re, long double &offset_im,                                  //fractal type, real and imaginary offsets
    bool &juliaMode, long double &julia_re, long double &julia_im,                                      //run in julia mode, julia offsets (i.e. setting c constant)
    long double &scalingFactor, int &max_iter, long double &bailout,                                    //zoom, maximum iterations, coloring mode, filename of config file for color palette
    int &colorMode, string &palette_filename, string &image_filename, bool &verbose){                   //filename of output image file, verbose output flag

    ///Parsing command line flags
    for(;;) {
        switch(getopt(argc, argv, "w:h:S:o:f:r:i:jA:B:s:t:b:T:p:c:vH")) {
        //------------------------------------------
        //Image related flags
        case 'w':       //Width
            imageWidth = stoul(optarg);
            continue;
        case 'h':       //Height
            imageHeight = stoul(optarg);
            continue;
        case 'S':       //Sector side
            maxSectorSide = stoul(optarg);
            continue;
        case 'o':       //Output filename
            image_filename = optarg;
            continue;
        //------------------------------------------
        //Fractal related flags
        case 'f':       //Fractal type
            fractal_type = stoi(optarg);
            continue;
        case 'r':       //Offset real
            offset_re = stold(optarg);
            continue;
        case 'i':       //Offset imaginary
            offset_im = stold(optarg);
            continue;
        case 'j':       //Run in Julia mode
            juliaMode = true;
            continue;
        case 'A':       //Real part of c in Julia mode
            julia_re = stold(optarg);
            continue;
        case 'B':       //Imag part of c in Julia mode
            julia_im = stold(optarg);
            continue;
        case 's':       //Scaling
            scalingFactor = stold(optarg);
            continue;
        case 't':       //Maximum iterations
            max_iter = stoi(optarg);
            continue;
        case 'b':
            bailout = stold(optarg);
            continue;
        case 'T':       //Maximum threads
            maxThreads = stoi(optarg);
            continue;
        //------------------------------------------
        //Colors related flags
        case 'p':       //Filename of config file for color palette
            palette_filename = optarg;
            continue;
        case 'c':       //Coloring mode
            colorMode = stoi(optarg);
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
    if(imageWidth > MAX_IMAGE_WIDTH){
        cerr << "ERROR: image width too large. Maximum allowed: " << MAX_IMAGE_WIDTH << endl;
        retval = 1;
    }
    if(imageHeight > MAX_IMAGE_HEIGHT){
        cerr << "ERROR: image height too high. Maximum allowed: " << MAX_IMAGE_HEIGHT << endl;
        retval = 1;
    }
    if(maxSectorSide > MAX_SECTOR_SIDE){
        cerr << "ERROR: sector side too large. Maximum allowed: " << MAX_SECTOR_SIDE << endl;
        retval = 1;
    }
    if(image_filename == ""){
        cerr << "ERROR: empty image name" << endl;
        retval = 1;
    }
    //-----------------------------------------------------------------------------------------------------
    //Checking fractal related variables
    if(fractal_type < 0 || fractal_type > MAX_FRACTAL_TYPE){
        cerr << "ERROR: fractal type out of range. Must be in [0," << MAX_FRACTAL_TYPE << "]" << endl;
        retval = 1;
    }
    if(!isfinite(offset_re)){
        cerr << "ERROR: real offset must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(offset_im)){
        cerr << "ERROR: imaginary offset must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(julia_re)){
        cerr << "ERROR: real offset of c in Julia mode must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(offset_im)){
        cerr << "ERROR: imaginary offset of c in Julia mode must be a finite number" << endl;
        retval = 1;
    }
    if(!isfinite(scalingFactor)){
        cerr << "ERROR: scaling must be a finite number" << endl;
        retval = 1;
    }
    if(scalingFactor <= 0){
        cerr << "ERROR: scaling must be a positive number" << endl;
        retval = 1;
    }
    if(max_iter < 0){
        cerr << "ERROR: maximum number of iterations must be a positive number" << endl;
        retval = 1;
    }
    if(bailout > 1.0e+2400l){
        cerr << "ERROR: bailout must be less than 1.0e+2400" << endl;
        retval = 1;
    }
    if(maxThreads < 0){
        cerr << "ERROR: maximum number of threads must be a positive number" << endl;
        retval = 1;
    }
    //-----------------------------------------------------------------------------------------------------
    //Checking colors related variables
    if(palette_filename == ""){
        cerr << "ERROR: empty palette config file name" << endl;
        retval = 1;
    }
    if(colorMode < 0 || colorMode > MAX_COLOR_MODE){
        cerr << "ERROR: color mode out of range. Must be in [0," << MAX_COLOR_MODE << "]" << endl;
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
    size_t imageWidth = 1920;
    size_t imageHeight = 1080;
    size_t maxSectorSide = 256;

    int fractal_type = 0;
    long double offset_re = 0;
    long double offset_im = 0;
    bool juliaMode = false;
    long double julia_re  = 0;
    long double julia_im  = 0;
    long double bailout   = 2;
    long double scalingFactor = 1;
    int max_iter = 2000;

    int colorMode = 3;

    //Name of the output file should be without the extension. It will be added later
    string imageName = {"fractal"};
    string paletteConfig = {"palette"};
    vector<png::rgb_pixel> palette = {};

    bool verboseOutput = false;

    {
        bool errors = parseArgv(argc, argv,
                                imageWidth, imageHeight, maxSectorSide,
                                fractal_type, offset_re, offset_im, juliaMode, julia_re, julia_im,
                                scalingFactor, max_iter, bailout,
                                colorMode, paletteConfig, imageName, verboseOutput);
        if(errors) return 1;

        //If the loading of the palette returns 1, it failed
        if(load_palette(palette, paletteConfig)){
            cerr << "WARN: palette could not be loaded correctly, defaulting to b/w" << endl;
            palette.clear();
            palette.push_back(png::rgb_pixel{0, 0, 0});         //Load black
            palette.push_back(png::rgb_pixel{255, 255, 255});   //Load white
        }
    }
    if(palette.size() == 0){cerr << "THIS SHOULDN'T HAVE HAPPENED" << endl; return 1;}

    //Declaring and initializing image after parsing user input, by assigning it the value of the computed fractal image
    image<rgb_pixel> fractalImage = calc_and_color( imageWidth, imageHeight, maxSectorSide,
                                                    fractal_type, offset_re, offset_im, juliaMode, julia_re, julia_im,
                                                    scalingFactor, max_iter, bailout,
                                                    colorMode, palette, verboseOutput);

    vcout << "                                        \r";
    vcout << "Writing image..." << endl;
    imageName = imageName + ".png";
    fractalImage.write(imageName);
    vcout << "Done!" << endl;
    return EXIT_SUCCESS;
}
