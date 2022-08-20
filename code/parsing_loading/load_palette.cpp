#include "liten.hpp"
#include "misc_functions.hpp"

#include <fstream>

using namespace std;
using namespace lf::internals;

//Function to load the color palette from file
//Returns 0 if the loading was succesful, 1 if some errors occurred
/*The config file of the color palette should have the following format
* uint8 uint8 uint8
* uint8 uint8 uint8
* uint8 uint8 uint8
* [...]
*
* Each line represents a color to be loaded in the color palette, the first row corresponds
* to the RED channel, the middle row to the GREEN channel and the last row to the BLUE channel
*/
int lf::load_palette(){
    //Clear color palette
    palette.clear();

    //Open input text file
    ifstream input_file(csettings.palette_filename, ifstream::in);
    if(!input_file.is_open()){
        print_warning("palette could not be loaded correctly, defaulting to b/w");
        palette.clear();
        palette.push_back(png::rgb_pixel{0, 0, 0});         //Load black
        palette.push_back(png::rgb_pixel{255, 255, 255});   //Load white
        return 1;
    }

    //Buffer to read the file line by line
    string line;
    while(getline(input_file, line)){
        //String stream to parse the line
        istringstream iss(line);
        //Temporary integers to store the values of the color, because the character 0 in the file gets interpreted as color 48,
        //so we have to store the reading from the file in an integer, so that it gets interpreted as a number
        int tmpred = 0, tmpgreen = 0, tmpblue = 0;

        //Read from iss and check that we read 3 color values. If not, return 1
        if(!(iss >> tmpred >> tmpgreen >> tmpblue)) return 1;

        palette.push_back(png::rgb_pixel(uint8_t(tmpred), uint8_t(tmpgreen), uint8_t(tmpblue)));
    }

    return 0;
}