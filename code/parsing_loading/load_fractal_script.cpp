#include "liten.hpp"

#include <fstream>

using namespace std;
using namespace lf::internals;

//Function to load the fractal script from filr
//Returns:
//  -1 if the file couldn't be openend
//   0 if the loading was succesful
//   1 if an invalid instruction was found
//   2 if an index of something exceeds 255
int lf::internals::load_custom_fractal_script(const string& script_filename){
    ifstream input_file(script_filename, std::ios::in);

    if(!input_file.is_open())
        return -1;

    return mrm::text_to_script(input_file, lf::internals::custom_fractal_script_constants, lf::internals::custom_fractal_script);
}