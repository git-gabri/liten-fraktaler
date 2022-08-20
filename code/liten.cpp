#include "liten.hpp"

#include <thread>

imagesettings_t             lf::internals::isettings{};                 //struct containing various image related settings
fractalsettings_t           lf::internals::fsettings{};                 //struct containing various fractal related settings
colorsettings_t             lf::internals::csettings{};                 //struct containing various color related settings
rendersettings_t            lf::internals::rsettings{};                 //struct containing various information on how to render the image
consolesettings_t           lf::internals::consettings{};               //struct containing various information on what the program should output on the console
std::vector<png::rgb_pixel> lf::internals::palette = {};                //color palette

//Implementing methods of the namespace "lf"
using namespace lf::internals;

//Initialization of the namespace
void lf::init(){
    rsettings.max_threads = std::thread::hardware_concurrency();
    if(rsettings.max_threads == 0) rsettings.max_threads = FALLBACK_NUM_THREADS;
}
