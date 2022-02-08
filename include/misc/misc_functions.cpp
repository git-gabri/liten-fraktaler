#include "misc_functions.hpp"
#include <iostream>

using namespace std;

//Conversion functions
string ftype_to_string(const ftype& f){
    string ret_str;

    switch (f){
        case ftype::mandelbrot:     ret_str = "Mandelbrot set";     break;
        case ftype::tippets:        ret_str = "Tippets Mandel";     break;
        case ftype::burnship:       ret_str = "Burning ship";       break;
        case ftype::mandelbar:      ret_str = "Mandelbar set";      break;
        case ftype::magnet1:        ret_str = "Magnet type I";      break;
        case ftype::magnet2:        ret_str = "Magnet type II";     break;
        case ftype::cactus:         ret_str = "Cactus";             break;
        case ftype::cactus0:        ret_str = "Cactus0";            break;
        case ftype::zubieta:        ret_str = "Zubieta";            break;
        case ftype::zubitheta:      ret_str = "Zubitheta";          break;
        case ftype::logmap:         ret_str = "Logistic map";       break;
        case ftype::unpolsquare:    ret_str = "Unpolished square";  break;
        case ftype::moth:           ret_str = "Moth";               break;
        case ftype::cave:           ret_str = "Cave";               break;
        case ftype::wankel:         ret_str = "Wankel engine";      break;
        case ftype::seaangel:       ret_str = "Sea Angel";          break;
        case ftype::smith:          ret_str = "Smith";              break;
        case ftype::spade:          ret_str = "Spade fractal";      break;

        case ftype::test0:          ret_str = "--TEST 0--";         break;
        case ftype::test1:          ret_str = "--TEST 1--";         break;
        case ftype::test2:          ret_str = "--TEST 2--";         break;
        case ftype::test3:          ret_str = "--TEST 3--";         break;
        case ftype::test4:          ret_str = "--TEST 4--";         break;
        case ftype::test5:          ret_str = "--TEST 5--";         break;
        case ftype::test6:          ret_str = "--TEST 6--";         break;
        case ftype::test7:          ret_str = "--TEST 7--";         break;
        case ftype::test8:          ret_str = "--TEST 8--";         break;
        case ftype::test9:          ret_str = "--TEST 9--";         break;

        default:
            ret_str = "Unknown/undocumented";
            break;
    }

    return ret_str;
}

string coloring_mode_to_strign(const coloring_mode& c){
    string ret_str;

    switch (c){
        case coloring_mode::binary:     ret_str = "binary (b/w)";   break;
        case coloring_mode::linear:     ret_str = "linear";         break;
        case coloring_mode::ln:         ret_str = "smooth (ln)";    break;
        case coloring_mode::sqrt:       ret_str = "sqrt(iter)";     break;
        case coloring_mode::cbrt:       ret_str = "cbrt(iter)";     break;
        case coloring_mode::scurve:     ret_str = "S curve";        break;
        case coloring_mode::lastangle:  ret_str = "atan2(im, re)";  break;
        default:
            ret_str = "Unknown/undocumented";
            break;
    }

    return ret_str;
}

//Function to print info about the current settings of the program
void print_info(const imagesettings_t& isettings,
                const fractalsettings_t& fsettings,
                const colorsettings_t& csettings,
                const rendersettings_t& rsettings,
                const consolesettings_t& consettings
              ){

    string fractal_name = ftype_to_string(fsettings.fractal_type);

    const string coloring_mode_name = coloring_mode_to_strign(csettings.cmode);

    cout << "Rendering      : " << fractal_name << (fsettings.julia_mode ? " (J)" : "") << " (" << rsettings.max_threads << " thread(s))" << endl;
    cout << "Image size     : " << isettings.image_width << " x " << isettings.image_height << endl;
    cout << "Sectors up to  : " << rsettings.max_sector_size << "x" << rsettings.max_sector_size << " = " << rsettings.max_sector_size * rsettings.max_sector_size << " pixels^2" << endl;

    cout << "Iterations     : " << fsettings.max_iter << endl;
    cout << "Bailout radius : " << fsettings.bailout << endl;
    cout << "Coloring mode  : " << coloring_mode_name << endl;
    cout << "Centered at    : (" << fsettings.offset_re << ", " << fsettings.offset_im << " i)" << endl;
    if(fsettings.julia_mode) cout << "Julia c const. : (" << fsettings.julia_re << ", " << fsettings.julia_im << " i)" << endl;

    //Auxiliary variable, used when calculating scaling from pixel space to complex plane, to prevent stretching
    const auto square_scale = min(isettings.image_width, isettings.image_height);
    //Top left corner
    const long double start_re = 4 * (- (long double)isettings.image_width / 2)  / ((long double)square_scale * fsettings.scaling_factor) + fsettings.offset_re;
    const long double start_im = 4 * (- (long double)isettings.image_height / 2) / ((long double)square_scale * fsettings.scaling_factor) + fsettings.offset_im;
    //Bottom right
    const long double stop_re  = 4 * ((long double)isettings.image_width / 2)  / ((long double)square_scale * fsettings.scaling_factor) + fsettings.offset_re;
    const long double stop_im  = 4 * ((long double)isettings.image_height / 2) / ((long double)square_scale * fsettings.scaling_factor) + fsettings.offset_im;
    cout << "Span  : (" << stop_re - start_re << " x " << stop_im - start_im << ")" << endl;
    cout << "Scale : " << fsettings.scaling_factor << endl;
}