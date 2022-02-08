#include "option_parsing.hpp"
#include "structs.hpp"
#include "help.hpp"
#include "misc_functions.hpp"

#include <list>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include <iostream>

using namespace std;

int string_to_st(const string& str, size_t& ull){
    size_t tmp;
    try{
        tmp = stoull(str);
    }
    catch(...){
        print_error("couldn't convert \"" + str + "\" to size_t");
        return 1;
    }

    ull = tmp;
    return 0;
}

int string_to_st(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, size_t& ull)
    {return (it == vec.end() ? 2 : string_to_st((*it), ull));}

int string_to_ld(const std::string& str, long double& ld){
    long double tmp;
    try{
        tmp = stold(str);
    }
    catch(...){
        print_error("couldn't convert \"" + str + "\" to long double");
        return 1;
    }

    ld = tmp;
    return 0;
}

int string_to_ld(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, long double& ld)
    {return (it == vec.end() ? 2 : string_to_ld(*it, ld));}

int string_to_int(const string& str, int& i){
    int tmp;
    try{
        tmp = stoi(str);
    }
    catch(...){
        print_error("couldn't convert \"" + str + "\" to int");
        return 1;
    }

    i = tmp;
    return 0;
}

int string_to_int(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, int& i)
    {return (it == vec.end() ? 2 : string_to_int(*it, i));}

int extract_n_numbers_from_vec(const vector<string>& stringvec, const size_t& n, vector<long double>& extracted_numbers){
    if(stringvec.size() < n){
        print_error("can't extract " + to_string(n) + " arguments from string list");
        return 1;
    }

    vector<long double> tmp;
    for(size_t i = 0; i < n; ++i){
        tmp.emplace_back();
        if(string_to_ld(stringvec[i], tmp.back()))
            return 2;
    }

    extracted_numbers = tmp;
    return 0;
}

void print_error(const string& message, ostream& os){
    os << "[ERROR]: " << message << "\n";
}

void print_warning(const string& message, ostream& os){
    os << "[WARN]: " << message << "\n";
}

//Function to parse a list of options, passed to the programs either from argv of by reading a configuration file
//
//Return values:
//-1 -> ERR :   unrecognized option
// 0 -> OK  :   all the parsing happened succesfully
// 1 -> OK  :   help was printed
// 2 -> ERR :   generic error
int parse_options(  vector<string> options,
                    imagesettings_t& isettings, fractalsettings_t& fsettings, colorsettings_t& csettings,
                    rendersettings_t& rsettings, consolesettings_t& consettings
                 ){
    
    //---------------------------------------------------------------------------------
    //START PARSING
    while(!options.empty()){
        //Load the first element in "options", which is the first string to be parsed
        const auto front_element = options.front();
        cmdline_option current_option = cmdline_option::unknown;

        //If the string is recognized as a valid flag to set a certain option, convert it to the relatice cmdline_option
        if(map_str_to_cmdlineopt.contains(front_element))
            current_option = map_str_to_cmdlineopt.at(front_element);

        //Variable used to know how many elements to pop from "options" after the parsing of the first element is finished
        size_t elements_to_pop = map_cmdlineopt_num_elem_to_pop.at(current_option);

        //Main switch case to differentiate the behaviour for the different options
        switch(current_option){
            //---------------------------------------------------------------------
            case cmdline_option::print_help:
                print_help(cout);
                return 1;
                break;
            
            //---------------------------------------------------------------------
            case cmdline_option::enable_verbose:
                consettings.verbose_output = true;
                break;

            //---------------------------------------------------------------------
            case cmdline_option::set_width:
            {   size_t tmp_width;
                if(string_to_st(options, options.begin() + 1, tmp_width)){
                    print_error("unspecified/specified image width is invalid");
                    return 2;
                }
                else
                    isettings.image_width = tmp_width;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_height:
            {   size_t tmp_height;
                if(string_to_st(options, options.begin() + 1, tmp_height)){
                    print_error("unspecified/specified image height is invalid");
                    return 2;
                }
                else
                    isettings.image_height = tmp_height;
            }  break;

            //---------------------------------------------------------------------
            case cmdline_option::set_output_image_filename:
                if(options.size() < 2){
                    print_error("unspecified/specified output image filename is invalid");
                    return 2;
                }
                else
                    isettings.image_name = *(options.begin() + 1);
                break;

            //---------------------------------------------------------------------
            case cmdline_option::set_sector_size:
            {   size_t tmp_secsize;
                if(string_to_st(options, options.begin() + 1, tmp_secsize)){
                    print_error("unspecified/specified sector size is invalid");
                    return 2;
                }
                else
                    rsettings.max_sector_size = tmp_secsize;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_max_threads:
            {   size_t tmp_max_threads;
                if(string_to_st(options, options.begin() + 1, tmp_max_threads)){
                    print_error("unspecified/specified maximum number of threads is invalid");
                    return 2;
                }
                else
                    rsettings.max_threads = tmp_max_threads;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_fractal:
            {   ftype tmp_fractal_type = ftype::unknown;
                if(options.size() < 2){
                    print_error("not enought arguments have been provided to set the fractal");
                    return 2;
                }

                const string tmp_fractal_str = *(options.begin() + 1);
                if(map_string_to_ftype.contains(tmp_fractal_str))
                    tmp_fractal_type = map_string_to_ftype.at(tmp_fractal_str);
                else{
                    print_error("unspecified/specified fractal is invalid");
                    return 2;
                }

                //If the selected fractal is one of the testing ones or is unknown a warning is printed on screen
                if(tmp_fractal_type == ftype::test0 || tmp_fractal_type == ftype::test1 || tmp_fractal_type == ftype::test2 ||
                   tmp_fractal_type == ftype::test3 || tmp_fractal_type == ftype::test4 || tmp_fractal_type == ftype::test5 ||
                   tmp_fractal_type == ftype::test6 || tmp_fractal_type == ftype::test7 || tmp_fractal_type == ftype::test8 ||
                   tmp_fractal_type == ftype::test9)
                    print_warning("selected fractal is a test fractal");

                fsettings.fractal_type = tmp_fractal_type;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_real_offset:
            {   long double tmp_real;
                if(string_to_ld(options, options.begin() + 1, tmp_real)){
                    print_error("unspecified/specified real part for offset from the center is invalid");
                    return 2;
                }
                else
                    fsettings.offset_re = tmp_real;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_imag_offset:
            {   long double tmp_imag;
                if(string_to_ld(options, options.begin() + 1, tmp_imag)){
                    print_error("unspecified/specified imaginary part for offset from the center is invalid");
                    return 2;
                }
                else
                    fsettings.offset_im = tmp_imag;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_scaling_factor:
            {   long double tmp_scale;
                if(string_to_ld(options, options.begin() + 1, tmp_scale)){
                    print_error("unspecified/specified scaling factor is invalid");
                    return 2;
                }
                else
                    fsettings.scaling_factor = tmp_scale;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_bailout_radius:
            {   long double tmp_bailout;
                if(string_to_ld(options, options.begin() + 1, tmp_bailout)){
                    print_error("unspecified/specified bailout radius is invalid");
                    return 2;
                }
                else
                    fsettings.bailout = tmp_bailout;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::enable_julia_mode:
                fsettings.julia_mode = true;
                break;

            //---------------------------------------------------------------------
            case cmdline_option::set_julia_real_offset:
            {   long double tmp_real;
                if(string_to_ld(options, options.begin() + 1, tmp_real)){
                    print_error("unspecified/specified real part for offset of Julia mode is invalid");
                    return 2;
                }
                else
                    fsettings.julia_re = tmp_real;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_julia_im_offset:
            {   long double tmp_imag;
                if(string_to_ld(options, options.begin() + 1, tmp_imag)){
                    print_error("unspecified/specified imaginary part for offset of Julia mode is invalid");
                    return 2;
                }
                else
                    fsettings.julia_im = tmp_imag;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_max_iterations:
            {   size_t tmp_max_iter;
                if(string_to_st(options, options.begin() + 1, tmp_max_iter)){
                    print_error("unspecified/specified maximum number of iterations is invalid");
                    return 2;
                }
                else
                    fsettings.max_iter = tmp_max_iter;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::set_palette_filename:
                if(options.size() < 2)
                    return 2;

                csettings.palette_filename = *(options.begin() + 1);
                break;

            //---------------------------------------------------------------------
            case cmdline_option::set_coloring_mode:
            {   coloring_mode tmp_cmode = coloring_mode::unknown;
                if(options.size() < 2){
                    print_error("not enought arguments have been provided to set the coloring mode");
                    return 2;
                }

                const string tmp_cmode_str = *(options.begin() + 1);
                if(map_string_to_coloring_mode.contains(tmp_cmode_str))
                    tmp_cmode = map_string_to_coloring_mode.at(tmp_cmode_str);
                else{
                    print_error("unspecified/specified coloring mode is invalid");
                    return 2;
                }

                csettings.cmode = tmp_cmode;
            }   break;

            //---------------------------------------------------------------------
            case cmdline_option::enable_crosshair:
                csettings.draw_crosshair = true;
                break;

            //---------------------------------------------------------------------
            case cmdline_option::unknown:
            //default:   this is omitted so that if an option is not switched the compiler gives a warning
                print_error("unrecognized flag \"" + front_element + "\"");
                return -1;
                break;
        }

        //Pop the correct number of elements from the front of the vector
        options.erase(options.begin(), options.begin() + elements_to_pop);
    }

    return 0;
}