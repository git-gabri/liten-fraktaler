#ifndef PARSE_OPTIONS_HPP_INCLUDED
#define PARSE_OPTIONS_HPP_INCLUDED

#include "liten.hpp"

#include <iostream>
#include <map>

//Utility functions
//Return values:
// 0 -> OK
//!0 -> some error occurred
int string_to_st(const std::string& str, size_t& ull);
int string_to_st(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, size_t& ull);
int string_to_ld(const std::string& str, long double& ld);
int string_to_ld(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, long double& ld);
int string_to_int(const std::string& str, int& i);
int string_to_int(const std::vector<std::string>& vec, const std::vector<std::string>::iterator& it, int& i);

int extract_n_numbers_from_vec(const std::vector<std::string>& stringvec, const size_t& n, std::vector<long double>& extracted_numbers);

//-------------------------------------------------------------------------------
//Command line option definition and maps used as lookup tables
enum class cmdline_option{
    print_help,
    enable_verbose,

    set_width, set_height,
    set_output_image_filename,

    set_sector_size,
    set_max_threads,
    set_renderer,
    set_mibc_max_iter,

    set_fractal,
    set_real_offset, set_imag_offset,
    set_scaling_factor,
    set_bailout_radius,
    enable_julia_mode, set_julia_real_offset, set_julia_im_offset,
    set_max_iterations,

    set_palette_filename,
    set_coloring_mode,
    enable_crosshair,
    unknown
};

const std::map<cmdline_option, size_t> map_cmdlineopt_num_elem_to_pop{
    {cmdline_option::print_help, 1},
    {cmdline_option::enable_verbose, 1},

    {cmdline_option::set_width, 2},
    {cmdline_option::set_height, 2},
    {cmdline_option::set_output_image_filename, 2},

    {cmdline_option::set_sector_size, 2},
    {cmdline_option::set_max_threads, 2},
    {cmdline_option::set_renderer, 2},
    {cmdline_option::set_mibc_max_iter, 2},

    {cmdline_option::set_fractal, 2},
    {cmdline_option::set_real_offset, 2},
    {cmdline_option::set_imag_offset, 2},
    {cmdline_option::set_scaling_factor, 2},
    {cmdline_option::set_bailout_radius, 2},
    {cmdline_option::enable_julia_mode, 1},
    {cmdline_option::set_julia_real_offset, 2},
    {cmdline_option::set_julia_im_offset, 2},
    {cmdline_option::set_max_iterations, 2},

    {cmdline_option::set_palette_filename, 2},
    {cmdline_option::set_coloring_mode, 2},
    {cmdline_option::enable_crosshair, 1},
    {cmdline_option::unknown, 1}
};

const std::map<std::string, cmdline_option> map_str_to_cmdlineopt{
    {"-H",              cmdline_option::print_help},
    {"--help",          cmdline_option::print_help},
    {"-v",              cmdline_option::enable_verbose},
    {"--verbose",       cmdline_option::enable_verbose},

    {"-w",              cmdline_option::set_width},
    {"--width",         cmdline_option::set_width},
    {"-h",              cmdline_option::set_height},
    {"--height",        cmdline_option::set_height},
    {"-o",              cmdline_option::set_output_image_filename},
    {"--output-image-filename", cmdline_option::set_output_image_filename},

    {"-S",              cmdline_option::set_sector_size},
    {"--sector-size",   cmdline_option::set_sector_size},
    {"-T",              cmdline_option::set_max_threads},
    {"--max-threads",   cmdline_option::set_max_threads},
    {"--renderer",      cmdline_option::set_renderer},
    {"--mibc-max-iter", cmdline_option::set_mibc_max_iter},

    {"-f",              cmdline_option::set_fractal},
    {"--fractal",       cmdline_option::set_fractal},
    {"-r",              cmdline_option::set_real_offset},
    {"--real",          cmdline_option::set_real_offset},
    {"-i",              cmdline_option::set_imag_offset},
    {"--imag",          cmdline_option::set_imag_offset},
    {"-s",              cmdline_option::set_scaling_factor},
    {"--scaling-factor",cmdline_option::set_scaling_factor},
    {"--zoom",          cmdline_option::set_scaling_factor},
    {"-b",              cmdline_option::set_bailout_radius},
    {"--bailout",       cmdline_option::set_bailout_radius},
    {"-j",              cmdline_option::enable_julia_mode},
    {"--julia-mode",    cmdline_option::enable_julia_mode},
    {"-A",              cmdline_option::set_julia_real_offset},
    {"--julia-real",    cmdline_option::set_julia_real_offset},
    {"-B",              cmdline_option::set_julia_im_offset},
    {"--julia-imag",    cmdline_option::set_julia_im_offset},
    {"-t",              cmdline_option::set_max_iterations},
    {"--max-iter",      cmdline_option::set_max_iterations},

    {"-p",              cmdline_option::set_palette_filename},
    {"--palette-filename", cmdline_option::set_palette_filename},
    {"-c",              cmdline_option::set_coloring_mode},
    {"--coloring-mode", cmdline_option::set_coloring_mode},
    {"-C",              cmdline_option::enable_crosshair},
    {"--crosshair",     cmdline_option::enable_crosshair}
};

const std::map<std::string, ftype> map_string_to_ftype{
    {"mandelbrot",  ftype::mandelbrot},
    {"tippets",     ftype::tippets},
    {"burnship",    ftype::burnship},
    {"mandelbar",   ftype::mandelbar},
    {"magnet1",     ftype::magnet1},
    {"magnet2",     ftype::magnet2},
    {"cactus",      ftype::cactus},
    {"cactus0",     ftype::cactus0},
    {"zubieta",     ftype::zubieta},
    {"zubitheta",   ftype::zubitheta},
    {"logmap",      ftype::logmap},
    {"unpolsquare", ftype::unpolsquare},
    {"moth",        ftype::moth},
    {"cave",        ftype::cave},
    {"wankel",      ftype::wankel},
    {"seaangel",    ftype::seaangel},
    {"smith",       ftype::smith},
    {"spade",       ftype::spade},

    {"__test0", ftype::test0},
    {"__test1", ftype::test1},
    {"__test2", ftype::test2},
    {"__test3", ftype::test3},
    {"__test4", ftype::test4},
    {"__test5", ftype::test5},
    {"__test6", ftype::test6},
    {"__test7", ftype::test7},
    {"__test8", ftype::test8},
    {"__test9", ftype::test9}
};

const std::map<std::string, rtype> map_string_to_rtype{
    {"basic",   rtype::basic},
    {"mibc",    rtype::mibc},

    {"__test0", rtype::test0},
    {"__test1", rtype::test1},
    {"__test2", rtype::test2},
    {"__test3", rtype::test3},
    {"__test4", rtype::test4},
    {"__test5", rtype::test5},
    {"__test6", rtype::test6},
    {"__test7", rtype::test7},
    {"__test8", rtype::test8},
    {"__test9", rtype::test9}
};

const std::map<std::string, coloring_mode> map_string_to_coloring_mode{
    {"binary",      coloring_mode::binary},
    {"linear",      coloring_mode::linear},
    {"ln",          coloring_mode::ln},
    {"sqrt",        coloring_mode::sqrt},
    {"cbrt",        coloring_mode::cbrt},
    {"scurve",      coloring_mode::scurve},
    {"lastangle",   coloring_mode::lastangle}
};

#endif