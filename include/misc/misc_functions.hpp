#ifndef MISC_FUNCTIONS_HPP_INCLUDED
#define MISC_FUNCTIONS_HPP_INCLUDED

#include <string>
#include "structs.hpp"

//Conversion functions
std::string ftype_to_string(const ftype& f);
std::string coloring_mode_to_strign(const coloring_mode& c);

//printinfo prototype
void print_info(const imagesettings_t& isettings,
                const fractalsettings_t& fsettings,
                const colorsettings_t& csettings,
                const rendersettings_t& rsettings,
                const consolesettings_t& consettings);

#endif // MISC_FUNCTIONS_HPP_INCLUDED