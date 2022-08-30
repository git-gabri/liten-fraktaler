#ifndef MISC_FUNCTIONS_HPP_INCLUDED
#define MISC_FUNCTIONS_HPP_INCLUDED

#include <iostream>
#include <string>
#include "structs.hpp"

//Conversion functions
std::string fractal_type_to_string(const ftype& f);
std::string renderer_type_to_string(const rtype& r);
std::string coloring_mode_to_string(const coloring_mode& c);

void print_error(const std::string& message, std::ostream& os = std::cerr);
void print_warning(const std::string& message, std::ostream& os = std::cerr);

#endif // MISC_FUNCTIONS_HPP_INCLUDED