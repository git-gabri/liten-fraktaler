#ifndef MISC_FUNCTIONS_HPP_INCLUDED
#define MISC_FUNCTIONS_HPP_INCLUDED

#include <iostream>
#include <string>
#include "structs.hpp"

//Conversion functions
std::string ftype_to_string(const ftype& f);
std::string coloring_mode_to_strign(const coloring_mode& c);

void print_error(const std::string& message, std::ostream& os = std::cerr);
void print_warning(const std::string& message, std::ostream& os = std::cerr);

#endif // MISC_FUNCTIONS_HPP_INCLUDED