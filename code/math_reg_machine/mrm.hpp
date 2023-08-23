#ifndef MATH_REG_MACHINE_HPP_INCLUDED
#define MATH_REG_MACHINE_HPP_INCLUDED

#define MRM_DEBUG_PRINT 0

#include <iostream>
#include <vector>

namespace mrm{
    enum class opcode;

    struct instruction;

    template<typename T>
    struct machine;

    template<typename T>
    int text_to_script(std::istream& input_stream, std::vector<T>& output_constants, std::vector<instruction>& output_script);
}

#include "mrm.ipp"

#endif