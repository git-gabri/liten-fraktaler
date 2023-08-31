#ifndef MATH_REG_MACHINE_HPP_INCLUDED
#define MATH_REG_MACHINE_HPP_INCLUDED

#define MRM_DEBUG_PRINT 0

#include <iostream>
#include <vector>
#include <array>
#include <complex>

namespace mrm{
    //---------------------------------------------------------------------------------------
    //Opcodes and instructions
    enum class opcode{add, sub, mul, div,
                      swap,
                      inv,    invre,  invim,
                      flip,   flipre, flipim,
                      abs,    absre,  absim,
                              nullre, nullim,
                      sin,    cos,    tan,
                      sinh,   cosh,   tanh,
                      pow,    pow2,   pow3,   log,
                      copy_ctr,   copy_ctrre,     copy_ctrim,     //Copy constant to register
                      copy_rtr,   copy_rtrre,     copy_rtrim      //Copy register to register
                    };

    struct instruction{
        mrm::opcode op;
        uint8_t result_reg;
        uint8_t op1_reg;
        uint8_t op2_reg;

        instruction(mrm::opcode _op, uint8_t res, uint8_t op1, uint8_t op2) : op(_op), result_reg(res), op1_reg(op1), op2_reg(op2) {}
        ~instruction() = default;
    };

    //---------------------------------------------------------------------------------------
    //Mathematical Register Machine
    struct machine{
        //Methods
        machine() = default;
        ~machine() = default;

        void load_script(const std::vector<mrm::instruction>& script){this->script = script;}
        std::complex<long double> execute_script(const std::vector<std::complex<long double>>& constants);

        //Members
        std::array<std::complex<long double>, 16> registers;
        std::vector<mrm::instruction> script;
    };

    //---------------------------------------------------------------------------------------
    //Text to script conversion
    int text_to_script(std::istream& input_stream, std::vector<std::complex<long double>>& output_constants, std::vector<instruction>& output_script);
}

#endif