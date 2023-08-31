#include "mrm.hpp"
#include <array>
#include <string>
#include <sstream>
#include <regex>

//--------------------------------------------------------------------------------------------------------------------
//Mathematical Register Machine
std::complex<long double> mrm::machine::execute_script(const std::vector<std::complex<long double>>& constants){
    for(auto pc = script.begin(); pc != script.end(); ++pc){
        const mrm::instruction current_instruction = *pc;

        switch(current_instruction.op){
            case mrm::opcode::add:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] + registers[current_instruction.op2_reg];
                break;
            case mrm::opcode::sub:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] - registers[current_instruction.op2_reg];
                break;
            case mrm::opcode::mul:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] * registers[current_instruction.op2_reg];
                break;
            case mrm::opcode::div:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] / registers[current_instruction.op2_reg];
                break;
            case mrm::opcode::swap:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].imag());
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].real());
                break;
            case mrm::opcode::inv:
                registers[current_instruction.result_reg] = std::complex<long double>{1, 0} / registers[current_instruction.op1_reg];
                break;
            case mrm::opcode::invre:
                registers[current_instruction.result_reg].real(1.0l / registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::invim:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(1.0l / registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::flip:
                registers[current_instruction.result_reg] = -1.0l * registers[current_instruction.op1_reg];
                break;
            case mrm::opcode::flipre:
                registers[current_instruction.result_reg].real(-1.0l * registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::flipim:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(-1.0l * registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::abs:
                registers[current_instruction.result_reg] = std::abs(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::absre:
                registers[current_instruction.result_reg].real(std::abs(registers[current_instruction.op1_reg].real()));
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::absim:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(std::abs(registers[current_instruction.op1_reg].imag()));
                break;
            case mrm::opcode::nullre:
                registers[current_instruction.result_reg].real(0.0l);
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::nullim:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].real());
                registers[current_instruction.result_reg].imag(0.0l);
                break;
            case mrm::opcode::sin:
                registers[current_instruction.result_reg] = std::sin(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::cos:
                registers[current_instruction.result_reg] = std::cos(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::tan:
                registers[current_instruction.result_reg] = std::tan(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::sinh:
                registers[current_instruction.result_reg] = std::sinh(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::cosh:
                registers[current_instruction.result_reg] = std::cosh(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::tanh:
                registers[current_instruction.result_reg] = std::tanh(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::pow:
                registers[current_instruction.result_reg] = std::pow(registers[current_instruction.op1_reg], registers[current_instruction.op2_reg]);
                break;
            case mrm::opcode::pow2:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] * registers[current_instruction.op1_reg];
                break;
            case mrm::opcode::pow3:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg] * registers[current_instruction.op1_reg] * registers[current_instruction.op1_reg];
                break;
            case mrm::opcode::log:
                registers[current_instruction.result_reg] = std::log(registers[current_instruction.op1_reg]);
                break;
            case mrm::opcode::copy_ctr:
                registers[current_instruction.result_reg] = constants[current_instruction.op1_reg];
                break;
            case mrm::opcode::copy_ctrre:
                registers[current_instruction.result_reg].real(constants[current_instruction.op1_reg].real());
                break;
            case mrm::opcode::copy_ctrim:
                registers[current_instruction.result_reg].imag(constants[current_instruction.op1_reg].imag());
                break;
            case mrm::opcode::copy_rtr:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg];
                break;
            case mrm::opcode::copy_rtrre:
                registers[current_instruction.result_reg].real(registers[current_instruction.op1_reg].real());
                break;
            case mrm::opcode::copy_rtrim:
                registers[current_instruction.result_reg].imag(registers[current_instruction.op1_reg].imag());
                break;
        }
    }

    return registers[0];
}

//--------------------------------------------------------------------------------------------------------------------
//Text to script conversion
int mrm::text_to_script(std::istream& input_stream, std::vector<std::complex<long double>>& output_constants, std::vector<instruction>& output_script){
    std::vector<std::complex<long double>> consts;
    std::vector<mrm::instruction> script;

    //Regexes
    std::regex constant_regex   (R"foo(^c\d+ = \S+$)foo");
    std::regex add_regex        (R"foo(^\d{1,2} = \d{1,2} \+ \d{1,2}$)foo");
    std::regex sub_regex        (R"foo(^\d{1,2} = \d{1,2} \- \d{1,2}$)foo");
    std::regex mul_regex        (R"foo(^\d{1,2} = \d{1,2} \* \d{1,2}$)foo");
    std::regex div_regex        (R"foo(^\d{1,2} = \d{1,2} \/ \d{1,2}$)foo");
    std::regex swap_regex       (R"foo(^\d{1,2} = swap \d{1,2}$)foo");
    std::regex inv_regex        (R"foo(^\d{1,2} = inv \d{1,2}$)foo");
    std::regex invre_regex      (R"foo(^\d{1,2} = invre \d{1,2}$)foo");
    std::regex invim_regex      (R"foo(^\d{1,2} = invim \d{1,2}$)foo");
    std::regex flip_regex       (R"foo(^\d{1,2} = flip \d{1,2}$)foo");
    std::regex flipre_regex     (R"foo(^\d{1,2} = flipre \d{1,2}$)foo");
    std::regex flipim_regex     (R"foo(^\d{1,2} = flipim \d{1,2}$)foo");
    std::regex abs_regex        (R"foo(^\d{1,2} = abs \d{1,2}$)foo");
    std::regex absre_regex      (R"foo(^\d{1,2} = absre \d{1,2}$)foo");
    std::regex absim_regex      (R"foo(^\d{1,2} = absim \d{1,2}$)foo");
    std::regex nullre_regex     (R"foo(^\d{1,2} = nullre \d{1,2}$)foo");
    std::regex nullim_regex     (R"foo(^\d{1,2} = nullim \d{1,2}$)foo");
    std::regex sin_regex        (R"foo(^\d{1,2} = sin \d{1,2}$)foo");
    std::regex cos_regex        (R"foo(^\d{1,2} = cos \d{1,2}$)foo");
    std::regex tan_regex        (R"foo(^\d{1,2} = tan \d{1,2}$)foo");
    std::regex sinh_regex       (R"foo(^\d{1,2} = sinh \d{1,2}$)foo");
    std::regex cosh_regex       (R"foo(^\d{1,2} = cosh \d{1,2}$)foo");
    std::regex tanh_regex       (R"foo(^\d{1,2} = tanh \d{1,2}$)foo");
    std::regex pow_regex        (R"foo(^\d{1,2} = \d{1,2} pow \d{1,2}$)foo");
    std::regex pow2_regex       (R"foo(^\d{1,2} = pow2 \d{1,2}$)foo");
    std::regex pow3_regex       (R"foo(^\d{1,2} = pow3 \d{1,2}$)foo");
    std::regex log_regex        (R"foo(^\d{1,2} = log \d{1,2}$)foo");
    std::regex log2_regex       (R"foo(^\d{1,2} = log2 \d{1,2}$)foo");
    std::regex copy_ctr_regex   (R"foo(^\d{1,2} = c\d+$)foo");
    std::regex copy_ctrre_regex (R"foo(^\d{1,2} = c\d+re$)foo");
    std::regex copy_ctrim_regex (R"foo(^\d{1,2} = c\d+im$)foo");
    std::regex copy_rtr_regex   (R"foo(^\d{1,2} = \d{1,2}$)foo");
    std::regex copy_rtrre_regex (R"foo(^\d{1,2} = \d{1,2}re$)foo");
    std::regex copy_rtrim_regex (R"foo(^\d{1,2} = \d{1,2}im$)foo");

    //Read and parse input line by line
    for(std::string line; std::getline(input_stream, line); ){
        //------------------------------------------------------------------
        //Split string into substrings separated by spaces
        std::vector<std::string> substrings;

        bool was_prev_char_space = true;
        for(size_t c = 0; c != line.length(); ++c){
            const bool is_curr_char_space = std::isspace(line[c]);

            //Beginning of new word: prev char is space, this char is not space
            if(was_prev_char_space && !is_curr_char_space){
                substrings.emplace_back(1, line[c]);
            }
            //Continuation of new word: prev char is not space, this char is not space
            if(!was_prev_char_space && !is_curr_char_space){
                substrings.back() += line[c];
            }
            //End of a word: prev char is not space, this char is space
            //Nothig to do
            //Continuation of a space: prev char is space, this char is space
            //Nothig to do

            //Update of boolean flag
            was_prev_char_space = is_curr_char_space;
        }

#if MRM_DEBUG_PRINT
        for(auto it = substrings.begin(); it != substrings.end(); ++it)
            std::cout << "\"" << *it << "\", ";
#endif

        //------------------------------------------------------------------
        //Convert the text information in constants and operations
        if(std::regex_match(line, constant_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "LOAD CONST" << std::endl;
#endif
            //sub[0] -> 'c' + constant id
            //sub[1] -> "="
            //sub[2] -> constant value
            substrings[0].erase(0, 1); //remove initial 'c'
            const auto const_id = std::stoul(substrings[0]);

            //Check if id is out of range of 8 bits
            if(const_id > 255) return 2;
            //If valid, resize
            consts.resize(std::max(consts.size(), const_id + 1));

            //Parse constant value
            std::complex<long double> const_value{};
            std::stringstream ss_sub2(substrings[2]);
            ss_sub2 >> const_value;
            consts[const_id] = const_value;
        } else
        if(std::regex_match(line, add_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "ADD" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> op1 reg
            //sub[3] -> "+"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::add, res, op1, op2);
        } else
        if(std::regex_match(line, sub_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "SUB" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> op1 reg
            //sub[3] -> "-"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::sub, res, op1, op2);
        } else
        if(std::regex_match(line, mul_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "MUL" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> op1 reg
            //sub[3] -> "*"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::mul, res, op1, op2);
        } else
        if(std::regex_match(line, div_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "DIV" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> op1 reg
            //sub[3] -> "/"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::div, res, op1, op2);
        } else
        if(std::regex_match(line, swap_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "SWAP REAL IMAGINARY" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "swap"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::swap, res, op1, op2);
        } else
        if(std::regex_match(line, inv_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "INVERT" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "inv"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::inv, res, op1, op2);
        } else
        if(std::regex_match(line, invre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "INVERT REAL" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "invre"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::invre, res, op1, op2);
        } else
        if(std::regex_match(line, invim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "INVERT IMAGINARY" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "invim"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::invim, res, op1, op2);
        } else
        if(std::regex_match(line, flip_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "FLIP SIGN" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "flip"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::flip, res, op1, op2);
        } else
        if(std::regex_match(line, flipre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "FLIP SIGN REAL" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "flipre"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::flipre, res, op1, op2);
        } else
        if(std::regex_match(line, flipim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "FLIP SIGN IMAGINARY" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "flipim"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::flipim, res, op1, op2);
        } else
        if(std::regex_match(line, abs_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "ABS" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "abs"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::abs, res, op1, op2);
        } else
        if(std::regex_match(line, absre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "ABS REAL" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "absre"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::absre, res, op1, op2);
        } else
        if(std::regex_match(line, absim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "ABS IMAGINARY" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "absim"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::absim, res, op1, op2);
        } else
        if(std::regex_match(line, nullre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "NULLIFY REAL" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "nullre"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::nullre, res, op1, op2);
        } else
        if(std::regex_match(line, nullim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "NULLIFY IMAGINARY" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "nullim"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::nullim, res, op1, op2);
        } else
        if(std::regex_match(line, sin_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "SIN" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "sin"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::sin, res, op1, op2);
        } else
        if(std::regex_match(line, cos_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COS" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "cos"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::cos, res, op1, op2);
        } else
        if(std::regex_match(line, tan_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "TAN" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "tan"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::tan, res, op1, op2);
        } else
        if(std::regex_match(line, sinh_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "SINH" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "sinh"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::sinh, res, op1, op2);
        } else
        if(std::regex_match(line, cosh_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COSH" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "cosh"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::cosh, res, op1, op2);
        } else
        if(std::regex_match(line, tanh_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "TANH" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "tanh"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::tanh, res, op1, op2);
        } else
        if(std::regex_match(line, pow_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "POW" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> op1 reg
            //sub[3] -> "pow"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::pow, res, op1, op2);
        } else
        if(std::regex_match(line, pow2_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "POW2" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "pow2"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::pow2, res, op1, op2);
        } else
        if(std::regex_match(line, pow3_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "POW3" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "pow3"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::pow3, res, op1, op2);
        } else
        if(std::regex_match(line, log_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "LOG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> "log"
            //sub[3] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[3]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::log, res, op1, op2);
        } else
        if(std::regex_match(line, copy_ctr_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY CONST TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> 'c' + constant id
            const auto res = std::stoul(substrings[0]);
            substrings[2].erase(0, 1);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_ctr, res, op1, op2);
        } else
        if(std::regex_match(line, copy_ctrre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY CONST REAL TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> 'c' + constant id + 're'
            const auto res = std::stoul(substrings[0]);
            substrings[2].erase(0, 1);
            substrings[2].pop_back();
            substrings[2].pop_back();
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_ctrre, res, op1, op2);
        } else
        if(std::regex_match(line, copy_ctrim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY CONST IMAGINARY TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> 'c' + constant id + 're'
            const auto res = std::stoul(substrings[0]);
            substrings[2].erase(0, 1);
            substrings[2].pop_back();
            substrings[2].pop_back();
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_ctrim, res, op1, op2);
        } else
        if(std::regex_match(line, copy_rtr_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY REG TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> src reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_rtr, res, op1, op2);
        } else
        if(std::regex_match(line, copy_rtrre_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY REG REAL TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> src reg + 're'
            substrings[2].pop_back();
            substrings[2].pop_back();
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_rtrre, res, op1, op2);
        } else
        if(std::regex_match(line, copy_rtrim_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY REG IMAGINARY TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "="
            //sub[2] -> src reg + 'im'
            substrings[2].pop_back();
            substrings[2].pop_back();
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_rtrim, res, op1, op2);
        } else
        {
#if MRM_DEBUG_PRINT
            std::cout << "ERROR!!!!" << std::endl;
#endif
            return 1;
        }
    }

    output_constants = consts;
    output_script = script;

    return 0;
}