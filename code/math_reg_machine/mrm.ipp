#ifndef MATH_REG_MACHINE_IPP_INCLUDED
#define MATH_REG_MACHINE_IPP_INCLUDED

#include "mrm.hpp"
#include <array>
#include <string>
#include <sstream>
#include <regex>

//--------------------------------------------------------------------------------------------------------------------
//Opcodes and instructions
enum class mrm::opcode{ add, sub, mul, div,
                        copy_ctr,   //Copy constant to register
                        copy_rtr    //Copy register to register
                      };

struct mrm::instruction{
    mrm::opcode op;
    uint8_t result_reg;
    uint8_t op1_reg;
    uint8_t op2_reg;

    instruction(mrm::opcode _op, uint8_t res, uint8_t op1, uint8_t op2) : op(_op), result_reg(res), op1_reg(op1), op2_reg(op2) {}
    ~instruction() = default;
};

//--------------------------------------------------------------------------------------------------------------------
//Mathematical Register Machine
template<typename T>
struct mrm::machine{
    //Methods
    machine() = default;
    ~machine() = default;

    void load_script(const std::vector<mrm::instruction>& script){this->script = script;}
    T execute_script(const std::vector<T>& constants);

    //Members
    std::array<T, 16> registers;
    std::vector<mrm::instruction> script;
};

template<typename T>
T mrm::machine<T>::execute_script(const std::vector<T>& constants){
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
            case mrm::opcode::copy_ctr:
                registers[current_instruction.result_reg] = constants[current_instruction.op1_reg];
                break;
            case mrm::opcode::copy_rtr:
                registers[current_instruction.result_reg] = registers[current_instruction.op1_reg];
                break;
        }
    }

    return registers[0];
}

//--------------------------------------------------------------------------------------------------------------------
//Text to script conversion
template<typename T>
int mrm::text_to_script(std::istream& input_stream, std::vector<T>& output_constants, std::vector<instruction>& output_script){
    std::vector<T> consts;
    std::vector<mrm::instruction> script;

    //Regexes
    std::regex constant_regex   (R"foo(^c\d+ <- \S+$)foo");
    std::regex add_regex        (R"foo(^\d{1,2} <- \d{1,2} \+ \d{1,2}$)foo");
    std::regex sub_regex        (R"foo(^\d{1,2} <- \d{1,2} \- \d{1,2}$)foo");
    std::regex mul_regex        (R"foo(^\d{1,2} <- \d{1,2} \* \d{1,2}$)foo");
    std::regex div_regex        (R"foo(^\d{1,2} <- \d{1,2} \/ \d{1,2}$)foo");
    std::regex copy_ctr_regex   (R"foo(^\d{1,2} <- c\d+$)foo");
    std::regex copy_rtr_regex   (R"foo(^\d{1,2} <- \d{1,2}+$)foo");

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
            //sub[1] -> "<-"
            //sub[2] -> constant value
            substrings[0].erase(0, 1); //remove initial 'c'
            const auto const_id = std::stoul(substrings[0]);

            //Check if id is out of range of 8 bits
            if(const_id > 255) return 2;
            //If valid, resize
            consts.resize(std::max(consts.size(), const_id + 1));

            //Parse constant value
            T const_value{};
            std::stringstream ss_sub2(substrings[2]);
            ss_sub2 >> const_value;
            consts[const_id] = const_value;
        } else
        if(std::regex_match(line, add_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "ADD" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "<-"
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
            //sub[1] -> "<-"
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
            //sub[1] -> "<-"
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
            //sub[1] -> "<-"
            //sub[2] -> op1 reg
            //sub[3] -> "/"
            //sub[4] -> op2 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = std::stoul(substrings[4]);
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::div, res, op1, op2);
        } else
        if(std::regex_match(line, copy_ctr_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY CONST TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "<-"
            //sub[2] -> 'c' + constant id
            const auto res = std::stoul(substrings[0]);
            substrings[2].erase(0, 1);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_ctr, res, op1, op2);
        } else
        if(std::regex_match(line, copy_rtr_regex)){
#if MRM_DEBUG_PRINT
            std::cout << "COPY REG TO REG" << std::endl;
#endif
            //sub[0] -> result reg
            //sub[1] -> "<-"
            //sub[2] -> op1 reg
            const auto res = std::stoul(substrings[0]);
            const auto op1 = std::stoul(substrings[2]);
            const auto op2 = 0;
            if((res > 255) || (op1 > 255) | (op2 > 255)) return 2;

            script.emplace_back(mrm::opcode::copy_rtr, res, op1, op2);
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

#endif