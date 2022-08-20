#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <iostream>
#include <chrono>

class stopwatch{
    private:
        std::chrono::high_resolution_clock::time_point tic_time;
        std::chrono::high_resolution_clock::time_point toc_time;
        double elapsed_time;

    public:
        std::chrono::high_resolution_clock::time_point tic() {
            tic_time = std::chrono::high_resolution_clock::now();
            return tic_time;
        }

        std::chrono::high_resolution_clock::time_point toc() {
            toc_time = std::chrono::high_resolution_clock::now();
            elapsed_time = std::chrono::duration<double>(toc_time - tic_time).count();
            return toc_time;
        }

        double elapsed() const {return elapsed_time;}

        friend std::ostream& operator<<(std::ostream& os, const stopwatch& sw){
            os << sw.elapsed_time;
            return os;
        }
};

#endif