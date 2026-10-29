#pragma once

#include <iomanip>
#include <iostream>
#include <string>

inline void printHeader(const std::string& title) {
    std::cout << "\n============================================================\n"
              << "  " << title << "\n"
              << "============================================================\n";
}

inline void printSubHeader(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n";
}
