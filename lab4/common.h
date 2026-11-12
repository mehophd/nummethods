#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

inline void printHeader(const std::string& title) {
    std::cout << "\n============================================================\n"
              << "  " << title << "\n"
              << "============================================================\n";
}

inline void printSubHeader(const std::string& title) {
    std::cout << "\n--- " << title << " ---\n";
}

inline std::string padRight(const std::string& s, int width) {
    int visible = 0;
    for (unsigned char ch : s)
        if ((ch & 0xC0) != 0x80) ++visible;  // байты-не-продолжения
    return s + std::string(width > visible ? width - visible : 0, ' ');
}

inline std::string tableHeader(
    const std::vector<std::pair<std::string, int>>& labels) {
    std::string out = " ";
    for (std::size_t i = 0; i < labels.size(); ++i) {
        out += padRight(labels[i].first, labels[i].second);
        if (i + 1 < labels.size()) out += "| ";
    }
    return out;
}

inline std::string tableSep(
    const std::vector<std::pair<std::string, int>>& labels) {
    std::string out = " ";
    for (std::size_t i = 0; i < labels.size(); ++i) {
        out += std::string(labels[i].second, '-');
        if (i + 1 < labels.size()) out += "+ ";
    }
    return out;
}

inline std::string fCell(double v, int width = 14, int prec = 9) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(prec) << std::setw(width) << v;
    return o.str();
}

inline std::string eCell(double v, int width = 11) {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << std::setw(width) << v;
    return o.str();
}

inline std::string xCell(double v, int width = 5) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(2) << std::setw(width) << v;
    return o.str();
}
