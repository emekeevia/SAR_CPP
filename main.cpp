#include <iostream>
#include <chrono>

#include "base_algorithm.h"

int main() {
    auto start = std::chrono::steady_clock::now();
    SAR();
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << duration.count() << endl;
    return 0;
}
