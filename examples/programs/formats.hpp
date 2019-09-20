#include <sys/stat.h>
#include <iostream>
#include <cmath>
#include <cstring>

void print_line_minus(int width) {
  std::cout << std::string(width, '-') << std::endl;
}

void print_line_plus(int width) {
  std::cout << std::string(width, '+') << std::endl;
}

void print_line_equal(int width) {
  std::cout << std::string(width, '=') << std::endl;
}

void print_line_star(int width) {
  std::cout << std::string(width, '*') << std::endl;
}

void print_line_dot(int width) {
  std::cout << std::string(width, '.') << std::endl;
}



