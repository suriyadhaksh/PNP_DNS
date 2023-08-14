#pragma once
#include <talyfem/talyfem.h>

#ifndef UTIL_FUNCTIONS_HPP
#define UTIL_FUNCTIONS_HPP

#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

inline std::string Suffix(std::string name, unsigned i) {
  std::ostringstream tmp;
  tmp << '_' << std::setw(4) << std::setfill('0') << i;
  std::string::size_type p = name.rfind('.');

  if (p != std::string::npos)
    name.insert(p, tmp.str());
  else
    name += tmp.str();

  return name;
}

inline std::string Suffix(std::string name, unsigned r, unsigned i) {
  std::ostringstream tmp;
  tmp << '_' << std::setw(3) << std::setfill('0') << r;
  tmp << '.' << std::setw(4) << std::setfill('0') << i;
  std::string::size_type p = name.rfind('.');

  if (p != std::string::npos)
    name.insert(p, tmp.str());
  else
    name += tmp.str();

  return name;
}

#endif
