/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified: 2023-06-15 13:58
#
# Filename: tic_toc.hpp
#
# Description:
#
************************************************/

#ifndef VSLAM_TIC_TOC_HPP_
#define VSLAM_TIC_TOC_HPP_

#include <chrono>
#include <cstdlib>
#include <ctime>

namespace magnav {
class TicToc {
public:
  TicToc() { tic(); }

  void tic() { start = std::chrono::system_clock::now(); }

  double toc() {
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count() * 1000;
  }

private:
  std::chrono::time_point<std::chrono::system_clock> start, end;
};
} // namespace magnav

#endif
