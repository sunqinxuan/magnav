/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-15 14:04
#
# Filename:		tolles_lawson.hpp
#
# Description:
#
************************************************/

#ifndef MAGNAV_MAGCMP_TL_HPP_
#define MAGNAV_MAGCMP_TL_HPP_

#include <ceres/ceres.h>
#include <eigen3/Eigen/Dense>
#include <unordered_set>
#include <vector>

#include "filter/filter.hpp"
#include "filter/butterworth.hpp"
#include "tools/message_print.hpp"

//#include "optimizer/camera_projection_factor.hpp"
//#include "optimizer/imu_factor.hpp"
//#include "optimizer/imu_preintegration.hpp"
//#include "optimizer/vertex_params.hpp"
//#include "tools/geometry.hpp"
//#include "tools/message_print.hpp"

namespace magnav {
class TollesLawson {
public:
  enum FDMscheme {
    BACKWARD,  // 1st derivative 1st-order backward difference
    FORWARD,   //  1st derivative 1st-order forward  difference
    CENTRAL,   //  1st derivative 2nd-order central  difference
    BACKWARD2, // 1st derivative 2nd-order backward difference
    FORWARD2,  // 1st derivative 2nd-order forward  difference
    FOURTH     //   4th derivative central difference
  };

  enum TLterm { PERMANENT, INDUCED, EDDY, FDM, BIAS };

  TollesLawson() {}

  /*
   * createMatrixA()
   * Create Tolles-Lawson `A` matrix using vector magnetometer measurements.
   * Optionally returns the magnitude and derivatives of total field.
   */
  bool createMatrixA(const std::vector<double> &Bx,
                     const std::vector<double> &By,
                     const std::vector<double> &Bz, std::vector<double> &Bt,
                     const std::unordered_set<TLterm> &terms = {PERMANENT,
                                                                INDUCED, EDDY},
                     const double Bt_scale = 50000.0);
  /*
   * createCoeff()
   * Create Tolles-Lawson coefficients using vector and scalar magnetometer
   * measurements and a bandpass, low-pass or high-pass filter.
   */
  double createCoeff(
      const std::vector<double> &Bx,
      const std::vector<double> &By, const std::vector<double> &Bz,
      const std::vector<double> &B, std::vector<double> &Bt,
      const std::unordered_set<TLterm> &terms = {PERMANENT, INDUCED, EDDY},
      const double lambda = 0.0, const double pass1 = 0.1,
      const double pass2 = 0.9, const double fs = 10.0, const int pole = 4,
      const int trim = 20, const double Bt_scale = 50000.0);

private:
  /*
   * fdm()
   * Finite difference method (FDM) on vector of input data.
   */
  bool fdm(std::vector<double> &dif, const std::vector<double> &x,
           const TollesLawson::FDMscheme scheme = CENTRAL);

bool bwbp_filter(std::vector<double> &y, const std::vector<double> &x,
    const double pass1, const double pass2, const int pole,
    const int trim);

void linear_regression(std::vector<double> &coeff, std::vector<double> &residual, const std::vector<double> &bb, const std::vector<std::vector<double>> &AA, const double lambda);

  //void getFilterCoeffAB(int n, double lowcut, double highcut, int fs,
  //                      std::vector<double> &acof_vec,
  //                      std::vector<double> &bcof_vec);

private:
  // A -  Tolles-Lawson `A` matrix
  std::vector<std::vector<double>> TL_A_;
  std::vector<std::vector<double>> TL_A_filt_;

	std::vector<double> TL_beta_;

  // [Bx, By, Bz] - vector magnetometer measurements [nT]
  // std::vector<double> &Bx;
  // std::vector<double> &By;
  // std::vector<double> &Bz;

  // B - scalar magnetometer measurements [nT]
  // std::vector<double> &B;

  // Bt - magnitude of vector magnetometer measurements or scalar magnetometer
  // measurements for modified Tolles-Lawson [nT]
  //	std::vector<double> &Bt;

  // B_dot - finite differences of total field vector [nT]
  std::vector<double> Bx_dot;
  std::vector<double> By_dot;
  std::vector<double> Bz_dot;
};
} // namespace magnav
#endif
