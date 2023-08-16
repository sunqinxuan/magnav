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
#include <vector>
#include <unordered_set>

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
	enum FDMscheme
	{
    BACKWARD,// 1st derivative 1st-order backward difference
    FORWARD,//  1st derivative 1st-order forward  difference
    CENTRAL,//  1st derivative 2nd-order central  difference
    BACKWARD2,//1st derivative 2nd-order backward difference
    FORWARD2,// 1st derivative 2nd-order forward  difference
    FOURTH//   4th derivative central difference
	};

	enum TLterm
	{
		PERMANENT,
		INDUCED,
		EDDY,
		FDM,
		BIAS
	};

  TollesLawson(){}

	/*
    create_TL_A(Bx, By, Bz;
                Bt       = sqrt.(Bx.^2+By.^2+Bz.^2),
                terms    = [:permanent,:induced,:eddy],
                Bt_scale = 50000,
                return_B = false)

Create Tolles-Lawson `A` matrix using vector magnetometer measurements. 
Optionally returns the magnitude and derivatives of total field.

**Arguments:**
- `Bx,By,Bz`: vector magnetometer measurements [nT]
- `Bt`:       (optional) magnitude of vector magnetometer measurements or scalar magnetometer measurements for modified Tolles-Lawson [nT]
- `terms`:    (optional) Tolles-Lawson terms to use {`:permanent`,`:induced`,`:eddy`,`:bias`}
- `Bt_scale`: (optional) scaling factor for induced and eddy current terms [nT]
- `return_B`: (optional) if true, also return `Bt` and `B_dot`

**Returns:**
- `A`:     Tolles-Lawson `A` matrix
- `Bt`:    if `return_B = true`, magnitude of total field measurements [nT]
- `B_dot`: if `return_B = true`, finite differences of total field vector [nT]
*/
bool createMatrixA(const std::vector<double> &Bx,const std::vector<double> &By,const std::vector<double> &Bz,
																 std::vector<double> &Bt,
																 const std::unordered_set<TLterm> &terms={PERMANENT,INDUCED,EDDY},
																 const double Bt_scale=50000.0,
																 const bool return_B=false);

	/*
    fdm(x::Vector; scheme::Symbol=:central)

Finite difference method (FDM) on vector of input data.

**Arguments:**
- `x`:      input data
- `scheme`: (optional) finite difference method scheme used
    - `backward`:  1st derivative 1st-order backward difference
    - `forward`:   1st derivative 1st-order forward  difference
    - `central`:   1st derivative 2nd-order central  difference
    - `backward2`: 1st derivative 2nd-order backward difference
    - `forward2`:  1st derivative 2nd-order forward  difference
    - `fourth`:    4th derivative central difference

**Returns:**
- `dif`: length of `x` finite differences
*/
bool fdm(std::vector<double> &dif, const std::vector<double> &x, const TollesLawson::FDMscheme scheme=CENTRAL);

private:
private:
};
} // namespace magnav
#endif
