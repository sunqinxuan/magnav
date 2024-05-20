/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-18 11:19
#
# Filename:		mag_compensation.hpp
#
# Description:
#
************************************************/

#ifndef MAGNAV_MAGCMP_HPP_
#define MAGNAV_MAGCMP_HPP_

#include <deque>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <string>

#include "mag_compensation/tolles_lawson.hpp"

namespace magnav {

class H5Data {
public:
  int N;
  double dt;

  std::vector<double> line;
  std::vector<double> flight;
  std::vector<double> tt;
  std::vector<double> utm_x;
  std::vector<double> utm_y;
  std::vector<double> utm_z;
  std::vector<double> msl;
  std::vector<double> lat;
  std::vector<double> lon;

  std::vector<double> baro;
  std::vector<double> radar;
  std::vector<double> topo;
  std::vector<double> dem;
  std::vector<double> drape;

  std::vector<double> ins_pitch;
  std::vector<double> ins_roll;
  std::vector<double> ins_yaw;
  std::vector<double> diurnal;

  std::vector<double> mag_1_c;
  std::vector<double> mag_1_lag;
  std::vector<double> mag_1_dc;
  std::vector<double> mag_1_igrf;
  std::vector<double> mag_1_uc;

  std::vector<double> mag_2_uc;
  std::vector<double> mag_3_uc;
  std::vector<double> mag_4_uc;
  std::vector<double> mag_5_uc;

  std::vector<double> flux_b_x;
  std::vector<double> flux_b_y;
  std::vector<double> flux_b_z;
  std::vector<double> flux_b_t;

  std::vector<double> flux_c_x;
  std::vector<double> flux_c_y;
  std::vector<double> flux_c_z;
  std::vector<double> flux_c_t;

  std::vector<double> flux_d_x;
  std::vector<double> flux_d_y;
  std::vector<double> flux_d_z;
  std::vector<double> flux_d_t;

  std::vector<double> ogs_mag;
  std::vector<double> ogs_alt;

  std::vector<double> ins_acc_x;
  std::vector<double> ins_acc_y;
  std::vector<double> ins_acc_z;
  std::vector<double> ins_wander;
  std::vector<double> ins_lat;
  std::vector<double> ins_lon;
  std::vector<double> ins_alt;
  std::vector<double> ins_vn;
  std::vector<double> ins_vw;
  std::vector<double> ins_vu;

  std::vector<double> pitch_rate;
  std::vector<double> roll_rate;
  std::vector<double> yaw_rate;
  std::vector<double> lgtl_acc;
  std::vector<double> ltrl_acc;
  std::vector<double> nrml_acc;

  std::vector<double> tas; // TODO: true_as?
  std::vector<double> pitot_p;
  std::vector<double> static_p;
  std::vector<double> total_p;

  std::vector<double> cur_com_1;
  std::vector<double> cur_ac_hi;
  std::vector<double> cur_ac_lo;
  std::vector<double> cur_tank;
  std::vector<double> cur_flap;
  std::vector<double> cur_strb;
  std::vector<double> cur_srvo_o;
  std::vector<double> cur_srvo_m;
  std::vector<double> cur_srvo_i;
  std::vector<double> cur_heat;
  std::vector<double> cur_acpwr;
  std::vector<double> cur_outpwr;
  std::vector<double> cur_bat_1;
  std::vector<double> cur_bat_2;

  std::vector<double> vol_acpwr;
  std::vector<double> vol_outpwr;
  std::vector<double> vol_bat_1;
  std::vector<double> vol_bat_2;
  std::vector<double> vol_res_p;
  std::vector<double> vol_res_n;
  std::vector<double> vol_back_p;
  std::vector<double> vol_back_n;
  std::vector<double> vol_gyro_1;
  std::vector<double> vol_gyro_2;
  std::vector<double> vol_acc_p;
  std::vector<double> vol_acc_n;
  std::vector<double> vol_block;
  std::vector<double> vol_back;
  std::vector<double> vol_srvo;
  std::vector<double> vol_cabt;
  std::vector<double> vol_fan;
};

class MagCompensation {
public:
  using Vector = std::vector<double>;
  using Matrix = std::vector<std::vector<double>>;

  MagCompensation();

  void compensate(H5Data &out_data, const H5Data &calib_data,
                  const H5Data &xyz_data);
  void compensateVector(H5Data &out_data, const H5Data &calib_data,
                        const H5Data &xyz_data);
  void compensateVector_Component(H5Data &out_data, const H5Data &calib_data,
                                  const H5Data &xyz_data);

private:
  void product(Vector &result, const Matrix &A, const Vector &b);
  void demean(Vector &result, const Vector &v);
  void substract(Vector &result, const Vector &a, const Vector &b);
  void split_vector(const Vector &pile, Vector &x, Vector &y, Vector &z);

private:
  std::shared_ptr<TollesLawson> tl_model_;
};
} // namespace magnav

#endif
