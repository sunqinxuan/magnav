/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-18 11:19
#
# Filename:		mag_compensation.cpp
#
# Description:
#
************************************************/

#include "mag_compensation/mag_compensation.hpp"

namespace magnav {

MagCompensation::MagCompensation() {
  tl_model_ = std::make_shared<TollesLawson>();
}

void MagCompensation::compensateVector_Component(H5Data &out_data,
                                                 const H5Data &calib_data,
                                                 const H5Data &xyz_data) {
  int i1 = -1, i2 = -1;
  for (size_t i = 0; i < calib_data.line.size(); i++) {
    if (calib_data.line[i] == 1002.02) {
      i1 = i;
      break;
    }
  }
  for (size_t i = calib_data.line.size() - 1; i > 0; i--) {
    if (calib_data.line[i] == 1002.02) {
      i2 = i;
      break;
    }
  }
  if (i1 < 0 || i2 < 0) {
    ERROR("the value 1002.02 not found in calib_data.line");
    return;
  }
  std::cout << "i1 = " << i1 << std::endl;
  std::cout << "i2 = " << i2 << std::endl;

  // std::vector<double> Bx, By, Bz, B3, B4, B5, Be;
  std::vector<double> Bbx, Bby, Bbz, Bcx, Bcy, Bcz, Bdx, Bdy, Bdz;
  std::vector<double> Bb, Bc, Bd;
  for (int i = i1; i <= i2; i++) {
    Bbx.push_back(calib_data.flux_b_x[i]);
    Bby.push_back(calib_data.flux_b_y[i]);
    Bbz.push_back(calib_data.flux_b_z[i]);
    Bcx.push_back(calib_data.flux_c_x[i]);
    Bcy.push_back(calib_data.flux_c_y[i]);
    Bcz.push_back(calib_data.flux_c_z[i]);
    Bdx.push_back(calib_data.flux_d_x[i]);
    Bdy.push_back(calib_data.flux_d_y[i]);
    Bdz.push_back(calib_data.flux_d_z[i]);
    Bb.push_back(calib_data.flux_b_t[i]);
    Bc.push_back(calib_data.flux_c_t[i]);
    Bd.push_back(calib_data.flux_d_t[i]);
  }

  std::vector<double> TL_coef_c, TL_coef_d;
  tl_model_->createCoeff_Vector(TL_coef_c, Bcx, Bcy, Bcz, Bbx, Bby, Bbz);
  tl_model_->createCoeff_Vector(TL_coef_d, Bdx, Bdy, Bdz, Bbx, Bby, Bbz);
  // TODO: createCoeff_Vector

  std::vector<double> BBcx, BBcy, BBcz, BBc;
  std::vector<double> BBdx, BBdy, BBdz, BBd;
  for (int i = 0; i < xyz_data.flux_c_x.size(); i++) {
    BBcx.push_back(xyz_data.flux_c_x[i]);
    BBcy.push_back(xyz_data.flux_c_y[i]);
    BBcz.push_back(xyz_data.flux_c_z[i]);
    BBc.push_back(xyz_data.flux_c_t[i]);
  }
  for (int i = 0; i < xyz_data.flux_d_x.size(); i++) {
    BBdx.push_back(xyz_data.flux_d_x[i]);
    BBdy.push_back(xyz_data.flux_d_y[i]);
    BBdz.push_back(xyz_data.flux_d_z[i]);
    BBd.push_back(xyz_data.flux_d_t[i]);
  }

  // compensation of the vector Ba;
  std::vector<std::vector<double>> TL_A_c_vec, TL_A_d_vec;
  tl_model_->createMatrixA_Vector(TL_A_c_vec, BBcx, BBcy, BBcz);
  tl_model_->createMatrixA_Vector(TL_A_d_vec, BBdx, BBdy, BBdz);
  // std::ofstream fp("debug.txt", std::ios::out);
  // for (int i = 0; i < TL_A_c_vec.size(); i++) {
  //  for (int j = 0; j < TL_A_c_vec[i].size(); j++) {
  //    fp << TL_A_c_vec[i][j] << " ";
  //  }
  //  fp << std::endl;
  //}
  // fp.close();

  std::vector<double> comp_part_c_vec, comp_part_d_vec;
  product(comp_part_c_vec, TL_A_c_vec, TL_coef_c);
  product(comp_part_d_vec, TL_A_d_vec, TL_coef_d);

  std::vector<double> comp_part_c_x, comp_part_c_y, comp_part_c_z;
  std::vector<double> comp_part_d_x, comp_part_d_y, comp_part_d_z;
  split_vector(comp_part_c_vec, comp_part_c_x, comp_part_c_y, comp_part_c_z);
  split_vector(comp_part_d_vec, comp_part_d_x, comp_part_d_y, comp_part_d_z);

  std::vector<double> flux_c_x_comp, flux_c_y_comp, flux_c_z_comp;
  std::vector<double> flux_d_x_comp, flux_d_y_comp, flux_d_z_comp;
  substract(flux_c_x_comp, BBcx, comp_part_c_x);
  substract(flux_c_y_comp, BBcy, comp_part_c_y);
  substract(flux_c_z_comp, BBcz, comp_part_c_z);
  substract(flux_d_x_comp, BBdx, comp_part_d_x);
  substract(flux_d_y_comp, BBdy, comp_part_d_y);
  substract(flux_d_z_comp, BBdz, comp_part_d_z);

  out_data.flux_c_x = flux_c_x_comp;
  out_data.flux_c_y = flux_c_y_comp;
  out_data.flux_c_z = flux_c_z_comp;
  out_data.flux_d_x = flux_d_x_comp;
  out_data.flux_d_y = flux_d_y_comp;
  out_data.flux_d_z = flux_d_z_comp;

  std::vector<double> flux_c_t_comp, flux_d_t_comp;
  for (int i = 0; i < flux_c_x_comp.size(); i++) {
    flux_c_t_comp.push_back(sqrt(flux_c_x_comp[i] * flux_c_x_comp[i] +
                                 flux_c_y_comp[i] * flux_c_y_comp[i] +
                                 flux_c_z_comp[i] * flux_c_z_comp[i]));
  }
  for (int i = 0; i < flux_d_x_comp.size(); i++) {
    flux_d_t_comp.push_back(sqrt(flux_d_x_comp[i] * flux_d_x_comp[i] +
                                 flux_d_y_comp[i] * flux_d_y_comp[i] +
                                 flux_d_z_comp[i] * flux_d_z_comp[i]));
  }

  out_data.flux_c_t = flux_c_t_comp;
  out_data.flux_d_t = flux_d_t_comp;
}

void MagCompensation::compensateVector(H5Data &out_data,
                                       const H5Data &calib_data,
                                       const H5Data &xyz_data) {
  int i1 = -1, i2 = -1;
  for (size_t i = 0; i < calib_data.line.size(); i++) {
    if (calib_data.line[i] == 1002.02) {
      i1 = i;
      break;
    }
  }
  for (size_t i = calib_data.line.size() - 1; i > 0; i--) {
    if (calib_data.line[i] == 1002.02) {
      i2 = i;
      break;
    }
  }
  if (i1 < 0 || i2 < 0) {
    ERROR("the value 1002.02 not found in calib_data.line");
    return;
  }
  std::cout << "i1 = " << i1 << std::endl;
  std::cout << "i2 = " << i2 << std::endl;

  // std::vector<double> Bx, By, Bz, B3, B4, B5, Be;
  std::vector<double> Bbx, Bby, Bbz, Bcx, Bcy, Bcz, Bdx, Bdy, Bdz;
  std::vector<double> Bb, Bc, Bd;
  for (int i = i1; i <= i2; i++) {
    // Bx.push_back(calib_data.flux_b_x[i]);
    // By.push_back(calib_data.flux_b_y[i]);
    // Bz.push_back(calib_data.flux_b_z[i]);
    // B3.push_back(calib_data.mag_3_uc[i]);
    // B4.push_back(calib_data.mag_4_uc[i]);
    // B5.push_back(calib_data.mag_5_uc[i]);
    // Be.push_back(calib_data.mag_1_uc[i]);
    Bbx.push_back(calib_data.flux_b_x[i]);
    Bby.push_back(calib_data.flux_b_y[i]);
    Bbz.push_back(calib_data.flux_b_z[i]);
    Bcx.push_back(calib_data.flux_c_x[i]);
    Bcy.push_back(calib_data.flux_c_y[i]);
    Bcz.push_back(calib_data.flux_c_z[i]);
    Bdx.push_back(calib_data.flux_d_x[i]);
    Bdy.push_back(calib_data.flux_d_y[i]);
    Bdz.push_back(calib_data.flux_d_z[i]);
    //		Bb.push_back(sqrt(Bbx[i] * Bbx[i] + Bby[i] * Bby[i] + Bbz[i] *
    // Bbz[i])); 		Bc.push_back(sqrt(Bcx[i] * Bcx[i] + Bcy[i] *
    // Bcy[i]
    // + Bcz[i]
    // * Bcz[i])); 		Bd.push_back(sqrt(Bdx[i] * Bdx[i] + Bdy[i] *
    // Bdy[i] + Bdz[i] * Bdz[i]));
    Bb.push_back(calib_data.flux_b_t[i]);
    Bc.push_back(calib_data.flux_c_t[i]);
    Bd.push_back(calib_data.flux_d_t[i]);
  }

  std::vector<double> Bt; // empty Bt;

  // std::vector<double> TL_coef_3, TL_coef_4, TL_coef_5;
  // tl_model_->createCoeff(TL_coef_3, Bx, By, Bz, B3, Be, Bt);
  // tl_model_->createCoeff(TL_coef_4, Bx, By, Bz, B4, Be, Bt);
  // tl_model_->createCoeff(TL_coef_5, Bx, By, Bz, B5, Be, Bt);
  std::vector<double> TL_coef_c, TL_coef_d;
  tl_model_->createCoeff(TL_coef_c, Bcx, Bcy, Bcz, Bc, Bb, Bt);
  tl_model_->createCoeff(TL_coef_d, Bdx, Bdy, Bdz, Bd, Bb, Bt);

  // std::vector<double> BBx, BBy, BBz;
  // for (int i = 0; i < xyz_data.flux_b_x.size(); i++) {
  //  BBx.push_back(xyz_data.flux_b_x[i]);
  //  BBy.push_back(xyz_data.flux_b_y[i]);
  //  BBz.push_back(xyz_data.flux_b_z[i]);
  //}
  std::vector<double> BBcx, BBcy, BBcz, BBc;
  std::vector<double> BBdx, BBdy, BBdz, BBd;
  for (int i = 0; i < xyz_data.flux_c_x.size(); i++) {
    BBcx.push_back(xyz_data.flux_c_x[i]);
    BBcy.push_back(xyz_data.flux_c_y[i]);
    BBcz.push_back(xyz_data.flux_c_z[i]);
    //		BBc.push_back(sqrt(BBcx[i] * BBcx[i] + BBcy[i] * BBcy[i] +
    // BBcz[i]
    //* BBcz[i]));
    BBc.push_back(xyz_data.flux_c_t[i]);
  }
  for (int i = 0; i < xyz_data.flux_d_x.size(); i++) {
    BBdx.push_back(xyz_data.flux_d_x[i]);
    BBdy.push_back(xyz_data.flux_d_y[i]);
    BBdz.push_back(xyz_data.flux_d_z[i]);
    //		BBd.push_back(sqrt(BBdx[i] * BBdx[i] + BBdy[i] * BBdy[i] +
    // BBdz[i]
    //* BBdz[i]));
    BBd.push_back(xyz_data.flux_d_t[i]);
  }

  // std::vector<std::vector<double>> TL_A;
  // tl_model_->createMatrixA(TL_A, BBx, BBy, BBz, Bt);
  std::vector<std::vector<double>> TL_A_c, TL_A_d;
  tl_model_->createMatrixA(TL_A_c, BBcx, BBcy, BBcz, Bt);
  tl_model_->createMatrixA(TL_A_d, BBdx, BBdy, BBdz, Bt);

  // std::vector<double> comp_part_3, comp_part_4, comp_part_5;
  // product(comp_part_3, TL_A, TL_coef_3);
  // product(comp_part_4, TL_A, TL_coef_4);
  // product(comp_part_5, TL_A, TL_coef_5);
  std::vector<double> comp_part_c, comp_part_d;
  product(comp_part_c, TL_A_c, TL_coef_c);
  product(comp_part_d, TL_A_d, TL_coef_d);

  // std::vector<double> dm_comp_part_3, dm_comp_part_4, dm_comp_part_5;
  // demean(dm_comp_part_3, comp_part_3);
  // demean(dm_comp_part_4, comp_part_4);
  // demean(dm_comp_part_5, comp_part_5);
  //
  // std::vector<double> mag_3_c, mag_4_c, mag_5_c;
  // substract(mag_3_c, xyz_data.mag_3_uc, dm_comp_part_3);
  // substract(mag_4_c, xyz_data.mag_4_uc, dm_comp_part_4);
  // substract(mag_5_c, xyz_data.mag_5_uc, dm_comp_part_5);
  //
  // std::vector<double> correction;
  // substract(correction, xyz_data.mag_1_dc, xyz_data.mag_1_igrf);
  //
  // substract(mag_3_c, mag_3_c, xyz_data.diurnal);
  // substract(mag_4_c, mag_4_c, xyz_data.diurnal);
  // substract(mag_5_c, mag_5_c, xyz_data.diurnal);
  //
  // substract(mag_3_c, mag_3_c, correction);
  // substract(mag_4_c, mag_4_c, correction);
  // substract(mag_5_c, mag_5_c, correction);

  // std::vector<double> mag_3_c, mag_4_c, mag_5_c;
  // substract(mag_3_c, xyz_data.mag_3_uc, comp_part_3);
  // substract(mag_4_c, xyz_data.mag_4_uc, comp_part_4);
  // substract(mag_5_c, xyz_data.mag_5_uc, comp_part_5);
  std::vector<double> flux_c_t_comp, flux_d_t_comp;
  substract(flux_c_t_comp, BBc, comp_part_c);
  substract(flux_d_t_comp, BBd, comp_part_d);

  out_data.tt = xyz_data.tt;
  // out_data.mag_1_igrf = xyz_data.mag_1_igrf;
  // out_data.mag_1_uc = xyz_data.mag_1_uc;
  out_data.flux_c_t = flux_c_t_comp;
  out_data.flux_d_t = flux_d_t_comp;

  // ******************************
  // compensation of the vector Ba;
  std::vector<std::vector<double>> TL_A_c_vec, TL_A_d_vec;
  tl_model_->createMatrixA_Vector(TL_A_c_vec, BBcx, BBcy, BBcz);
  tl_model_->createMatrixA_Vector(TL_A_d_vec, BBdx, BBdy, BBdz);
  // std::ofstream fp("debug.txt", std::ios::out);
  // for (int i = 0; i < TL_A_c_vec.size(); i++) {
  //  for (int j = 0; j < TL_A_c_vec[i].size(); j++) {
  //    fp << TL_A_c_vec[i][j] << " ";
  //  }
  //  fp << std::endl;
  //}
  // fp.close();

  std::vector<double> comp_part_c_vec, comp_part_d_vec;
  product(comp_part_c_vec, TL_A_c_vec, TL_coef_c);
  product(comp_part_d_vec, TL_A_d_vec, TL_coef_d);

  std::vector<double> comp_part_c_x, comp_part_c_y, comp_part_c_z;
  std::vector<double> comp_part_d_x, comp_part_d_y, comp_part_d_z;
  split_vector(comp_part_c_vec, comp_part_c_x, comp_part_c_y, comp_part_c_z);
  split_vector(comp_part_d_vec, comp_part_d_x, comp_part_d_y, comp_part_d_z);

  std::vector<double> flux_c_x_comp, flux_c_y_comp, flux_c_z_comp;
  std::vector<double> flux_d_x_comp, flux_d_y_comp, flux_d_z_comp;
  substract(flux_c_x_comp, BBcx, comp_part_c_x);
  substract(flux_c_y_comp, BBcy, comp_part_c_y);
  substract(flux_c_z_comp, BBcz, comp_part_c_z);
  substract(flux_d_x_comp, BBdx, comp_part_d_x);
  substract(flux_d_y_comp, BBdy, comp_part_d_y);
  substract(flux_d_z_comp, BBdz, comp_part_d_z);

  out_data.flux_c_x = flux_c_x_comp;
  out_data.flux_c_y = flux_c_y_comp;
  out_data.flux_c_z = flux_c_z_comp;
  out_data.flux_d_x = flux_d_x_comp;
  out_data.flux_d_y = flux_d_y_comp;
  out_data.flux_d_z = flux_d_z_comp;
}

void MagCompensation::compensate(H5Data &out_data, const H5Data &calib_data,
                                 const H5Data &xyz_data) {
  int i1 = -1, i2 = -1;
  for (size_t i = 0; i < calib_data.line.size(); i++) {
    if (calib_data.line[i] == 1002.02) {
      i1 = i;
      break;
    }
  }
  for (size_t i = calib_data.line.size() - 1; i > 0; i--) {
    if (calib_data.line[i] == 1002.02) {
      i2 = i;
      break;
    }
  }
  if (i1 < 0 || i2 < 0) {
    ERROR("the value 1002.02 not found in calib_data.line");
    return;
  }
  std::cout << "i1 = " << i1 << std::endl;
  std::cout << "i2 = " << i2 << std::endl;

  std::vector<double> Bx, By, Bz, B3, B4, B5, Be;
  for (int i = i1; i <= i2; i++) {
    Bx.push_back(calib_data.flux_c_x[i]);
    By.push_back(calib_data.flux_c_y[i]);
    Bz.push_back(calib_data.flux_c_z[i]);
    B3.push_back(calib_data.mag_3_uc[i]);
    B4.push_back(calib_data.mag_4_uc[i]);
    B5.push_back(calib_data.mag_5_uc[i]);
    Be.push_back(calib_data.mag_1_uc[i]);
  }

  std::vector<double> Bt; // empty Bt;

  std::vector<double> TL_coef_3, TL_coef_4, TL_coef_5;
  tl_model_->createCoeff(TL_coef_3, Bx, By, Bz, B3, Be, Bt);
  tl_model_->createCoeff(TL_coef_4, Bx, By, Bz, B4, Be, Bt);
  tl_model_->createCoeff(TL_coef_5, Bx, By, Bz, B5, Be, Bt);

  std::vector<double> BBx, BBy, BBz;
  for (int i = 0; i < xyz_data.flux_b_x.size(); i++) {
    BBx.push_back(xyz_data.flux_c_x[i]);
    BBy.push_back(xyz_data.flux_c_y[i]);
    BBz.push_back(xyz_data.flux_c_z[i]);
  }

  std::vector<std::vector<double>> TL_A;
  tl_model_->createMatrixA(TL_A, BBx, BBy, BBz, Bt);

  std::vector<double> comp_part_3, comp_part_4, comp_part_5;
  product(comp_part_3, TL_A, TL_coef_3);
  product(comp_part_4, TL_A, TL_coef_4);
  product(comp_part_5, TL_A, TL_coef_5);

  // std::vector<double> dm_comp_part_3, dm_comp_part_4, dm_comp_part_5;
  // demean(dm_comp_part_3, comp_part_3);
  // demean(dm_comp_part_4, comp_part_4);
  // demean(dm_comp_part_5, comp_part_5);
  //
  // std::vector<double> mag_3_c, mag_4_c, mag_5_c;
  // substract(mag_3_c, xyz_data.mag_3_uc, dm_comp_part_3);
  // substract(mag_4_c, xyz_data.mag_4_uc, dm_comp_part_4);
  // substract(mag_5_c, xyz_data.mag_5_uc, dm_comp_part_5);
  //
  // std::vector<double> correction;
  // substract(correction, xyz_data.mag_1_dc, xyz_data.mag_1_igrf);
  //
  // substract(mag_3_c, mag_3_c, xyz_data.diurnal);
  // substract(mag_4_c, mag_4_c, xyz_data.diurnal);
  // substract(mag_5_c, mag_5_c, xyz_data.diurnal);
  //
  // substract(mag_3_c, mag_3_c, correction);
  // substract(mag_4_c, mag_4_c, correction);
  // substract(mag_5_c, mag_5_c, correction);

  std::vector<double> mag_3_c, mag_4_c, mag_5_c;
  substract(mag_3_c, xyz_data.mag_3_uc, comp_part_3);
  substract(mag_4_c, xyz_data.mag_4_uc, comp_part_4);
  substract(mag_5_c, xyz_data.mag_5_uc, comp_part_5);

  out_data.tt = xyz_data.tt;
  // out_data.mag_1_igrf = xyz_data.mag_1_igrf;
  // out_data.mag_1_uc = xyz_data.mag_1_uc;
  out_data.mag_3_uc = mag_3_c;
  out_data.mag_4_uc = mag_4_c;
  out_data.mag_5_uc = mag_5_c;

  // ******************************
  // compensation of the vector Ba;

  /*
std::vector<double> BBcx, BBcy, BBcz, BBc;
std::vector<double> BBdx, BBdy, BBdz, BBd;
for (int i = 0; i < xyz_data.flux_c_x.size(); i++) {
BBcx.push_back(xyz_data.flux_c_x[i]);
BBcy.push_back(xyz_data.flux_c_y[i]);
BBcz.push_back(xyz_data.flux_c_z[i]);
//		BBc.push_back(sqrt(BBcx[i] * BBcx[i] + BBcy[i] * BBcy[i] +
// BBcz[i]
//* BBcz[i]));
BBc.push_back(xyz_data.flux_c_t[i]);
}
for (int i = 0; i < xyz_data.flux_d_x.size(); i++) {
BBdx.push_back(xyz_data.flux_d_x[i]);
BBdy.push_back(xyz_data.flux_d_y[i]);
BBdz.push_back(xyz_data.flux_d_z[i]);
//		BBd.push_back(sqrt(BBdx[i] * BBdx[i] + BBdy[i] * BBdy[i] +
// BBdz[i]
//* BBdz[i]));
BBd.push_back(xyz_data.flux_d_t[i]);
}

std::vector<std::vector<double>> TL_A_c, TL_A_d;
tl_model_->createMatrixA(TL_A_c, BBcx, BBcy, BBcz, Bt);
tl_model_->createMatrixA(TL_A_d, BBdx, BBdy, BBdz, Bt);

std::vector<double> comp_part_c, comp_part_d;
product(comp_part_c, TL_A_c, TL_coef_3);
product(comp_part_d, TL_A_d, TL_coef_3);

std::vector<double> flux_c_t_comp, flux_d_t_comp;
substract(flux_c_t_comp, BBc, comp_part_c);
substract(flux_d_t_comp, BBd, comp_part_d);

out_data.tt = xyz_data.tt;
// out_data.mag_1_igrf = xyz_data.mag_1_igrf;
// out_data.mag_1_uc = xyz_data.mag_1_uc;
out_data.flux_c_t = flux_c_t_comp;
out_data.flux_d_t = flux_d_t_comp;

std::vector<std::vector<double>> TL_A_c_vec, TL_A_d_vec;
tl_model_->createMatrixA_Vector(TL_A_c_vec, BBcx, BBcy, BBcz);
tl_model_->createMatrixA_Vector(TL_A_d_vec, BBdx, BBdy, BBdz);
// std::ofstream fp("debug.txt", std::ios::out);
// for (int i = 0; i < TL_A_c_vec.size(); i++) {
//  for (int j = 0; j < TL_A_c_vec[i].size(); j++) {
//    fp << TL_A_c_vec[i][j] << " ";
//  }
//  fp << std::endl;
//}
// fp.close();

std::vector<double> comp_part_c_vec, comp_part_d_vec;
product(comp_part_c_vec, TL_A_c_vec, TL_coef_3);
product(comp_part_d_vec, TL_A_d_vec, TL_coef_3);

std::vector<double> comp_part_c_x, comp_part_c_y, comp_part_c_z;
std::vector<double> comp_part_d_x, comp_part_d_y, comp_part_d_z;
split_vector(comp_part_c_vec, comp_part_c_x, comp_part_c_y, comp_part_c_z);
split_vector(comp_part_d_vec, comp_part_d_x, comp_part_d_y, comp_part_d_z);

std::vector<double> flux_c_x_comp, flux_c_y_comp, flux_c_z_comp;
std::vector<double> flux_d_x_comp, flux_d_y_comp, flux_d_z_comp;
substract(flux_c_x_comp, BBcx, comp_part_c_x);
substract(flux_c_y_comp, BBcy, comp_part_c_y);
substract(flux_c_z_comp, BBcz, comp_part_c_z);
substract(flux_d_x_comp, BBdx, comp_part_d_x);
substract(flux_d_y_comp, BBdy, comp_part_d_y);
substract(flux_d_z_comp, BBdz, comp_part_d_z);

out_data.flux_c_x = flux_c_x_comp;
out_data.flux_c_y = flux_c_y_comp;
out_data.flux_c_z = flux_c_z_comp;
out_data.flux_d_x = flux_d_x_comp;
out_data.flux_d_y = flux_d_y_comp;
out_data.flux_d_z = flux_d_z_comp;
  */
}

void MagCompensation::product(Vector &result, const Matrix &A,
                              const Vector &b) {
  if (A.size() != b.size()) {
    ERROR("[MagCompensation][product] dimensions do not match!");
    return;
  }
  int dim = A.size();
  int N = A[0].size();
  result.resize(N);
  for (int i = 0; i < N; i++) {
    result[i] = 0.0;
    for (int j = 0; j < dim; j++) {
      result[i] += A[j][i] * b[j];
    }
  }
}

void MagCompensation::demean(Vector &result, const Vector &v) {
  //	result=v;
  //	return;

  double mean = 0.0;
  for (size_t i = 0; i < v.size(); i++) {
    mean += v[i];
  }
  mean /= v.size();

  result.resize(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    result[i] = v[i] - mean;
  }
}

void MagCompensation::substract(Vector &result, const Vector &a,
                                const Vector &b) {
  if (a.size() != b.size()) {
    ERROR("[MagCompensation][substract] input vectors are not of the same "
          "dimension!");
    return;
  }
  result.resize(a.size());
  for (size_t i = 0; i < a.size(); i++) {
    result[i] = a[i] - b[i];
  }
}

void MagCompensation::split_vector(const Vector &pile, Vector &x, Vector &y,
                                   Vector &z) {
  if (pile.size() % 3 != 0 || pile.size() == 0) {
    ERROR("[MagCompensation][split_vector] wrong input!");
    return;
  }
  x.clear();
  y.clear();
  z.clear();
  for (size_t i = 0; i < pile.size(); i++) {
    if (i % 3 == 0) {
      x.push_back(pile[i]);
    } else if (i % 3 == 1) {
      y.push_back(pile[i]);
    } else if (i % 3 == 2) {
      z.push_back(pile[i]);
    }
  }
}
} // namespace magnav
