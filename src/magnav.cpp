/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-16 14:14
#
# Filename:		magnav.cpp
#
# Description:
#
************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <ceres/ceres.h>

#include "H5Cpp.h"
#include "mag_compensation/mag_compensation.hpp"
// using namespace H5;
using namespace magnav;

// const H5std_string FILE_NAME("/home/sun/magnav/data/Flt1002_train.h5");
// const H5std_string DATASET_NAME("flux_b_x");
// const int NX_SUB = 3; // hyperslab dimensions
// const int NY_SUB = 4;
// const int NX = 7; // output buffer dimensions
// const int NY = 7;
// const int NZ = 3;
// const int RANK_OUT = 3;

void readH5Data(const int N, const H5::H5File &file, const std::string field,
                std::vector<double> &data) {
  // H5::H5File file(filename, H5F_ACC_RDONLY);
  H5::DataSet dataset = file.openDataSet(field);

  data.resize(N);
  dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
}

void writeH5Data(
    H5::H5File &file, const std::string field,
    const std::vector<double>
        &data) //, const H5::DataType &datatype, const H5::DataSpace &dataspace)
{
  hsize_t dimsf[1]; // dataset dimensions
  dimsf[0] = data.size();
  H5::DataSpace dataspace(1, dimsf);
  H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

  H5::DataSet dataset = file.createDataSet(field, datatype, dataspace);
  dataset.write(data.data(), datatype);
}

void getFlightData(const std::string filename, H5Data &h5data) {
  H5::H5File file(filename, H5F_ACC_RDONLY);
  std::cout << "opening file: " << filename << std::endl;
  H5::DataSet data_N = file.openDataSet("N");
  data_N.read(&h5data.N, H5::PredType::NATIVE_INT64);
  // std::cout << "N=" << h5data.N << std::endl;
  H5::DataSet data_dt = file.openDataSet("dt");
  data_dt.read(&h5data.dt, H5::PredType::NATIVE_DOUBLE);
  // std::cout << "dt=" << h5data.dt << std::endl;

  readH5Data(h5data.N, file, "line", h5data.line);
  readH5Data(h5data.N, file, "flight", h5data.flight);
  readH5Data(h5data.N, file, "tt", h5data.tt);
  readH5Data(h5data.N, file, "utm_x", h5data.utm_x);
  readH5Data(h5data.N, file, "utm_y", h5data.utm_y);
  readH5Data(h5data.N, file, "utm_z", h5data.utm_z);
  readH5Data(h5data.N, file, "msl", h5data.msl);
  readH5Data(h5data.N, file, "lat", h5data.lat);
  readH5Data(h5data.N, file, "lon", h5data.lon);

  readH5Data(h5data.N, file, "baro", h5data.baro);
  readH5Data(h5data.N, file, "radar", h5data.radar);
  readH5Data(h5data.N, file, "topo", h5data.topo);
  readH5Data(h5data.N, file, "dem", h5data.dem);
  readH5Data(h5data.N, file, "drape", h5data.drape);

  readH5Data(h5data.N, file, "ins_pitch", h5data.ins_pitch);
  readH5Data(h5data.N, file, "ins_roll", h5data.ins_roll);
  readH5Data(h5data.N, file, "ins_yaw", h5data.ins_yaw);
  readH5Data(h5data.N, file, "diurnal", h5data.diurnal);

  readH5Data(h5data.N, file, "mag_1_c", h5data.mag_1_c);
  readH5Data(h5data.N, file, "mag_1_lag", h5data.mag_1_lag);
  readH5Data(h5data.N, file, "mag_1_dc", h5data.mag_1_dc);
  readH5Data(h5data.N, file, "mag_1_igrf", h5data.mag_1_igrf);
  readH5Data(h5data.N, file, "mag_1_uc", h5data.mag_1_uc);

  readH5Data(h5data.N, file, "mag_2_uc", h5data.mag_2_uc);
  readH5Data(h5data.N, file, "mag_3_uc", h5data.mag_3_uc);
  readH5Data(h5data.N, file, "mag_4_uc", h5data.mag_4_uc);
  readH5Data(h5data.N, file, "mag_5_uc", h5data.mag_5_uc);

  readH5Data(h5data.N, file, "flux_b_x", h5data.flux_b_x);
  readH5Data(h5data.N, file, "flux_b_y", h5data.flux_b_y);
  readH5Data(h5data.N, file, "flux_b_z", h5data.flux_b_z);
  readH5Data(h5data.N, file, "flux_b_t", h5data.flux_b_t);

  readH5Data(h5data.N, file, "flux_c_x", h5data.flux_c_x);
  readH5Data(h5data.N, file, "flux_c_y", h5data.flux_c_y);
  readH5Data(h5data.N, file, "flux_c_z", h5data.flux_c_z);
  readH5Data(h5data.N, file, "flux_c_t", h5data.flux_c_t);

  readH5Data(h5data.N, file, "flux_d_x", h5data.flux_d_x);
  readH5Data(h5data.N, file, "flux_d_y", h5data.flux_d_y);
  readH5Data(h5data.N, file, "flux_d_z", h5data.flux_d_z);
  readH5Data(h5data.N, file, "flux_d_t", h5data.flux_d_t);

  readH5Data(h5data.N, file, "ogs_mag", h5data.ogs_mag);
  readH5Data(h5data.N, file, "ogs_alt", h5data.ogs_alt);

  readH5Data(h5data.N, file, "ins_acc_x", h5data.ins_acc_x);
  readH5Data(h5data.N, file, "ins_acc_y", h5data.ins_acc_y);
  readH5Data(h5data.N, file, "ins_acc_z", h5data.ins_acc_z);
  readH5Data(h5data.N, file, "ins_wander", h5data.ins_wander);
  readH5Data(h5data.N, file, "ins_lat", h5data.ins_lat);
  readH5Data(h5data.N, file, "ins_lon", h5data.ins_lon);
  readH5Data(h5data.N, file, "ins_alt", h5data.ins_alt);
  readH5Data(h5data.N, file, "ins_vn", h5data.ins_vn);
  readH5Data(h5data.N, file, "ins_vw", h5data.ins_vw);
  readH5Data(h5data.N, file, "ins_vu", h5data.ins_vu);

  readH5Data(h5data.N, file, "pitch_rate", h5data.pitch_rate);
  readH5Data(h5data.N, file, "roll_rate", h5data.roll_rate);
  readH5Data(h5data.N, file, "yaw_rate", h5data.yaw_rate);
  readH5Data(h5data.N, file, "lgtl_acc", h5data.lgtl_acc);
  readH5Data(h5data.N, file, "ltrl_acc", h5data.ltrl_acc);
  readH5Data(h5data.N, file, "nrml_acc", h5data.nrml_acc);

  readH5Data(h5data.N, file, "tas", h5data.tas);
  readH5Data(h5data.N, file, "pitot_p", h5data.pitot_p);
  readH5Data(h5data.N, file, "static_p", h5data.static_p);
  readH5Data(h5data.N, file, "total_p", h5data.total_p);

  readH5Data(h5data.N, file, "cur_com_1", h5data.cur_com_1);
  readH5Data(h5data.N, file, "cur_ac_hi", h5data.cur_ac_hi);
  readH5Data(h5data.N, file, "cur_ac_lo", h5data.cur_ac_lo);
  readH5Data(h5data.N, file, "cur_tank", h5data.cur_tank);
  readH5Data(h5data.N, file, "cur_flap", h5data.cur_flap);
  readH5Data(h5data.N, file, "cur_strb", h5data.cur_strb);
  readH5Data(h5data.N, file, "cur_srvo_o", h5data.cur_srvo_o);
  readH5Data(h5data.N, file, "cur_srvo_m", h5data.cur_srvo_m);
  readH5Data(h5data.N, file, "cur_srvo_i", h5data.cur_srvo_i);
  readH5Data(h5data.N, file, "cur_heat", h5data.cur_heat);
  readH5Data(h5data.N, file, "cur_acpwr", h5data.cur_acpwr);
  readH5Data(h5data.N, file, "cur_outpwr", h5data.cur_outpwr);
  readH5Data(h5data.N, file, "cur_bat_1", h5data.cur_bat_1);
  readH5Data(h5data.N, file, "cur_bat_2", h5data.cur_bat_2);

  readH5Data(h5data.N, file, "vol_acpwr", h5data.vol_acpwr);
  readH5Data(h5data.N, file, "vol_outpwr", h5data.vol_outpwr);
  readH5Data(h5data.N, file, "vol_bat_1", h5data.vol_bat_1);
  readH5Data(h5data.N, file, "vol_bat_2", h5data.vol_bat_2);
  readH5Data(h5data.N, file, "vol_res_p", h5data.vol_res_p);
  readH5Data(h5data.N, file, "vol_res_n", h5data.vol_res_n);
  readH5Data(h5data.N, file, "vol_back_p", h5data.vol_back_p);
  readH5Data(h5data.N, file, "vol_back_n", h5data.vol_back_n);
  readH5Data(h5data.N, file, "vol_gyro_1", h5data.vol_gyro_1);
  readH5Data(h5data.N, file, "vol_gyro_2", h5data.vol_gyro_2);
  readH5Data(h5data.N, file, "vol_acc_p", h5data.vol_acc_p);
  readH5Data(h5data.N, file, "vol_acc_n", h5data.vol_acc_n);
  readH5Data(h5data.N, file, "vol_block", h5data.vol_block);
  readH5Data(h5data.N, file, "vol_back", h5data.vol_back);
  readH5Data(h5data.N, file, "vol_srvo", h5data.vol_srvo);
  readH5Data(h5data.N, file, "vol_cabt", h5data.vol_cabt);
  readH5Data(h5data.N, file, "vol_fan", h5data.vol_fan);
}

int main1(void) {
  H5Data h5data1002, h5data1003;

  getFlightData("/home/sun/magnav/data/Flt1002_train.h5", h5data1002);
  getFlightData("/home/sun/magnav/data/Flt1003_train.h5", h5data1003);

  MagCompensation mag_comp;
  H5Data out_data;
  // mag_comp.compensate(out_data, h5data1002, h5data1003);
  // mag_comp.compensateVector(out_data, h5data1002, h5data1003);
  mag_comp.compensateVector_Component(out_data, h5data1002, h5data1003);

  H5::H5File file("data_TL.h5", H5F_ACC_TRUNC);
  // hsize_t dimsf[1]; // dataset dimensions
  // dimsf[0] = out_data.tt.size();
  //// dimsf[1] = NY;
  // H5::DataSpace dataspace(1, dimsf);
  //// IntType datatype(PredType::NATIVE_INT);
  // H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

  writeH5Data(file, "tt", out_data.tt);             //,datatype, dataspace);
  writeH5Data(file, "mag_1_uc", out_data.mag_1_uc); //,datatype, dataspace);
  writeH5Data(file, "mag_3_c", out_data.mag_3_uc);  //,datatype, dataspace);
  writeH5Data(file, "mag_4_c", out_data.mag_4_uc);  //,datatype, dataspace);
  writeH5Data(file, "mag_5_c", out_data.mag_5_uc);  //,datatype, dataspace);
  writeH5Data(file, "flux_c_t", out_data.flux_c_t); //,datatype, dataspace);
  writeH5Data(file, "flux_c_x", out_data.flux_c_x); //,datatype, dataspace);
  writeH5Data(file, "flux_c_y", out_data.flux_c_y); //,datatype, dataspace);
  writeH5Data(file, "flux_c_z", out_data.flux_c_z); //,datatype, dataspace);
  writeH5Data(file, "flux_d_t", out_data.flux_d_t); //,datatype, dataspace);
  writeH5Data(file, "flux_d_x", out_data.flux_d_x); //,datatype, dataspace);
  writeH5Data(file, "flux_d_y", out_data.flux_d_y); //,datatype, dataspace);
  writeH5Data(file, "flux_d_z", out_data.flux_d_z); //,datatype, dataspace);

  // H5::DataSet dataset = file.createDataSet("tt", datatype, dataspace);
  // dataset.write(out_data.tt.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // dataset = file.createDataSet("mag_3_c", datatype, dataspace);
  // dataset.write(out_data.mag_3_uc.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // dataset = file.createDataSet("mag_4_c", datatype, dataspace);
  // dataset.write(out_data.mag_4_uc.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // dataset = file.createDataSet("mag_5_c", datatype, dataspace);
  // dataset.write(out_data.mag_5_uc.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // dataset = file.createDataSet("flux_c_t", datatype, dataspace);
  // dataset.write(out_data.flux_c_t.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // dataset = file.createDataSet("flux_d_t", datatype, dataspace);
  // dataset.write(out_data.flux_d_t.data(), H5::PredType::NATIVE_DOUBLE);

  return 0;
}

bool readH5DataSet(const H5::H5File &file, const std::string field,
                   std::vector<double> &data, std::vector<int> &dims_out) {
  H5::DataSet dataset = file.openDataSet(field);
  /*
   * Get dataspace of the dataset.
   */
  H5::DataSpace dataspace = dataset.getSpace();

  /*
   * Get the number of dimensions in the dataspace.
   */
  int rank = dataspace.getSimpleExtentNdims();
  std::cout << "rank = " << rank << std::endl;

  /*
   * Get the dimension size of each dimension in the dataspace and
   * display them.
   */
  hsize_t dims[rank];
  dataspace.getSimpleExtentDims(dims, NULL);
  hsize_t data_size = 1;
  dims_out.clear();
  for (int i = 0; i < rank; i++) {
    std::cout << "dims[" << i << "] = " << dims[i] << std::endl;
    data_size *= dims[i];
    dims_out.push_back(dims[i]);
  }

  data.resize(data_size);
  dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
  std::cout << "data.size() = " << data.size() << std::endl;

  return true;
}

double getAnomalyValue(const std::vector<double> &map, const double x,
                       const double y) {
  return 0;
}

class CostFunctionCreator
{
public:
  CostFunctionCreator(const Eigen::Matrix3d &coeff_D_tilde_inv, const Eigen::Vector3d &coeff_o_hat, const Eigen::Matrix3d &orientation_r_c, const Eigen::Vector3d &mag_r,const Eigen::Vector3d  &mag_c)
      : coeff_D_tilde_inv_(coeff_D_tilde_inv), coeff_o_hat_(coeff_o_hat), orientation_r_c_(orientation_r_c), mag_r_(mag_r), mag_c_(mag_c)
  {
  }

  template <typename T>
  bool operator()(const T *const quaternion, T *residual_ptr) const
  {
    //Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_i(position_i);
    //Eigen::Map<const Eigen::Quaternion<T>> q_i(orientation_i);
//
    //Eigen::Map<const Eigen::Matrix<T, 3, 1>> p_j(position_j);
    //Eigen::Map<const Eigen::Quaternion<T>> q_j(orientation_j);
//
    //Eigen::Quaternion<T> q_i_inv = q_i.conjugate();
    //Eigen::Quaternion<T> q_ij = q_i_inv * q_j;
    //Eigen::Matrix<T, 3, 1> p_ij = q_i_inv * (p_j - p_i);
//
    //Eigen::Quaternion<T> q_ij_meas(pose_ij_meas_.linear().template cast<T>());
    //Eigen::Quaternion<T> delta_q = q_ij_meas * q_ij.conjugate();
//
    //Eigen::Map<Eigen::Matrix<T, 6, 1>> residual(residual_ptr);
    //Eigen::Map<Eigen::Matrix<T, 3, 1>> residual_trs(residual_ptr);
    //Eigen::Map<Eigen::Matrix<T, 3, 1>> residual_rot(residual_ptr + 3);
//
    //residual_trs = p_ij - pose_ij_meas_.translation().template cast<T>();
    //residual_rot = T(2.0) * delta_q.vec();
    //residual.applyOnTheLeft(sqrt_information_.template cast<T>());

    return true;
  }

  static ceres::CostFunction *Create(const Eigen::Matrix3d &coeff_D_tilde_inv, const Eigen::Vector3d &coeff_o_hat, const Eigen::Matrix3d &orientation_r_c, const Eigen::Vector3d &mag_r,const Eigen::Vector3d  &mag_c)
  {
    return new ceres::AutoDiffCostFunction<CostFunctionCreator, 3, 4>(new CostFunctionCreator(coeff_D_tilde_inv, coeff_o_hat, orientation_r_c, mag_r,mag_c));
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
	Eigen::Matrix3d coeff_D_tilde_inv_;
	Eigen::Vector3d coeff_o_hat_;
	Eigen::Matrix3d orientation_r_c_; // rotation from b_k to b_k+1;
	Eigen::Vector3d mag_r_, mag_c_;
};



/*
    euler2dcm(roll, pitch, yaw)

Converts a (Euler) roll-pitch-yaw (`X`-`Y`-`Z`) right-handed body to navigation
frame rotation (or the opposite rotation), to a DCM (direction cosine matrix).
Yaw is synonymous with azimuth and heading here.
If frame 1 is rotated to frame 2, then the returned DCM, when pre-multiplied,
rotates a vector in frame 1 into frame 2. There are 2 use cases:

1) With `order = :body2nav`, the body frame is rotated in the standard
-roll, -pitch, -yaw sequence to the navigation frame. For example, if v1 is
a 3x1 vector in the body frame [nose, right wing, down], then that vector
rotated into the navigation frame [north, east, down] would be v2 = dcm * v1.

2) With `order = :nav2body`, the navigation frame is rotated in the standard
yaw, pitch, roll sequence to the body frame. For example, if v1 is a 3x1
vector in the navigation frame [north, east, down], then that vector rotated
into the body frame [nose, right wing, down] would be v2 = dcm * v1.

Reference: Titterton & Weston, Strapdown Inertial Navigation Technology, 2004,
Section 3.6 (pg. 36-41 & 537).

**Arguments:**
- `roll`:  length roll  angle [rad], right-handed rotation about x-axis
- `pitch`: length pitch angle [rad], right-handed rotation about y-axis
- `yaw`:   length yaw   angle [rad], right-handed rotation about z-axis

**Returns:**
- `dcm`: `3` x `3` direction cosine matrix [-]
*/
Eigen::Matrix3d euler2dcm(double roll, double pitch, double yaw)
{
    double cr = cos(roll);
    double sr = sin(roll);
    double cp = cos(pitch);
    double sp = sin(pitch);
    double cy = cos(yaw);
    double sy = sin(yaw);

		Eigen::Matrix3d dcm=Eigen::Matrix3d::Zero();
		dcm(0, 0) =  cp * cy;
		dcm(0, 1) = -cr * sy + sr * sp * cy;
		dcm(0, 2) =  sr * sy + cr * sp * cy;
		dcm(1, 0) =  cp * sy;
		dcm(1, 1) =  cr * cy + sr * sp * sy;
		dcm(1, 2) = -sr * cy + cr * sp * sy;
		dcm(2, 0) = -sp;
		dcm(2, 1) =  sr * cp;
		dcm(2, 2) =  cr * cp;

		return dcm;
}


int main(void) {
  H5Data calib_data;

  getFlightData("/home/sun/magnav/data/Flt1002_train.h5", calib_data);

  std::vector<double> lines = {1002.02, 1002.20};
  std::vector<double> idx1_lines(lines.size(), -1);
  std::vector<double> idx2_lines(lines.size(), -1);
  for (size_t l = 0; l < lines.size(); l++) {
    for (size_t i = 0; i < calib_data.line.size(); i++) {
      if (calib_data.line[i] == lines[l]) {
        idx1_lines[l] = i;
        break;
      }
    }
    for (size_t i = calib_data.line.size() - 1; i > 0; i--) {
      if (calib_data.line[i] == lines[l]) {
        idx2_lines[l] = i;
        break;
      }
    }
    if (idx1_lines[l] < 0 || idx2_lines[l] < 0) {
      ERROR("the value ", lines[l], " not found in calib_data.line");
      return -1;
    }
    std::cout << "line " << lines[l] << " from " << idx1_lines[l] << " to "
              << idx2_lines[l] << std::endl;
  }

  std::vector<double> x_m, y_m, z_m; // mag_x, mag_y, mag_z;
	std::vector<double> ins_pitch, ins_roll, ins_yaw; 
  for (size_t l = 0; l < idx1_lines.size(); l++) {
    int i1 = idx1_lines[l];
    int i2 = idx2_lines[l];
    // load flux_c data
    // load ins_pitch, ins_roll, ins_yaw
    for (size_t i = i1; i <= i2; i++) {
      x_m.push_back(calib_data.flux_c_x[i]);
      y_m.push_back(calib_data.flux_c_y[i]);
      z_m.push_back(calib_data.flux_c_z[i]);
			ins_pitch.push_back(calib_data.ins_pitch[i]);
			ins_roll.push_back(calib_data.ins_roll[i]);
			ins_yaw.push_back(calib_data.ins_yaw[i]);
    }
  }

  // read anomaly map file;
  /*
H5::H5File anomaly_map_file("/home/sun/magnav/data/Canada_MAG_RES_200m.hdf5",
                        H5F_ACC_RDONLY);

std::vector<double> map;
std::vector<int> dims_map;
std::cout << std::endl
      << "reading Canada_MAG_RES_200m.hdf5 --- 'map' " << std::endl;
readH5DataSet(anomaly_map_file, "map", map, dims_map);

// std::ofstream fp("map.txt", std::ios::out);
// for (int i = 0; i < dims_map[0]; i++) {
//  for (int j = 0; j < dims_map[1]; j++) {
//    fp << map[i * dims_map[1] + j] << " ";
//  }
//  fp << std::endl;
//}
  //fp.close();

std::vector<double> xx;
std::vector<int> dims_xx;
std::cout << std::endl
      << "reading Canada_MAG_RES_200m.hdf5 --- 'xx' " << std::endl;
readH5DataSet(anomaly_map_file, "xx", xx, dims_xx);

std::vector<double> yy;
std::vector<int> dims_yy;
std::cout << std::endl
      << "reading Canada_MAG_RES_200m.hdf5 --- 'yy' " << std::endl;
readH5DataSet(anomaly_map_file, "yy", yy, dims_yy);
  */




  //*****************************************************
  // ellipsoid fitting;
  std::vector<double> x(x_m);
  std::vector<double> y(y_m);
  std::vector<double> z(z_m);

  const int N = x.size();
  Eigen::MatrixXd design_matrix(10, N);
  for (int i = 0; i < N; i++) {
    design_matrix(0, i) = x[i] * x[i];
    design_matrix(1, i) = y[i] * y[i];
    design_matrix(2, i) = z[i] * z[i];
    design_matrix(3, i) = 2 * y[i] * z[i];
    design_matrix(4, i) = 2 * x[i] * z[i];
    design_matrix(5, i) = 2 * x[i] * y[i];
    design_matrix(6, i) = 2 * x[i];
    design_matrix(7, i) = 2 * y[i];
    design_matrix(8, i) = 2 * z[i];
    design_matrix(9, i) = 1.0;
  }
  // std::ofstream fp("design_matrix.txt", std::ios::out);
  // for (int i = 0; i < design_matrix.rows(); i++) {
  //  for (int j = 0; j < design_matrix.cols(); j++) {
  //    fp << std::fixed << design_matrix(i, j) << " ";
  //  }
  //  fp << std::endl;
  //}
  // fp.close();

  // Constraint kJ > I^2
  // Ellipsoid if k = 4
  const int k = 4;

  // Eqn(7)
  Eigen::Matrix<double, 6, 6> C1;
  C1 << -1, 0.5 * k - 1, 0.5 * k - 1, 0, 0, 0, 0.5 * k - 1, -1, 0.5 * k - 1, 0,
      0, 0, 0.5 * k - 1, 0.5 * k - 1, -1, 0, 0, 0, 0, 0, 0, -k, 0, 0, 0, 0, 0,
      0, -k, 0, 0, 0, 0, 0, 0, -k;

  // Eqn(11)
  Eigen::MatrixXd S = design_matrix * design_matrix.transpose();
  Eigen::MatrixXd S11 = S.block(0, 0, 6, 6); // 6X6
  Eigen::MatrixXd S12 = S.block(0, 6, 6, 4); // 6X4
  Eigen::MatrixXd S21 = S.block(6, 0, 4, 6); // 4X6
  Eigen::MatrixXd S22 = S.block(6, 6, 4, 4); // 4X4

  // Output the matrices (for debugging purposes)
  // cout << "design matrix:\n" << design_matrix << endl;
  // cout << "S matrix:\n" << S << endl;
  // cout << "S11 matrix:\n" << S11 << endl;
  // cout << "S12 matrix:\n" << S12 << endl;
  // cout << "S21 matrix:\n" << S21 << endl;
  // cout << "S22 matrix:\n" << S22 << endl;

  // Eqn(14) and Eqn(15)
  Eigen::MatrixXd M = C1.inverse() * (S11 - S12 * (S22.inverse() * S21));
  // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_m(M);

  // Eigen::MatrixXd evec = es_m.eigenvectors();
  // Eigen::VectorXd eval = es_m.eigenvalues();

  Eigen::EigenSolver<Eigen::MatrixXd> es_m(M);
  Eigen::MatrixXd evec = es_m.eigenvectors().real();
  Eigen::VectorXd eval = es_m.eigenvalues().real();

  cout << "M:\n" << M << endl;
  cout << "M_evec:\n" << evec << endl;
  cout << "M_eval:\n" << eval.transpose() << endl;

  // Find the column index of the maximum eigenvalue
  int max_column_index;
  eval.maxCoeff(&max_column_index);

  cout << "max_column_index = " << max_column_index << endl;

  // Get the corresponding eigenvector
  Eigen::VectorXd u1 = evec.col(max_column_index);
  Eigen::VectorXd u2 = -(S22.inverse() * S21) * u1;

  // Concatenate u1 and u2 into a single vector u
  Eigen::VectorXd u(u1.size() + u2.size());
  u << u1, u2;

  // Output the results (for debugging purposes)
  cout << "u1 vector:\n" << u1.transpose() << endl;
  cout << "u2 vector:\n" << u2.transpose() << endl;
  cout << "u vector:\n" << u.transpose() << endl;
  //*****************************************************





  //*****************************************************
  // computation of the compensation model coefficients;
  double mag_earth_intensity = 54093.9956380105; // nT

  Eigen::Matrix<double, 10, 1> ellipsoid_coeffs(u);
  // ellipsoid_coeffs << -0.662597125572740, -0.527259510484714,
  //    -0.531301167116744, -0.00961069105634865, -0.00641802453979663,
  //    -0.0234087843548756, -7178.33182399735 * 2.0, -1396.12980160834 * 2.0,
  //    1742.16959639974 * 2.0, 1235304461.23394;

  double a = ellipsoid_coeffs[0];
  double b = ellipsoid_coeffs[1];
  double c = ellipsoid_coeffs[2];
  double f = ellipsoid_coeffs[3];
  double g = ellipsoid_coeffs[4];
  double h = ellipsoid_coeffs[5];
  double p = ellipsoid_coeffs[6];
  double q = ellipsoid_coeffs[7];
  double r = ellipsoid_coeffs[8];
  double d = ellipsoid_coeffs[9];

  Eigen::Matrix3d As_hat;
  As_hat << a, h, g, h, b, f, g, f, c;
  Eigen::Vector3d bs_hat;
  bs_hat << p, q, r;
  double cs_hat = d;
  std::cout << "As_hat = " << std::endl << As_hat << std::endl;
  std::cout << "As_hat.inverse() = " << std::endl
            << As_hat.inverse() << std::endl;
  std::cout << "bs_hat = " << bs_hat.transpose() << std::endl;
  std::cout << "cs_hat = " << cs_hat << std::endl;

  double den = bs_hat.dot(As_hat.inverse() * bs_hat) - cs_hat;
  double alpha = mag_earth_intensity * mag_earth_intensity / den;
  std::cout << "alpha = " << alpha << std::endl;

  Eigen::Vector3d o_hat = -As_hat.inverse() * bs_hat;
  std::cout << "o_hat = " << o_hat.transpose() << std::endl;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(As_hat);
  Eigen::Matrix3d As_evec = es.eigenvectors();
  Eigen::Vector3d As_eval = es.eigenvalues(); // increasing order;

  Eigen::Matrix3d sqrt_As_eval;
  sqrt_As_eval << sqrt(fabs(As_eval[0])), 0, 0, 0, sqrt(fabs(As_eval[1])), 0, 0,
      0, sqrt(fabs(As_eval[2]));

  Eigen::Matrix3d D_tilde_inv =
      sqrt(fabs(alpha)) * sqrt_As_eval * As_evec.transpose();
  std::cout << "D_tilde_inv = " << std::endl << D_tilde_inv << std::endl;
  //*****************************************************

  // Eigen::Matrix3d D_tilde;
  // D_tilde << 1.19555366931592, 0.205168980701854, 0.0710780020693875,
  //    0.152967675172370, -0.575201075233501, -0.912626201443581,
  //    -0.118363856179626, 0.891188576390446, -0.581528856436894;
  //
  // Eigen::Vector3d o_hat;
  // o_hat << -10788.1970344684, -2231.81509399570, 3449.75299982293;


  //*****************************************************
  // estimate initial value of orthogonal matrix R=V*R^{mb};


  //*****************************************************






  //*****************************************************
  // ceres optimization problem;
	ceres::Problem problem;
	ceres::LossFunction *loss_function=new ceres::HuberLoss(1.0);
  ceres::LocalParameterization *quat_param = new ceres::EigenQuaternionParameterization;

	Eigen::Matrix3d rotation_rc;
	Eigen::Vector3d mag_r,mag_c;
	mag_r(0)=x_m[0];
	mag_r(1)=y_m[0];
	mag_r(2)=z_m[0];
	mag_c(0)=x_m[1];
	mag_c(1)=y_m[1];
	mag_c(2)=z_m[1];

	ceres::CostFunction *cost_function=CostFunctionCreator::Create(D_tilde_inv,o_hat,rotation_rc,mag_r,mag_c);

	Eigen::Quaterniond quat;
	problem.AddResidualBlock(cost_function,loss_function,quat.coeffs().data());
	problem.SetParameterization(quat.coeffs().data(),quat_param);

  ceres::Solver::Summary summary;
  ceres::Solver::Options options;
  options.max_num_iterations = 10;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.minimizer_progress_to_stdout = true;

  //ceres::Solve(options, &problem, &summary);
  //if (minimizer_progress_to_stdout_) std::cout << summary.FullReport() << std::endl;

  //return summary.IsSolutionUsable();
  //*****************************************************

  return 0;
}
