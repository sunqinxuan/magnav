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
// using std::cout;
// using std::endl;

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
  std::cout << "N=" << h5data.N << std::endl;
  H5::DataSet data_dt = file.openDataSet("dt");
  data_dt.read(&h5data.dt, H5::PredType::NATIVE_DOUBLE);
  std::cout << "dt=" << h5data.dt << std::endl;

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

int main(void) {
  H5Data h5data1002, h5data1003;

  getFlightData("/home/sun/magnav/data/Flt1002_train.h5", h5data1002);
  getFlightData("/home/sun/magnav/data/Flt1003_train.h5", h5data1003);

  MagCompensation mag_comp;
  H5Data out_data;
  //mag_comp.compensate(out_data, h5data1002, h5data1003);
  //mag_comp.compensateVector(out_data, h5data1002, h5data1003);
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
