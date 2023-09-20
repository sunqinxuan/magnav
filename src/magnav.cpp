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

/*
using namespace magnav;
using namespace std;

enum TLterm { PERMANENT, INDUCED, EDDY, FDM, BIAS };

void print(const auto &set) {
  for (const auto &elem : set)
    std::cout << elem << "; ";
  std::cout << '\n';
}


int main1(int argc, char *argv[]) {

  std::vector<double> b_coeff;
  //	b_coeff.push_back(4.0);
  //	b_coeff.push_back(5.0);
  //	b_coeff.push_back(6.0);
  vectord a_coeff;
  //	a_coeff.push_back(1.0);
  //	a_coeff.push_back(2.0);
  //	a_coeff.push_back(3.0);

  //coeff_ab(4, 0.1, 0.9, 10, a_coeff, b_coeff);

  b_coeff = {0.4328, 0, -1.7314, 0, 2.5971, 0, -1.7314, 0, 0.4328};

  a_coeff = {1.0000, 0.0000,  -2.3695, -0.0000, 2.3140,
             0.0000, -1.0547, -0.0000, 0.1874};

  cout << endl << "a= " << endl;
  print(a_coeff);
  cout << endl;

  cout << endl << "b= " << endl;
  print(b_coeff);
  cout << endl;

  vectord input_signal;
  for (int i = 0; i < 30; i++)
    input_signal.push_back(i);

  cout << endl << "x= " << endl;
  print(input_signal);
  cout << endl;

  //	input_signal.push_back(1);
  //	input_signal.push_back(2);
  //	input_signal.push_back(3);
  //	input_signal.push_back(4);
  //	input_signal.push_back(5);
  //	input_signal.push_back(6);
  //	input_signal.push_back(7);
  //	input_signal.push_back(8);

  vectord y_filtfilt_ori;
  y_filtfilt_ori.push_back(-6731884.25000000);
  y_filtfilt_ori.push_back(7501778.75000000);
  y_filtfilt_ori.push_back(-2757230.25000000);
  y_filtfilt_ori.push_back(-662443.250000000);
  y_filtfilt_ori.push_back(1360955.75000000);
  y_filtfilt_ori.push_back(-686678.250000000);
  y_filtfilt_ori.push_back(4135.75000000000);
  y_filtfilt_ori.push_back(227147.750000000);

  // vectord y_filter_out; vectord zi = { 0 };
  // filter(b_coeff, a_coeff, input_signal, y_filter_out, zi);
  // Assert::IsTrue(compare(y_filter_out, y_filter_ori, 0.0001));
  // assert(compare(y_filter_out, y_filter_ori, 0.0001));

  vectord y_filtfilt_out;
  filtfilt(b_coeff, a_coeff, input_signal, y_filtfilt_out);
  // Assert::IsTrue(compare(y_filtfilt_out, y_filtfilt_ori, 0.0001));
  // assert(compare(y_filter_out, y_filter_ori, 0.0001));

  std::cout << "y_filtfilt_out: " << endl;
  print(y_filtfilt_out);

  return 0;
}

int main2()
{
//	ifstream ifile;
//	ifile.open("E:\\HRdataset\\butterworth\\input.txt");

        vector<double> input, output;
        double fps = 10;

        const int N = 30;
        //
        //for (int i = 0; i < 1000; i++)
        //{
        //	float x;
        //	ifile >> x;
        //	input.push_back(x);
        //}

        //Frequency bands is a vector of values - Lower Frequency Band and
Higher Frequency Band

        //First value is lower cutoff and second value is higher cutoff
//	double FrequencyBands[2] = { 1.5/fps*2, 2.5/fps*2 };//these values are
as a ratio of f/fs, where fs is sampling rate, and f is cutoff frequency double
FrequencyBands[2] = { 0.1,0.9 };//these values are as a ratio of f/fs, where fs
is sampling rate, and f is cutoff frequency
        //and therefore should lie in the range [0 1]
        //Filter Order

        int FiltOrd = 4;

        //Pixel Time Series

        //Create the variables for the numerator and denominator coefficients
        vector<double> a;
        vector<double> b;
        //Pass Numerator Coefficients and Denominator Coefficients arrays into
function, will return the same

        vector<double> x(N);
        vector<double> y(N);

        for (int i = 0; i < N; i++)
        {
//		ifile >> x[i];
                x[i]=i;
        }

        //is A in matlab function and the numbers are correct
        a = ComputeDenCoeffs(FiltOrd, FrequencyBands[0], FrequencyBands[1]);
//	a.resize(a.size()-1);
        cout<<"a="<<endl;
        print(a);
//	for (int k = 0; k<a.size(); k++)
//	{
//		printf("DenC is: %lf\n", a[k]);
//	}

        b = ComputeNumCoeffs(FiltOrd, FrequencyBands[0], FrequencyBands[1], a);
        cout<<"b="<<endl;
        print(b);
//	for (int k = 0; k<b.size(); k++)
//	{
//		printf("NumC is: %lf\n", b[k]);
//	}

//	y = filter(x,b,a);
  filtfilt(b, a, x, y);

        cout<<"y="<<endl;
        print(y);

        a.resize(a.size()-1);
  filtfilt(b, a, x, y);

        cout<<"y="<<endl;
        print(y);

//	ofstream ofile;
//	ofile.open("E:\\HRdataset\\butterworth\\output.txt");
//	for (int i = 0; i < N; i++)
//	{
//		ofile << y[i] << endl;
//	}
//	ofile.close();

        cout<<"x="<<endl;
        print(x);

        int n=5;
        x.erase(x.begin(),x.begin()+n);
        x.erase(x.end()-n,x.end());

        cout<<"x trim ="<<endl;
        print(x);

   Eigen::MatrixXf A(3,3);
   Eigen::VectorXf B(3);
   A << 1,2,3,  4,5,6,  7,8,10;
         Eigen::MatrixXf AA=A.transpose()*A+4.0*Eigen::MatrixXf::Identity(3,3);
   B << 3, 3, 4;
   std::cout << "Here is the matrix A:\n" << AA << std::endl;
   std::cout << "Here is the vector B:\n" << B << std::endl;
   Eigen::Vector3f X = AA.colPivHouseholderQr().solve(B);
   std::cout << "The solution is:\n" << X << std::endl;
   Eigen::Vector3f XX = AA.llt().solve(B);
   std::cout << "The solution is:\n" << XX << std::endl;

        return 0;
}
*/

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
  // H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);
  // DataSet dataset = file.openDataSet(DATASET_NAME);
  //
  // DataSpace dataspace = dataset.getSpace();
  // int dimNums = dataspace.getSimpleExtentNdims();
  // cout << "dimNums= " << dimNums << endl;
  // hsize_t *dims = new hsize_t[dimNums];
  // dataspace.getSimpleExtentDims(dims);
  // for (int i = 0; i < dimNums; i++) {
  //  cout << i << "\t" << dims[i] << endl;
  //}
  //
  // int vds = 1;
  // for (size_t i = 0; i < dimNums; i++) {
  //  vds = vds * dims[i];
  //}
  // std::vector<double> data;
  // data.resize(vds);
  // dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
  //
  // std::ofstream fp("/home/sun/magnav/h5data.txt", std::ios::out);
  // for (int i = 0; i < data.size(); i++) {
  //  fp << data[i] << endl;
  //}
  // fp.close();

  H5Data h5data1002, h5data1003;

  getFlightData("/home/sun/magnav/data/Flt1002_train.h5", h5data1002);
  getFlightData("/home/sun/magnav/data/Flt1003_train.h5", h5data1003);

  // std::ofstream fp("/home/sun/magnav/h5data.txt", std::ios::out);
  // for (int i = 0; i < h5data1002.line.size(); i++) {
  //  fp << i << "\t" << h5data1002.line[i] << std::endl;
  //}
  // fp.close();

  MagCompensation mag_comp;
  H5Data out_data;
  mag_comp.compensate(out_data, h5data1002, h5data1003);

  //using namespace H5;
  //const H5std_string FILE_NAME("SDS.h5");
  //const H5std_string DATASET_NAME("IntArray");
  //const int NX = 5; // dataset dimensions
  //const int NY = 6;
  //const int RANK = 2;
  /*
   * Data initialization.
   */
  //int i, j;
  //int data[NX][NY]; // buffer for data to write
  //for (j = 0; j < NX; j++) {
  //  for (i = 0; i < NY; i++)
  //    data[j][i] = i + j;
  //}
  /*
   * 0 1 2 3 4 5
   * 1 2 3 4 5 6
   * 2 3 4 5 6 7
   * 3 4 5 6 7 8
   * 4 5 6 7 8 9
   */
  /*
   * Create a new file using H5F_ACC_TRUNC access,
   * default file creation properties, and default file
   * access properties.
   */
	H5::H5File file("data_TL.h5", H5F_ACC_TRUNC);

  /*
   * Define the size of the array and create the data space for fixed
   *
   * size dataset.
   */
  hsize_t dimsf[1]; // dataset dimensions
  dimsf[0] = out_data.mag_1_igrf.size();
  //dimsf[1] = NY;
	H5::DataSpace dataspace(1, dimsf);

  /*
   * Define datatype for the data in the file.
   * We will store little endian INT numbers.
   */
  //IntType datatype(PredType::NATIVE_INT);
	H5::DataType datatype(H5::PredType::NATIVE_DOUBLE);

	H5::DataSet dataset = file.createDataSet("tt", datatype, dataspace);
  dataset.write(out_data.tt.data(), H5::PredType::NATIVE_DOUBLE);

	dataset = file.createDataSet("slg", datatype, dataspace);
  dataset.write(out_data.mag_1_igrf.data(), H5::PredType::NATIVE_DOUBLE);

	dataset = file.createDataSet("mag_3_c", datatype, dataspace);
  dataset.write(out_data.mag_3_uc.data(), H5::PredType::NATIVE_DOUBLE);

	dataset = file.createDataSet("mag_4_c", datatype, dataspace);
  dataset.write(out_data.mag_4_uc.data(), H5::PredType::NATIVE_DOUBLE);

	dataset = file.createDataSet("mag_5_c", datatype, dataspace);
  dataset.write(out_data.mag_5_uc.data(), H5::PredType::NATIVE_DOUBLE);

  return 0;
}
