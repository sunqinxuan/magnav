/***********************************************
#
# Author: Sun Qinxuan
#
# Email: sunqinxuan@outlook.com
#
# Last modified:	2023-08-15 14:04
#
# Filename:		tolles_lawson.cpp
#
# Description:
#
************************************************/

#include "mag_compensation/tolles_lawson.hpp"

namespace magnav {

bool TollesLawson::createMatrixA_Vector(
    std::vector<std::vector<double>> &TL_A_, const std::vector<double> &Bx,
    const std::vector<double> &By, const std::vector<double> &Bz,
    const std::unordered_set<TLterm> &terms) {
  if (!(Bx.size() == By.size() && Bx.size() == Bz.size())) {
    ERROR("[TollesLawson][createMatrixA] Bx By Bz sizes not equal!");
    return false;
  }
  if (terms.empty()) {
    ERROR("[TollesLawson][createMatrixA] terms empty!");
    return false;
  }

  int N = Bx.size();

  std::vector<double> Bx_dot, By_dot, Bz_dot;
  fdm(Bx_dot, Bx);
  fdm(By_dot, By);
  fdm(Bz_dot, Bz);

  TL_A_.clear();

  if (terms.find(PERMANENT) != terms.end()) {
    TL_A_.resize(3);
    TL_A_[0].clear();
    TL_A_[1].clear();
    TL_A_[2].clear();
    for (int i = 0; i < N; i++) {
      TL_A_[0].push_back(1);
      TL_A_[0].push_back(0);
      TL_A_[0].push_back(0);

      TL_A_[1].push_back(0);
      TL_A_[1].push_back(1);
      TL_A_[1].push_back(0);

      TL_A_[2].push_back(0);
      TL_A_[2].push_back(0);
      TL_A_[2].push_back(1);
    }
  }

  if (terms.find(INDUCED) != terms.end()) {
    int idx = TL_A_.size();
    TL_A_.resize(idx + 6);
    TL_A_[idx].clear();
    TL_A_[idx + 1].clear();
    TL_A_[idx + 2].clear();
    TL_A_[idx + 3].clear();
    TL_A_[idx + 4].clear();
    TL_A_[idx + 5].clear();
    for (int i = 0; i < N; i++) {
      TL_A_[idx].push_back(Bx[i]);
      TL_A_[idx].push_back(0);
      TL_A_[idx].push_back(0);

      TL_A_[idx + 1].push_back(By[i]);
      TL_A_[idx + 1].push_back(Bx[i]);
      TL_A_[idx + 1].push_back(0);

      TL_A_[idx + 2].push_back(Bz[i]);
      TL_A_[idx + 2].push_back(0);
      TL_A_[idx + 2].push_back(Bx[i]);

      TL_A_[idx + 3].push_back(0);
      TL_A_[idx + 3].push_back(By[i]);
      TL_A_[idx + 3].push_back(0);

      TL_A_[idx + 4].push_back(0);
      TL_A_[idx + 4].push_back(Bz[i]);
      TL_A_[idx + 4].push_back(By[i]);

      TL_A_[idx + 5].push_back(0);
      TL_A_[idx + 5].push_back(0);
      TL_A_[idx + 5].push_back(Bz[i]);
    }
  }

  if (terms.find(EDDY) != terms.end()) {
    int idx = TL_A_.size();
    TL_A_.resize(idx + 9);
    TL_A_[idx].clear();
    TL_A_[idx + 1].clear();
    TL_A_[idx + 2].clear();
    TL_A_[idx + 3].clear();
    TL_A_[idx + 4].clear();
    TL_A_[idx + 5].clear();
    TL_A_[idx + 6].clear();
    TL_A_[idx + 7].clear();
    TL_A_[idx + 8].clear();
    for (int i = 0; i < N; i++) {
      TL_A_[idx].push_back(Bx_dot[i]);
      TL_A_[idx].push_back(0);
      TL_A_[idx].push_back(0);

      TL_A_[idx + 1].push_back(By_dot[i]);
      TL_A_[idx + 1].push_back(0);
      TL_A_[idx + 1].push_back(0);

      TL_A_[idx + 2].push_back(Bz_dot[i]);
      TL_A_[idx + 2].push_back(0);
      TL_A_[idx + 2].push_back(0);

      TL_A_[idx + 3].push_back(0);
      TL_A_[idx + 3].push_back(Bx_dot[i]);
      TL_A_[idx + 3].push_back(0);

      TL_A_[idx + 4].push_back(0);
      TL_A_[idx + 4].push_back(By_dot[i]);
      TL_A_[idx + 4].push_back(0);

      TL_A_[idx + 5].push_back(0);
      TL_A_[idx + 5].push_back(Bz_dot[i]);
      TL_A_[idx + 5].push_back(0);

      TL_A_[idx + 6].push_back(0);
      TL_A_[idx + 6].push_back(0);
      TL_A_[idx + 6].push_back(Bx_dot[i]);

      TL_A_[idx + 7].push_back(0);
      TL_A_[idx + 7].push_back(0);
      TL_A_[idx + 7].push_back(By_dot[i]);

      TL_A_[idx + 8].push_back(0);
      TL_A_[idx + 8].push_back(0);
      TL_A_[idx + 8].push_back(Bz_dot[i]);
    }
  }

  return true;
}

double TollesLawson::createCoeff_Vector(
    std::vector<double> &TL_beta_, const std::vector<double> &Bx,
    const std::vector<double> &By, const std::vector<double> &Bz,
    const std::vector<double> &Bex, const std::vector<double> &Bey,
    const std::vector<double> &Bez, const std::unordered_set<TLterm> &terms) {
  if (!(Bx.size() == By.size() && Bx.size() == Bz.size())) {
    ERROR("[TollesLawson][createCoeff] Bx By Bz sizes not equal!");
    return false;
  }
  if (terms.empty()) {
    ERROR("[TollesLawson][createCoeff] terms empty!");
    return false;
  }

  // if (Bt.empty()) {
  //  Bt.resize(Bx.size());
  //  for (size_t i = 0; i < Bx.size(); i++) {
  //    Bt[i] = sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]);
  //  }
  //} else {
  //  if (Bx.size() != Bt.size()) {
  //    ERROR("[TollesLawson][createCoeff] Bx Bt sizes not equal!");
  //    return false;
  //  }
  //}

  // bool perform_filter;
  // if ((pass1 > 0 && pass1 < fs / 2.0) || (pass2 > 0 && pass2 < fs / 2.0)) {
  //  perform_filter = true;
  //} else {
  //  perform_filter = false;
  //  INFO("not filtering (or trimming) Tolles-Lawson data.");
  //}

  std::vector<std::vector<double>> TL_A_;
  if (!createMatrixA_Vector(TL_A_, Bx, By, Bz, terms)) {
    ERROR("[TollesLawson][createCoeff] error calculating matrix TL_A_!");
    return false;
  }

  std::vector<double> Bt, Be;
  stack_vector(Bt, Bx, By, Bz);
  stack_vector(Be, Bex, Bey, Bez);

  std::vector<std::vector<double>> TL_A_filt_;
  std::vector<double> B_filt;

  for (int i = 0; i < Bt.size(); i++) {
    B_filt.push_back(Bt[i] - Be[i]);
  }
  TL_A_filt_ = TL_A_;

  // linear regression to get TL coefficients;
  std::vector<double> residual;
  linear_regression(TL_beta_, residual, B_filt, TL_A_filt_, 0);

  // compute TL fit error variance;
  double sum = std::accumulate(std::begin(residual), std::end(residual), 0.0);
  double mean = sum / residual.size();
  double variance = 0.0;
  for (size_t i = 0; i < residual.size(); i++) {
    variance = variance + pow(residual[i] - mean, 2);
  }
  variance = variance / residual.size();
  INFO("TL fit error residual: ", mean);
  INFO("TL fit error variance: ", variance, "\n");

  return variance;
}

bool TollesLawson::createMatrixA(std::vector<std::vector<double>> &TL_A_,
                                 const std::vector<double> &Bx,
                                 const std::vector<double> &By,
                                 const std::vector<double> &Bz,
                                 std::vector<double> &Bt,
                                 const std::unordered_set<TLterm> &terms,
                                 const double Bt_scale) {
  if (!(Bx.size() == By.size() && Bx.size() == Bz.size())) {
    ERROR("[TollesLawson][createMatrixA] Bx By Bz sizes not equal!");
    return false;
  }
  if (terms.empty()) {
    ERROR("[TollesLawson][createMatrixA] terms empty!");
    return false;
  }

  Bt.clear(); // TODO
  if (Bt.empty()) {
    Bt.resize(Bx.size());
    for (size_t i = 0; i < Bx.size(); i++) {
      Bt[i] = sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]);
    }
  } else {
    if (Bx.size() != Bt.size()) {
      ERROR("[TollesLawson][createMatrixA] Bx Bt sizes not equal!");
      return false;
    }
  }

  int N = Bx.size();
  std::vector<double> Bx_hat(N), By_hat(N), Bz_hat(N);
  for (int i = 0; i < N; i++) {
    Bx_hat[i] = Bx[i] / Bt[i];
    By_hat[i] = By[i] / Bt[i];
    Bz_hat[i] = Bz[i] / Bt[i];
  }

  std::vector<double> Bx_dot, By_dot, Bz_dot;
  // Bx_dot.clear();
  // By_dot.clear();
  // Bz_dot.clear();
  fdm(Bx_dot, Bx);
  fdm(By_dot, By);
  fdm(Bz_dot, Bz);

  // modified according to lw-kong/MagNav/src/create_TL_Amat.jl
  double Bt_sum = 0;
  for (int i = 0; i < N; i++) {
    Bt_sum += Bt[i];
  }
  Bt_sum = Bt_sum / double(N);
  double Bt_mean = (Bt_sum) > 1e-6 ? Bt_sum : 1e-6;

  // if use Bt_scale
  if (false) {
    Bt_mean = Bt_scale;
  }

  std::vector<double> Bx_hat_Bx(N), Bx_hat_By(N), Bx_hat_Bz(N), By_hat_By(N),
      By_hat_Bz(N), Bz_hat_Bz(N), Bx_hat_Bx_dot(N), Bx_hat_By_dot(N),
      Bx_hat_Bz_dot(N), By_hat_Bx_dot(N), By_hat_By_dot(N), By_hat_Bz_dot(N),
      Bz_hat_Bx_dot(N), Bz_hat_By_dot(N), Bz_hat_Bz_dot(N);
  for (int i = 0; i < N; i++) {
    Bx_hat_Bx[i] = Bx_hat[i] * Bx[i]; // / Bt_mean;
    Bx_hat_By[i] = Bx_hat[i] * By[i]; // / Bt_mean;
    Bx_hat_Bz[i] = Bx_hat[i] * Bz[i]; // / Bt_mean;
    By_hat_By[i] = By_hat[i] * By[i]; // / Bt_mean;
    By_hat_Bz[i] = By_hat[i] * Bz[i]; // / Bt_mean;
    Bz_hat_Bz[i] = Bz_hat[i] * Bz[i]; // / Bt_mean;

    Bx_hat_Bx_dot[i] = Bx_hat[i] * Bx_dot[i]; // / Bt_mean;
    Bx_hat_By_dot[i] = Bx_hat[i] * By_dot[i]; // / Bt_mean;
    Bx_hat_Bz_dot[i] = Bx_hat[i] * Bz_dot[i]; // / Bt_mean;
    By_hat_Bx_dot[i] = By_hat[i] * Bx_dot[i]; // / Bt_mean;
    By_hat_By_dot[i] = By_hat[i] * By_dot[i]; // / Bt_mean;
    By_hat_Bz_dot[i] = By_hat[i] * Bz_dot[i]; // / Bt_mean;
    Bz_hat_Bx_dot[i] = Bz_hat[i] * Bx_dot[i]; // / Bt_mean;
    Bz_hat_By_dot[i] = Bz_hat[i] * By_dot[i]; // / Bt_mean;
    Bz_hat_Bz_dot[i] = Bz_hat[i] * Bz_dot[i]; // / Bt_mean;
  }

  TL_A_.clear();

  if (terms.find(PERMANENT) != terms.end()) {
    TL_A_.resize(3);
    TL_A_[0] = Bx_hat;
    TL_A_[1] = By_hat;
    TL_A_[2] = Bz_hat;
  }

  if (terms.find(INDUCED) != terms.end()) {
    int idx = TL_A_.size();
    TL_A_.resize(idx + 6);
    TL_A_[idx] = Bx_hat_Bx;
    TL_A_[idx + 1] = Bx_hat_By;
    TL_A_[idx + 2] = Bx_hat_Bz;
    TL_A_[idx + 3] = By_hat_By;
    TL_A_[idx + 4] = By_hat_Bz;
    TL_A_[idx + 5] = Bz_hat_Bz;
  }

  if (terms.find(EDDY) != terms.end()) {
    int idx = TL_A_.size();
    TL_A_.resize(idx + 9);
    TL_A_[idx] = Bx_hat_Bx_dot;
    TL_A_[idx + 1] = Bx_hat_By_dot;
    TL_A_[idx + 2] = Bx_hat_Bz_dot;
    TL_A_[idx + 3] = By_hat_Bx_dot;
    TL_A_[idx + 4] = By_hat_By_dot;
    TL_A_[idx + 5] = By_hat_Bz_dot;
    TL_A_[idx + 6] = Bz_hat_Bx_dot;
    TL_A_[idx + 7] = Bz_hat_By_dot;
    TL_A_[idx + 8] = Bz_hat_Bz_dot;
  }

  if (terms.find(FDM) != terms.end()) {
    int idx = TL_A_.size();
    TL_A_.resize(idx + 3);
    TL_A_[idx] = Bx_dot;
    TL_A_[idx + 1] = By_dot;
    TL_A_[idx + 2] = Bz_dot;
  }

  if (terms.find(BIAS) != terms.end()) {
    int idx = TL_A_.size();
    std::vector<double> tmp(N, 1);
    TL_A_.resize(idx + 1, tmp);
  }

  return true;
}

double TollesLawson::createCoeff(
    std::vector<double> &TL_beta_, const std::vector<double> &Bx,
    const std::vector<double> &By, const std::vector<double> &Bz,
    const std::vector<double> &B, const std::vector<double> &Be,
    std::vector<double> &Bt, const std::unordered_set<TLterm> &terms,
    const double lambda, const double pass1, const double pass2,
    const double fs, const int pole, const int trim, const double Bt_scale) {
  if (!(Bx.size() == By.size() && Bx.size() == Bz.size() &&
        Bx.size() == B.size())) {
    ERROR("[TollesLawson][createCoeff] Bx By Bz B sizes not equal!");
    return false;
  }
  if (terms.empty()) {
    ERROR("[TollesLawson][createCoeff] terms empty!");
    return false;
  }

  if (Bt.empty()) {
    Bt.resize(Bx.size());
    for (size_t i = 0; i < Bx.size(); i++) {
      Bt[i] = sqrt(Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i]);
    }
  } else {
    if (Bx.size() != Bt.size()) {
      ERROR("[TollesLawson][createCoeff] Bx Bt sizes not equal!");
      return false;
    }
  }

  bool perform_filter;
  if ((pass1 > 0 && pass1 < fs / 2.0) || (pass2 > 0 && pass2 < fs / 2.0)) {
    perform_filter = true;
  } else {
    perform_filter = false;
    INFO("not filtering (or trimming) Tolles-Lawson data.");
  }

  std::vector<std::vector<double>> TL_A_;
  if (!createMatrixA(TL_A_, Bx, By, Bz, Bt, terms, Bt_scale)) {
    ERROR("[TollesLawson] error calculating matrix TL_A_!");
    return false;
  }

  std::vector<std::vector<double>> TL_A_filt_;
  std::vector<double> B_filt;

  if (false) {
    if (perform_filter) {
      // bandpass filter;
      // filter columns of matrix A and the measurements B
      // and trim edges;
      TL_A_filt_.resize(TL_A_.size());
      for (size_t i = 0; i < TL_A_.size(); i++) {
        bwbp_filter(TL_A_filt_[i], TL_A_[i], pass1, pass2, pole, trim);
      }
      bwbp_filter(B_filt, B, pass1, pass2, pole, trim);
    } else {
      TL_A_filt_ = TL_A_;
      B_filt = B;
    }
  } else {
    for (int i = 0; i < B.size(); i++) {
      B_filt.push_back(B[i] - Be[i]);
    }
    TL_A_filt_ = TL_A_;
  }

  // linear regression to get TL coefficients;
  std::vector<double> residual;
  linear_regression(TL_beta_, residual, B_filt, TL_A_filt_, lambda);

  // compute TL fit error variance;
  double sum = std::accumulate(std::begin(residual), std::end(residual), 0.0);
  double mean = sum / residual.size();
  double variance = 0.0;
  for (size_t i = 0; i < residual.size(); i++) {
    variance = variance + pow(residual[i] - mean, 2);
  }
  variance = variance / residual.size();
  INFO("TL fit error residual: ", mean);
  INFO("TL fit error variance: ", variance, "\n");

  return variance;
}

// sum((itr .- mean(itr)).^2) / (length(itr) - 1)

bool TollesLawson::fdm(std::vector<double> &dif, const std::vector<double> &x,
                       const TollesLawson::FDMscheme scheme) {
  // N = length(x)
  int N = x.size();
  //		std::vector<double> dif(N);
  dif.resize(N);

  if (scheme == BACKWARD && N > 1) {
    dif[0] = x[1] - x[0];
    for (int i = 1; i < N; i++) {
      dif[i] = x[i] - x[i - 1];
    }
    return true;
  } else if (scheme == FORWARD && N > 1) {
    dif[N - 1] = x[N - 1] - x[N - 2];
    for (int i = 0; i < N - 1; i++) {
      dif[i] = x[i + 1] - x[i];
    }
    return true;
  } else if (scheme == CENTRAL && N > 2) {
    dif[0] = x[1] - x[0];
    dif[N - 1] = x[N] - x[N - 1];
    for (int i = 1; i < N - 1; i++) {
      dif[i] = 0.5 * (x[i + 1] - x[i - 1]);
    }
    return true;
  } else if (scheme == BACKWARD2 && N > 3) {
    dif[0] = x[1] - x[0];
    dif[1] = x[2] - x[1];
    for (int i = 2; i < N; i++) {
      dif[i] = 0.5 * (3.0 * x[i] - 4.0 * x[i - 1] + x[i - 2]);
    }
    return true;
  } else if (scheme == FORWARD2 && N > 3) {
    dif[N - 1] = x[N - 1] - x[N - 2];
    dif[N - 2] = x[N - 2] - x[N - 3];
    for (int i = 0; i < N - 2; i++) {
      dif[i] = 0.5 * (-x[i + 2] + 4.0 * x[i + 1] - 3.0 * x[i]);
    }
    return true;
  } else if (scheme == FOURTH && N > 4) {
    dif[0] = 0;
    dif[1] = 0;
    dif[N - 1] = 0;
    dif[N - 2] = 0;
    for (int i = 2; i < N - 2; i++) {
      dif[i] =
          (x[i - 2] - 4.0 * x[i - 1] + 6.0 * x[i] - 4.0 * x[i + 1] + x[i + 2]) /
          16.0;
    }
    return true;
  } else {
    std::fill(dif.begin(), dif.end(), 0.0);
    return false;
  }
}

bool TollesLawson::bwbp_filter(std::vector<double> &y,
                               const std::vector<double> &x, const double pass1,
                               const double pass2, const int pole,
                               const int trim) {
  std::vector<double> bwbpfilter_A = ComputeDenCoeffs(pole, pass1, pass2);
  std::vector<double> bwbpfilter_B =
      ComputeNumCoeffs(pole, pass1, pass2, bwbpfilter_A);

  if (x.size() <= 3 * bwbpfilter_A.size() ||
      x.size() <= 3 * bwbpfilter_B.size()) {
    ERROR("[TollesLawson][bwbp_filter] input data too short!");
    return false;
  }

  // TODO
  // to be consistent with the resultes of 'butter' function in matlab;
  bwbpfilter_A.resize(bwbpfilter_A.size() - 1);

  y.clear();
  filtfilt(bwbpfilter_B, bwbpfilter_A, x, y);

  y.erase(y.begin(), y.begin() + trim);
  y.erase(y.end() - trim, y.end());

  return true;
}

void TollesLawson::linear_regression(std::vector<double> &coeff,
                                     std::vector<double> &residual,
                                     const std::vector<double> &bb,
                                     const std::vector<std::vector<double>> &AA,
                                     const double lambda) {
  const int N = bb.size();
  const int M = AA.size();
  Eigen::MatrixXd A(N, M);
  Eigen::VectorXd b(N);
  for (int i = 0; i < N; i++) {
    b(i) = bb[i];
    for (int j = 0; j < M; j++) {
      A(i, j) = AA[j][i];
    }
  }
  Eigen::MatrixXd ATA =
      A.transpose() * A; // + lambda * Eigen::MatrixXd::Identity(M, M);
  // std::cout << "ATA = " << std::endl << ATA << std::endl;
  Eigen::MatrixXd ATb = A.transpose() * b;
  // std::cout << "ATb = " << std::endl << ATb.transpose() << std::endl;
  Eigen::VectorXd x = ATA.llt().solve(ATb);
  // std::cout << "x = " << std::endl << x.transpose() << std::endl;

  coeff.resize(M);
  for (int i = 0; i < M; i++) {
    coeff[i] = x(i);
  }

  // linear regression fit error;
  Eigen::VectorXd delta = b - A * x;
  residual.resize(N);
  for (int i = 0; i < N; i++) {
    residual[i] = delta(i);
  }
}

void TollesLawson::stack_vector(std::vector<double> &result,
                                const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<double> &z) {
  if (!(x.size() == y.size() && x.size() == z.size())) {
    ERROR("[TollesLawson][stack_vector] x y z sizes not equal!");
    return;
  }
  result.clear();
  for (size_t i = 0; i < x.size(); i++) {
    result.push_back(x[i]);
    result.push_back(y[i]);
    result.push_back(z[i]);
  }
}

//
// void TollesLawson::getFilterCoeffAB(int n, double lowcut, double highcut,
//                                    int fs, std::vector<double> &acof_vec,
//                                    std::vector<double> &bcof_vec) {
//
//  double nyq = 0.5 * fs;
//  double f1f = lowcut / nyq;
//  double f2f = highcut / nyq;
//  double sf; // scaling factor
//
//  double *acof;
//  int *bcof;
//  /* calculate the d coefficients */
//  acof = dcof_bwbp(n, f1f, f2f);
//  if (acof == NULL) {
//    perror("Unable to calculate d coefficients");
//  }
//
//  /* calculate the c coefficients */
//  bcof = ccof_bwbp(n);
//  if (bcof == NULL) {
//    perror("Unable to calculate c coefficients");
//  }
//
//  sf = sf_bwbp(n, f1f, f2f); /* scaling factor for the c coefficients */
//
//  /* create the filter coefficient file */
//  for (int i = 0; i <= 2 * n; ++i) {
//    bcof_vec.push_back((double)bcof[i] * sf);
//  }
//
//  for (int i = 0; i <= 2 * n; ++i)
//    acof_vec.push_back(acof[i]);
//}

} // namespace magnav
