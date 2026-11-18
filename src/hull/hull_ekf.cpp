#include "localization_v2/hull_fusion.hpp"

#include <cmath>

namespace naiad
{
namespace localization_v2
{
namespace
{

constexpr int kN = 20;
constexpr int kM = 4;
constexpr int kWheelScale[kM] = {0, 1, 3, 4};
constexpr double kPi = 3.1415926535897932;
constexpr double kDeg2Rad = kPi / 180.0;
constexpr double kMeasVarRear = 0.1 * 0.1;
constexpr double kMeasVarSmall = 0.05 * 0.05;
// 仿真 N1 FOG：ARW 0.001 °/√h，陀螺零偏 0.01 °/h，加表零偏 0.0002 m/s²，相关时间 1 h。
constexpr double kFogArwDegPerSqrtH = 0.001;
constexpr double kFogGyroBiasDegPerH = 0.01;
constexpr double kFogAccelBias = 0.0002;
constexpr double kFogBiasTauSec = 3600.0;

double meas_var(int meas)
{
  return meas >= 2 ? kMeasVarSmall : kMeasVarRear;
}

/// 杨元喜 IGG3。λ 乘在 R 上，与 rtkpos.cc 的 igg3_lambda 相同。|vt|≥k1 时用 1e6 等价于剔除。
double igg3_lambda(double vt, double k0, double k1)
{
  vt = std::fabs(vt);
  if (vt <= k0) {
    return 1.0;
  }
  if (k1 <= k0 || vt >= k1) {
    return 1.0e6;
  }
  const double span = k1 - vt;
  if (span < 1.0e-6) {
    return 1.0e6;
  }
  const double slope = (k1 - k0) / span;
  return (vt / k0) * slope * slope;
}

void zero(double * a, int n)
{
  for (int i = 0; i < n; ++i) {
    a[i] = 0.0;
  }
}

void mul3(const double * a, const double * b, double * c)
{
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      c[row + 3 * col] = a[row] * b[3 * col] + a[row + 3] * b[1 + 3 * col] +
                         a[row + 6] * b[2 + 3 * col];
    }
  }
}

void mul_vec3(const double * a, const double * x, double * y)
{
  for (int row = 0; row < 3; ++row) {
    y[row] = a[row] * x[0] + a[row + 3] * x[1] + a[row + 6] * x[2];
  }
}

void transpose3(const double * a, double * b)
{
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      b[row + 3 * col] = a[col + 3 * row];
    }
  }
}

void mul_nn(const double * a, const double * b, double * c)
{
  for (int col = 0; col < kN; ++col) {
    for (int row = 0; row < kN; ++row) {
      double sum = 0.0;
      for (int k = 0; k < kN; ++k) {
        sum += a[row + kN * k] * b[k + kN * col];
      }
      c[row + kN * col] = sum;
    }
  }
}

void transpose_n(const double * a, double * b)
{
  for (int col = 0; col < kN; ++col) {
    for (int row = 0; row < kN; ++row) {
      b[row + kN * col] = a[col + kN * row];
    }
  }
}

void right_perturb(const double * c, const double * dphi, double * out)
{
  const double dx = dphi[0];
  const double dy = dphi[1];
  const double dz = dphi[2];
  const double m[9] = {1.0, dz, -dy, -dz, 1.0, dx, dy, -dx, 1.0};
  mul3(c, m, out);
}

void orthonormalize(double * c)
{
  double c0[3] = {c[0], c[1], c[2]};
  double c1[3] = {c[3], c[4], c[5]};
  double n0 = std::sqrt(c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2]);
  if (n0 < 1e-15) {
    return;
  }
  for (int i = 0; i < 3; ++i) {
    c0[i] /= n0;
  }
  const double dot = c1[0] * c0[0] + c1[1] * c0[1] + c1[2] * c0[2];
  for (int i = 0; i < 3; ++i) {
    c1[i] -= dot * c0[i];
  }
  double n1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
  if (n1 < 1e-15) {
    return;
  }
  for (int i = 0; i < 3; ++i) {
    c1[i] /= n1;
  }
  const double c2[3] = {
    c0[1] * c1[2] - c0[2] * c1[1], c0[2] * c1[0] - c0[0] * c1[2],
    c0[0] * c1[1] - c0[1] * c1[0]};
  for (int i = 0; i < 3; ++i) {
    c[i] = c0[i];
    c[3 + i] = c1[i];
    c[6 + i] = c2[i];
  }
}

bool invert_meas(const double * a, double * inv)
{
  constexpr int kAug = 2 * kM;
  double m[kM][kAug];
  for (int row = 0; row < kM; ++row) {
    for (int col = 0; col < kM; ++col) {
      m[row][col] = a[row * kM + col];
      m[row][kM + col] = row == col ? 1.0 : 0.0;
    }
  }
  for (int col = 0; col < kM; ++col) {
    int pivot = col;
    for (int row = col + 1; row < kM; ++row) {
      if (std::fabs(m[row][col]) > std::fabs(m[pivot][col])) {
        pivot = row;
      }
    }
    if (std::fabs(m[pivot][col]) < 1e-12) {
      return false;
    }
    if (pivot != col) {
      for (int j = 0; j < kAug; ++j) {
        const double tmp = m[col][j];
        m[col][j] = m[pivot][j];
        m[pivot][j] = tmp;
      }
    }
    const double diag = m[col][col];
    for (int j = 0; j < kAug; ++j) {
      m[col][j] /= diag;
    }
    for (int row = 0; row < kM; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = m[row][col];
      for (int j = 0; j < kAug; ++j) {
        m[row][j] -= factor * m[col][j];
      }
    }
  }
  for (int row = 0; row < kM; ++row) {
    for (int col = 0; col < kM; ++col) {
      inv[row * kM + col] = m[row][kM + col];
    }
  }
  return true;
}

void body_velocity(const double * c_bs, const double * v_s, double * v_b)
{
  for (int axis = 0; axis < 3; ++axis) {
    v_b[axis] = c_bs[3 * axis] * v_s[0] + c_bs[3 * axis + 1] * v_s[1] +
                c_bs[3 * axis + 2] * v_s[2];
  }
}

/// 一个轮子的前向观测。gain = 1+s_i，比例因子列写 v_x。姿态三列为 0。
void fill_wheel_row(
  double * h, double * jacobian, int row_index, int scale_index, const double * c_bs,
  const double * v_b, const double * scale)
{
  const double gain = 1.0 + scale[scale_index];
  h[row_index] = gain * v_b[0];
  double * row = jacobian + row_index * kN;
  row[3] = gain * c_bs[0];
  row[4] = gain * c_bs[1];
  row[5] = gain * c_bs[2];
  row[15 + scale_index] = v_b[0];
}

/// H 行主序 4×20。顺序是大轮右、大轮左、小轮右、小轮左，对应 s_0、s_1、s_3、s_4。前轮不在观测里。
void wheel_observation(
  const double * c_bs, const double * v_s, const double * scale, double * h, double * jacobian)
{
  double v_b[3];
  body_velocity(c_bs, v_s, v_b);
  zero(h, kM);
  zero(jacobian, kM * kN);
  for (int row = 0; row < kM; ++row) {
    fill_wheel_row(h, jacobian, row, kWheelScale[row], c_bs, v_b, scale);
  }
}

constexpr double kNhcVarY = 0.02 * 0.02;
constexpr double kNhcVarZ = 0.02 * 0.02;

/// 模式 A 的一条 NHC。axis 为 1（侧向，后轮中心）或 2（法向，IMU）。量测为 0。
/// 两行姿态三列都是 0，只改对应的速度。侧向预测是 e_y·(v_imu + ω×l)，零偏列是 [l×] 的右向行。
void nhc_row(
  const double * c_bs, const double * v_b, int axis, const double * omega_b,
  const double * lever_frd, double * row, double * prediction)
{
  zero(row, kN);
  *prediction = v_b[axis];
  if (axis == 1) {
    const double wxr_y = omega_b[2] * lever_frd[0] - omega_b[0] * lever_frd[2];
    *prediction += wxr_y;
    row[12] = lever_frd[2];
    row[14] = -lever_frd[0];
  }
  row[3] = c_bs[3 * axis];
  row[4] = c_bs[1 + 3 * axis];
  row[5] = c_bs[2 + 3 * axis];
}

}  // namespace

void HullFusion::reset_error_covariance()
{
  zero(P_, kN * kN);
  const double att = std::pow(5.0 * kDeg2Rad, 2.0);
  const double vel = 0.5 * 0.5;
  const double pos = 1.0;
  const double acc_bias = kFogAccelBias * kFogAccelBias;
  const double gyro_bias_rad = kFogGyroBiasDegPerH * kDeg2Rad / 3600.0;
  const double gyro_bias = gyro_bias_rad * gyro_bias_rad;
  const double scale = 0.05 * 0.05;
  for (int i = 0; i < 3; ++i) {
    P_[i + kN * i] = att;
    P_[(3 + i) + kN * (3 + i)] = vel;
    P_[(6 + i) + kN * (6 + i)] = pos;
    P_[(9 + i) + kN * (9 + i)] = acc_bias;
    P_[(12 + i) + kN * (12 + i)] = gyro_bias;
  }
  for (int i = 0; i < 5; ++i) {
    P_[(15 + i) + kN * (15 + i)] = scale;
  }
}

void HullFusion::set_wheel_aiding(bool enabled)
{
  if (enabled && !wheel_aiding_) {
    reset_error_covariance();
  }
  wheel_aiding_ = enabled;
}

void HullFusion::set_update_hz(double wheel_hz, double nhc_hz)
{
  wheel_period_sec_ = 1.0 / wheel_hz;
  nhc_period_sec_ = 1.0 / nhc_hz;
}

void HullFusion::set_wheel_igg3(double k0, double k1)
{
  igg3_k0_ = k0;
  igg3_k1_ = k1;
}

bool HullFusion::wheel_igg3_holds()
{
  if (std::fabs(igg3_lambda(0.0, 5.0, 15.0) - 1.0) > 1e-12) {
    return false;
  }
  if (std::fabs(igg3_lambda(5.0, 5.0, 15.0) - 1.0) > 1e-12) {
    return false;
  }
  // |vt|=10：(|vt|/k0)*((k1-k0)/(k1-|vt|))^2 = 2 * 2^2 = 8
  if (std::fabs(igg3_lambda(10.0, 5.0, 15.0) - 8.0) > 1e-9) {
    return false;
  }
  if (std::fabs(igg3_lambda(-10.0, 5.0, 15.0) - 8.0) > 1e-9) {
    return false;
  }
  if (igg3_lambda(15.0, 5.0, 15.0) < 1.0e5 || igg3_lambda(20.0, 5.0, 15.0) < 1.0e5) {
    return false;
  }
  return true;
}

void HullFusion::propagate_error(double dt, const double * f_b)
{
  double c_se_t[9];
  transpose3(Cse_, c_se_t);
  double c_bs[9];
  mul3(c_se_t, Cbe_, c_bs);

  double phi[kN * kN];
  zero(phi, kN * kN);
  for (int i = 0; i < kN; ++i) {
    phi[i + kN * i] = 1.0;
  }
  // δφ 是 C_b^s 的右扰动。积分先减零偏、反馈再加回去，所以 φ̇=-δbg，δv̇=-C[f×]δφ-C δba。
  const double sk[9] = {
    0.0, f_b[2], -f_b[1], -f_b[2], 0.0, f_b[0], f_b[1], -f_b[0], 0.0};
  double c_sk[9];
  mul3(c_bs, sk, c_sk);
  for (int col = 0; col < 3; ++col) {
    phi[col + kN * (12 + col)] += -dt;
    for (int row = 0; row < 3; ++row) {
      phi[(3 + row) + kN * (9 + col)] += -c_bs[row + 3 * col] * dt;
      phi[(3 + row) + kN * col] += -c_sk[row + 3 * col] * dt;
    }
    phi[(6 + col) + kN * (3 + col)] += dt;
  }

  double tmp[kN * kN];
  double phi_t[kN * kN];
  mul_nn(phi, P_, tmp);
  transpose_n(phi, phi_t);
  mul_nn(tmp, phi_t, P_);

  const double arw = kFogArwDegPerSqrtH * kDeg2Rad / 60.0;
  const double q_att = arw * arw * dt;
  const double q_vel = 0.05 * 0.05 * dt;
  const double gyro_bias_rad = kFogGyroBiasDegPerH * kDeg2Rad / 3600.0;
  const double q_ba = (2.0 * kFogAccelBias * kFogAccelBias / kFogBiasTauSec) * dt;
  const double q_bg = (2.0 * gyro_bias_rad * gyro_bias_rad / kFogBiasTauSec) * dt;
  const double q_s = 1e-8 * dt;
  for (int i = 0; i < 3; ++i) {
    P_[i + kN * i] += q_att;
    P_[(3 + i) + kN * (3 + i)] += q_vel;
    P_[(9 + i) + kN * (9 + i)] += q_ba;
    P_[(12 + i) + kN * (12 + i)] += q_bg;
  }
  for (int i = 0; i < 5; ++i) {
    P_[(15 + i) + kN * (15 + i)] += q_s;
  }
  for (int col = 0; col < kN; ++col) {
    for (int row = col + 1; row < kN; ++row) {
      const double sym = 0.5 * (P_[row + kN * col] + P_[col + kN * row]);
      P_[row + kN * col] = sym;
      P_[col + kN * row] = sym;
    }
  }
}

void HullFusion::feedback_error(const double * dx)
{
  double corrected[9];
  right_perturb(Cbe_, dx, corrected);
  for (int i = 0; i < 9; ++i) {
    Cbe_[i] = corrected[i];
  }
  orthonormalize(Cbe_);

  double dv_e[3];
  double dp_e[3];
  mul_vec3(Cse_, dx + 3, dv_e);
  mul_vec3(Cse_, dx + 6, dp_e);
  for (int i = 0; i < 3; ++i) {
    ve_[i] += dv_e[i];
    re_[i] += dp_e[i];
    ba_[i] += dx[9 + i];
    bg_[i] += dx[12 + i];
  }
  for (int i = 0; i < 5; ++i) {
    scale_[i] += dx[15 + i];
  }
}

bool HullFusion::update_wheels(double t_sec, const WheelSpeedAtFog & speed)
{
  if (!armed_ || !wheel_aiding_ || !speed.valid) {
    return false;
  }
  if (wheel_sec_ < 0.0) {
    wheel_sec_ = t_sec;
    return false;
  }
  // 周期来自 wheel_fusion_hz。留 1 ms，避免时间戳的舍入把整拍丢掉。
  if (t_sec < wheel_sec_ + wheel_period_sec_ - 1e-3) {
    return false;
  }
  wheel_sec_ = t_sec;
  double c_se_t[9];
  transpose3(Cse_, c_se_t);
  double c_bs[9];
  mul3(c_se_t, Cbe_, c_bs);
  double v_s[3];
  mul_vec3(c_se_t, ve_, v_s);

  double h[kM];
  double jacobian[kM * kN];
  wheel_observation(c_bs, v_s, scale_, h, jacobian);
  // 四个轮子各自补到 IMU 后的前向速度。顺序与 kWheelScale 相同。
  const double z[kM] = {
    speed.motor_right[0], speed.motor_left[0], speed.small_right[0], speed.small_left[0]};

  double pht[kN * kM];
  for (int meas = 0; meas < kM; ++meas) {
    for (int row = 0; row < kN; ++row) {
      double sum = 0.0;
      for (int col = 0; col < kN; ++col) {
        sum += P_[row + kN * col] * jacobian[meas * kN + col];
      }
      pht[row + kN * meas] = sum;
    }
  }
  double r_var[kM];
  for (int meas = 0; meas < kM; ++meas) {
    r_var[meas] = meas_var(meas);
  }
  double innov_cov[kM * kM];
  for (int row = 0; row < kM; ++row) {
    for (int col = 0; col < kM; ++col) {
      double sum = row == col ? r_var[row] : 0.0;
      for (int k = 0; k < kN; ++k) {
        sum += jacobian[row * kN + k] * pht[k + kN * col];
      }
      innov_cov[row * kM + col] = sum;
    }
  }
  // 预测新息方差用膨胀前的 R。标准化残差相对 INS 预测速度，再把 λ 乘回该轮的 R。
  if (igg3_k0_ > 0.0) {
    for (int meas = 0; meas < kM; ++meas) {
      const double qii = innov_cov[meas * kM + meas];
      if (!(qii > 1.0e-12)) {
        continue;
      }
      const double vt = std::fabs(z[meas] - h[meas]) / std::sqrt(qii);
      const double lambda = igg3_lambda(vt, igg3_k0_, igg3_k1_);
      const double inflated = r_var[meas] * lambda;
      innov_cov[meas * kM + meas] += inflated - r_var[meas];
      r_var[meas] = inflated;
    }
  }
  double inv_s[kM * kM];
  if (!invert_meas(innov_cov, inv_s)) {
    return false;
  }
  double gain[kN * kM];
  for (int meas = 0; meas < kM; ++meas) {
    for (int row = 0; row < kN; ++row) {
      double sum = 0.0;
      for (int k = 0; k < kM; ++k) {
        sum += pht[row + kN * k] * inv_s[k * kM + meas];
      }
      gain[row + kN * meas] = sum;
    }
  }
  double dx[kN];
  for (int row = 0; row < kN; ++row) {
    double sum = 0.0;
    for (int meas = 0; meas < kM; ++meas) {
      sum += gain[row + kN * meas] * (z[meas] - h[meas]);
    }
    dx[row] = sum;
  }

  double kh[kN * kN];
  for (int col = 0; col < kN; ++col) {
    for (int row = 0; row < kN; ++row) {
      double sum = row == col ? 1.0 : 0.0;
      for (int meas = 0; meas < kM; ++meas) {
        sum -= gain[row + kN * meas] * jacobian[meas * kN + col];
      }
      kh[row + kN * col] = sum;
    }
  }
  double tmp[kN * kN];
  double kh_t[kN * kN];
  mul_nn(kh, P_, tmp);
  transpose_n(kh, kh_t);
  mul_nn(tmp, kh_t, P_);
  for (int col = 0; col < kN; ++col) {
    for (int row = 0; row < kN; ++row) {
      double extra = 0.0;
      for (int meas = 0; meas < kM; ++meas) {
        extra += gain[row + kN * meas] * r_var[meas] * gain[col + kN * meas];
      }
      P_[row + kN * col] += extra;
    }
  }
  for (int col = 0; col < kN; ++col) {
    for (int row = col + 1; row < kN; ++row) {
      const double sym = 0.5 * (P_[row + kN * col] + P_[col + kN * row]);
      P_[row + kN * col] = sym;
      P_[col + kN * row] = sym;
    }
  }

  feedback_error(dx);
  project_state();
  return true;
}

bool HullFusion::apply_scalar_update(
  const double * row, double prediction, double measurement, double variance)
{
  double ph[kN];
  for (int state = 0; state < kN; ++state) {
    double sum = 0.0;
    for (int col = 0; col < kN; ++col) {
      sum += P_[state + kN * col] * row[col];
    }
    ph[state] = sum;
  }
  double innov_var = variance;
  for (int state = 0; state < kN; ++state) {
    innov_var += row[state] * ph[state];
  }
  if (!(innov_var > 1e-12)) {
    return false;
  }
  const double innov = measurement - prediction;
  double dx[kN];
  for (int state = 0; state < kN; ++state) {
    dx[state] = ph[state] / innov_var * innov;
  }

  double kh[kN * kN];
  for (int col = 0; col < kN; ++col) {
    for (int state = 0; state < kN; ++state) {
      kh[state + kN * col] = (state == col ? 1.0 : 0.0) - ph[state] / innov_var * row[col];
    }
  }
  double tmp[kN * kN];
  double kh_t[kN * kN];
  mul_nn(kh, P_, tmp);
  transpose_n(kh, kh_t);
  mul_nn(tmp, kh_t, P_);
  for (int col = 0; col < kN; ++col) {
    for (int state = 0; state < kN; ++state) {
      P_[state + kN * col] += (ph[state] / innov_var) * variance * (ph[col] / innov_var);
    }
  }
  for (int col = 0; col < kN; ++col) {
    for (int state = col + 1; state < kN; ++state) {
      const double sym = 0.5 * (P_[state + kN * col] + P_[col + kN * state]);
      P_[state + kN * col] = sym;
      P_[col + kN * state] = sym;
    }
  }

  feedback_error(dx);
  project_state();
  return true;
}

void HullFusion::update_nhc(double t_sec, const double * omega_b, const double * lever_frd)
{
  if (!armed_ || !wheel_aiding_) {
    return;
  }
  if (nhc_sec_ < 0.0) {
    nhc_sec_ = t_sec;
    return;
  }
  // 周期来自 nhc_hz。留 1 ms，避免时间戳的舍入把整拍丢掉。
  if (t_sec < nhc_sec_ + nhc_period_sec_ - 1e-3) {
    return;
  }
  nhc_sec_ = t_sec;

  const int axes[2] = {1, 2};
  const double variance[2] = {kNhcVarY, kNhcVarZ};
  for (int i = 0; i < 2; ++i) {
    double c_se_t[9];
    transpose3(Cse_, c_se_t);
    double c_bs[9];
    mul3(c_se_t, Cbe_, c_bs);
    double v_s[3];
    double v_b[3];
    mul_vec3(c_se_t, ve_, v_s);
    body_velocity(c_bs, v_s, v_b);
    const int axis = axes[i];
    double row[kN];
    double prediction = 0.0;
    nhc_row(c_bs, v_b, axis, omega_b, lever_frd, row, &prediction);
    apply_scalar_update(row, prediction, 0.0, variance[i]);
  }
}

bool HullFusion::wheel_observation_holds()
{
  const double c_bs[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const double v_s[3] = {1.0, 2.0, 3.0};
  const double scale0[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double h0[kM];
  double jacobian[kM * kN];
  wheel_observation(c_bs, v_s, scale0, h0, jacobian);

  for (int row = 0; row < kM; ++row) {
    if (std::fabs(h0[row] - v_s[0]) > 1e-12) {
      return false;
    }
    for (int col = 0; col < 3; ++col) {
      if (std::fabs(jacobian[row * kN + col]) > 1e-12) {
        return false;
      }
    }
    for (int col = 6; col < 15; ++col) {
      if (std::fabs(jacobian[row * kN + col]) > 1e-12) {
        return false;
      }
    }
    for (int col = 0; col < 5; ++col) {
      const double expect = col == kWheelScale[row] ? v_s[0] : 0.0;
      if (std::fabs(jacobian[row * kN + 15 + col] - expect) > 1e-12) {
        return false;
      }
    }
  }

  const double eps = 1e-6;
  for (int axis = 0; axis < 3; ++axis) {
    double v_step[3] = {v_s[0], v_s[1], v_s[2]};
    v_step[axis] += eps;
    double h_step[kM];
    double unused[kM * kN];
    wheel_observation(c_bs, v_step, scale0, h_step, unused);
    for (int row = 0; row < kM; ++row) {
      const double numeric = (h_step[row] - h0[row]) / eps;
      if (std::fabs(numeric - jacobian[row * kN + 3 + axis]) > 1e-6) {
        return false;
      }
    }
  }

  const double scale[5] = {0.1, -0.2, 0.0, 0.3, -0.1};
  double h_scaled[kM];
  double scaled[kM * kN];
  wheel_observation(c_bs, v_s, scale, h_scaled, scaled);
  for (int row = 0; row < kM; ++row) {
    const double gain = 1.0 + scale[kWheelScale[row]];
    if (std::fabs(h_scaled[row] - gain * v_s[0]) > 1e-12) {
      return false;
    }
    for (int col = 0; col < 15; ++col) {
      if (std::fabs(scaled[row * kN + col] - gain * jacobian[row * kN + col]) > 1e-9) {
        return false;
      }
    }
    for (int col = 0; col < 5; ++col) {
      if (std::fabs(scaled[row * kN + 15 + col] - jacobian[row * kN + 15 + col]) > 1e-12) {
        return false;
      }
    }
  }
  return true;
}

bool HullFusion::nhc_observation_holds()
{
  const double c_bs[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const double v_s[3] = {1.0, 2.0, 3.0};
  const double omega0[3] = {0.0, 0.0, 0.0};
  const double lever0[3] = {0.0, 0.0, 0.0};
  double v_b[3];
  body_velocity(c_bs, v_s, v_b);
  const double eps = 1e-6;
  for (int axis = 1; axis <= 2; ++axis) {
    double row[kN];
    double prediction = 0.0;
    nhc_row(c_bs, v_b, axis, omega0, lever0, row, &prediction);
    if (std::fabs(prediction - v_b[axis]) > 1e-12) {
      return false;
    }
    for (int col = 6; col < kN; ++col) {
      if (std::fabs(row[col]) > 1e-12) {
        return false;
      }
    }
    const double expect_v[3] = {axis == 1 ? 0.0 : 0.0, axis == 1 ? 1.0 : 0.0, axis == 1 ? 0.0 : 1.0};
    for (int col = 0; col < 3; ++col) {
      if (std::fabs(row[3 + col] - expect_v[col]) > 1e-12) {
        return false;
      }
    }
    for (int col = 0; col < 3; ++col) {
      double v_step[3] = {v_s[0], v_s[1], v_s[2]};
      v_step[col] += eps;
      double v_step_b[3];
      body_velocity(c_bs, v_step, v_step_b);
      double row_step[kN];
      double pred_step = 0.0;
      nhc_row(c_bs, v_step_b, axis, omega0, lever0, row_step, &pred_step);
      if (std::fabs((pred_step - prediction) / eps - row[3 + col]) > 1e-6) {
        return false;
      }
      if (std::fabs(row[col]) > 1e-12) {
        return false;
      }
    }
  }

  const double rolled[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0};
  body_velocity(rolled, v_s, v_b);
  for (int axis = 1; axis <= 2; ++axis) {
    double row[kN];
    double prediction = 0.0;
    nhc_row(rolled, v_b, axis, omega0, lever0, row, &prediction);
    for (int col = 0; col < 3; ++col) {
      double v_step[3] = {v_s[0], v_s[1], v_s[2]};
      v_step[col] += eps;
      double v_step_b[3];
      body_velocity(rolled, v_step, v_step_b);
      double row_step[kN];
      double pred_step = 0.0;
      nhc_row(rolled, v_step_b, axis, omega0, lever0, row_step, &pred_step);
      if (std::fabs((pred_step - prediction) / eps - row[3 + col]) > 1e-6) {
        return false;
      }
      if (std::fabs(row[col]) > 1e-12) {
        return false;
      }
    }
  }

  body_velocity(c_bs, v_s, v_b);
  const double omega[3] = {0.1, -0.2, 0.3};
  const double lever[3] = {0.2, 0.0, -0.1};
  double row_lat[kN];
  double pred_lat = 0.0;
  nhc_row(c_bs, v_b, 1, omega, lever, row_lat, &pred_lat);
  const double wxr_y = omega[2] * lever[0] - omega[0] * lever[2];
  if (std::fabs(pred_lat - (v_b[1] + wxr_y)) > 1e-12) {
    return false;
  }
  for (int col = 0; col < 3; ++col) {
    if (std::fabs(row_lat[col]) > 1e-12) {
      return false;
    }
  }
  const double expect_bg[3] = {lever[2], 0.0, -lever[0]};
  for (int col = 0; col < 3; ++col) {
    if (std::fabs(row_lat[12 + col] - expect_bg[col]) > 1e-12) {
      return false;
    }
    double omega_step[3] = {omega[0], omega[1], omega[2]};
    omega_step[col] += eps;
    double row_step[kN];
    double pred_step = 0.0;
    nhc_row(c_bs, v_b, 1, omega_step, lever, row_step, &pred_step);
    const double numeric = (pred_step - pred_lat) / eps;
    const double analytic =
      col == 0 ? -lever[2] : (col == 2 ? lever[0] : 0.0);
    if (std::fabs(numeric - analytic) > 1e-6) {
      return false;
    }
  }
  for (int col = 6; col < 12; ++col) {
    if (std::fabs(row_lat[col]) > 1e-12) {
      return false;
    }
  }
  for (int col = 15; col < kN; ++col) {
    if (std::fabs(row_lat[col]) > 1e-12) {
      return false;
    }
  }
  double row_nor[kN];
  double pred_nor = 0.0;
  nhc_row(c_bs, v_b, 2, omega, lever, row_nor, &pred_nor);
  if (std::fabs(pred_nor - v_b[2]) > 1e-12) {
    return false;
  }
  for (int col = 12; col < 15; ++col) {
    if (std::fabs(row_nor[col]) > 1e-12) {
      return false;
    }
  }
  return true;
}

}  // namespace localization_v2
}  // namespace naiad
