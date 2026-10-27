#include "localization_v2/hull_wheel_ekf.hpp"

#include <cmath>

namespace naiad
{
namespace localization_v2
{
namespace
{

constexpr int kNx = HullWheelEkf::kNx;
constexpr int kNv = HullWheelEkf::kNv;
constexpr double kPi = 3.1415926535897932;
constexpr double kDeg2Rad = kPi / 180.0;
constexpr double kR = 0.1 * 0.1;

void mat_vec3(const double * a, const double * x, double * y)
{
  for (int row = 0; row < 3; ++row) {
    y[row] = a[row] * x[0] + a[row + 3] * x[1] + a[row + 6] * x[2];
  }
}

void skew(const double * a, double * s)
{
  s[0] = 0.0;
  s[1] = a[2];
  s[2] = -a[1];
  s[3] = -a[2];
  s[4] = 0.0;
  s[5] = a[0];
  s[6] = a[1];
  s[7] = -a[0];
  s[8] = 0.0;
}

void identity(double * a, int n)
{
  for (int i = 0; i < n * n; ++i) {
    a[i] = 0.0;
  }
  for (int i = 0; i < n; ++i) {
    a[i * n + i] = 1.0;
  }
}

void mul(const double * a, const double * b, double * c, int n)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += a[i * n + k] * b[k * n + j];
      }
      c[i * n + j] = sum;
    }
  }
}

void mul_abt(const double * a, const double * b, double * c, int n)
{
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += a[i * n + k] * b[j * n + k];
      }
      c[i * n + j] = sum;
    }
  }
}

void symmetrize(double * a, int n)
{
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double v = 0.5 * (a[i * n + j] + a[j * n + i]);
      a[i * n + j] = v;
      a[j * n + i] = v;
    }
  }
}

bool invert5(const double * s, double * inv)
{
  double a[5][10];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      a[i][j] = s[i * 5 + j];
      a[i][j + 5] = i == j ? 1.0 : 0.0;
    }
  }
  for (int col = 0; col < 5; ++col) {
    int pivot = col;
    double best = std::fabs(a[col][col]);
    for (int row = col + 1; row < 5; ++row) {
      const double v = std::fabs(a[row][col]);
      if (v > best) {
        best = v;
        pivot = row;
      }
    }
    if (!(best > 1.0e-12)) {
      return false;
    }
    if (pivot != col) {
      for (int j = 0; j < 10; ++j) {
        const double tmp = a[col][j];
        a[col][j] = a[pivot][j];
        a[pivot][j] = tmp;
      }
    }
    const double div = a[col][col];
    for (int j = col; j < 10; ++j) {
      a[col][j] /= div;
    }
    for (int row = 0; row < 5; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = a[row][col];
      for (int j = col; j < 10; ++j) {
        a[row][j] -= factor * a[col][j];
      }
    }
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      inv[i * 5 + j] = a[i][j + 5];
    }
  }
  return true;
}

}  // namespace

HullWheelEkf::HullWheelEkf()
{
  reset();
}

void HullWheelEkf::reset()
{
  for (int i = 0; i < kNx * kNx; ++i) {
    p_[i] = 0.0;
  }
  const double att = 5.0 * kDeg2Rad;
  const double bg = 0.1 * kDeg2Rad;
  const double var[kNx] = {
    att * att, att * att, att * att, 0.5 * 0.5, 0.5 * 0.5, 0.5 * 0.5, 1.0, 1.0, 1.0,
    0.05 * 0.05, 0.05 * 0.05, 0.05 * 0.05, bg * bg, bg * bg, bg * bg, 0.05 * 0.05, 0.05 * 0.05,
    0.05 * 0.05, 0.05 * 0.05, 0.05 * 0.05};
  for (int i = 0; i < kNx; ++i) {
    p_[i * kNx + i] = var[i];
  }
}

void HullWheelEkf::propagate(double dt, const double * c_bs, const double * f_b)
{
  if (!(dt > 0.0) || dt > 1.0) {
    return;
  }
  double cf[3];
  mat_vec3(c_bs, f_b, cf);
  double cf_skew[9];
  skew(cf, cf_skew);

  double phi[kNx * kNx];
  identity(phi, kNx);
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      phi[row * kNx + (12 + col)] += -c_bs[row + 3 * col] * dt;
      phi[(3 + row) * kNx + col] += -cf_skew[row + 3 * col] * dt;
      phi[(3 + row) * kNx + (9 + col)] += c_bs[row + 3 * col] * dt;
    }
    phi[(6 + col) * kNx + (3 + col)] += dt;
  }

  double temp[kNx * kNx];
  double next[kNx * kNx];
  mul(phi, p_, temp, kNx);
  mul_abt(temp, phi, next, kNx);

  const double arw = 0.01 * kDeg2Rad / 60.0;
  const double q_att = arw * arw;
  const double q_vel = 0.05 * 0.05;
  const double q_ba = 1.0e-8;
  const double q_bg = 1.0e-12;
  const double q_s = 1.0e-8;
  for (int i = 0; i < 3; ++i) {
    next[i * kNx + i] += q_att * dt;
    next[(3 + i) * kNx + (3 + i)] += q_vel * dt;
    next[(9 + i) * kNx + (9 + i)] += q_ba * dt;
    next[(12 + i) * kNx + (12 + i)] += q_bg * dt;
  }
  for (int i = 0; i < kNv; ++i) {
    next[(15 + i) * kNx + (15 + i)] += q_s * dt;
  }
  symmetrize(next, kNx);
  for (int i = 0; i < kNx * kNx; ++i) {
    p_[i] = next[i];
  }
}

bool HullWheelEkf::update(
  const double * c_bs, const double * v_s, const double * scale, const double * z, double * dx)
{
  const double u[3] = {c_bs[0], c_bs[1], c_bs[2]};
  const double v_bx = u[0] * v_s[0] + u[1] * v_s[1] + u[2] * v_s[2];
  const double h_phi[3] = {
    u[1] * v_s[2] - u[2] * v_s[1], u[2] * v_s[0] - u[0] * v_s[2], u[0] * v_s[1] - u[1] * v_s[0]};
  const double h_v[3] = {-u[0], -u[1], -u[2]};

  double h[kNv * kNx];
  double innov[kNv];
  for (int i = 0; i < kNv * kNx; ++i) {
    h[i] = 0.0;
  }
  for (int i = 0; i < kNv; ++i) {
    for (int k = 0; k < 3; ++k) {
      h[i * kNx + k] = h_phi[k];
      h[i * kNx + 3 + k] = h_v[k];
    }
    h[i * kNx + 15 + i] = -v_bx;
    innov[i] = z[i] - (1.0 + scale[i]) * v_bx;
  }

  double hp[kNv * kNx];
  for (int i = 0; i < kNv; ++i) {
    for (int j = 0; j < kNx; ++j) {
      double sum = 0.0;
      for (int k = 0; k < kNx; ++k) {
        sum += h[i * kNx + k] * p_[k * kNx + j];
      }
      hp[i * kNx + j] = sum;
    }
  }
  double s[kNv * kNv];
  for (int i = 0; i < kNv; ++i) {
    for (int j = 0; j < kNv; ++j) {
      double sum = i == j ? kR : 0.0;
      for (int k = 0; k < kNx; ++k) {
        sum += hp[i * kNx + k] * h[j * kNx + k];
      }
      s[i * kNv + j] = sum;
    }
  }
  double sinv[kNv * kNv];
  if (!invert5(s, sinv)) {
    return false;
  }

  double pht[kNx * kNv];
  for (int i = 0; i < kNx; ++i) {
    for (int j = 0; j < kNv; ++j) {
      pht[i * kNv + j] = hp[j * kNx + i];
    }
  }
  double k[kNx * kNv];
  for (int i = 0; i < kNx; ++i) {
    for (int j = 0; j < kNv; ++j) {
      double sum = 0.0;
      for (int m = 0; m < kNv; ++m) {
        sum += pht[i * kNv + m] * sinv[m * kNv + j];
      }
      k[i * kNv + j] = sum;
    }
  }
  for (int i = 0; i < kNx; ++i) {
    double sum = 0.0;
    for (int j = 0; j < kNv; ++j) {
      sum += k[i * kNv + j] * innov[j];
    }
    dx[i] = sum;
  }

  double kh[kNx * kNx];
  for (int i = 0; i < kNx; ++i) {
    for (int j = 0; j < kNx; ++j) {
      double sum = i == j ? 1.0 : 0.0;
      for (int m = 0; m < kNv; ++m) {
        sum -= k[i * kNv + m] * h[m * kNx + j];
      }
      kh[i * kNx + j] = sum;
    }
  }
  double temp[kNx * kNx];
  double next[kNx * kNx];
  mul(kh, p_, temp, kNx);
  mul_abt(temp, kh, next, kNx);
  for (int i = 0; i < kNx; ++i) {
    for (int j = 0; j < kNx; ++j) {
      double sum = 0.0;
      for (int m = 0; m < kNv; ++m) {
        sum += k[i * kNv + m] * kR * k[j * kNv + m];
      }
      next[i * kNx + j] += sum;
    }
  }
  symmetrize(next, kNx);
  for (int i = 0; i < kNx * kNx; ++i) {
    p_[i] = next[i];
  }
  return true;
}

}  // namespace localization_v2
}  // namespace naiad
