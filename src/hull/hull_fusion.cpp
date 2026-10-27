#include "localization_v2/hull_fusion.hpp"

#include <cmath>

namespace naiad
{
namespace localization_v2
{
namespace
{

constexpr double kPi = 3.1415926535897932;
constexpr double kDeg2Rad = kPi / 180.0;
constexpr double kRad2Deg = 180.0 / kPi;
constexpr double kOmegaEarth = 7.2921151467e-5;
constexpr double kEarthA = 6378137.0;
constexpr double kEarthF = 1.0 / 298.257223563;
constexpr double kMaxDt = 1.0;

double norm3(const double * v)
{
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void cross3(const double * a, const double * b, double * c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void mat_mul(const double * a, const double * b, double * c)
{
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      c[row + 3 * col] = a[row] * b[3 * col] + a[row + 3] * b[1 + 3 * col] +
                         a[row + 6] * b[2 + 3 * col];
    }
  }
}

void mat_vec(const double * a, const double * x, double * y)
{
  for (int row = 0; row < 3; ++row) {
    y[row] = a[row] * x[0] + a[row + 3] * x[1] + a[row + 6] * x[2];
  }
}

void transpose(const double * a, double * b)
{
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      b[row + 3 * col] = a[col + 3 * row];
    }
  }
}

/// 话题轴前-上-右，算法机体系前-右-下。列主序，v_frd = C * v_fur。
void fur_to_frd(const double * fur, double * frd)
{
  frd[0] = fur[0];
  frd[1] = fur[2];
  frd[2] = -fur[1];
}

void ned_to_ecef(double lat, double lon, double * c_ne)
{
  const double sp = std::sin(lat);
  const double cp = std::cos(lat);
  const double sl = std::sin(lon);
  const double cl = std::cos(lon);
  c_ne[0] = -sp * cl;
  c_ne[1] = -sp * sl;
  c_ne[2] = cp;
  c_ne[3] = -sl;
  c_ne[4] = cl;
  c_ne[5] = 0.0;
  c_ne[6] = -cp * cl;
  c_ne[7] = -cp * sl;
  c_ne[8] = -sp;
}

void llh_to_ecef(double lat, double lon, double h, double * re)
{
  const double e2 = kEarthF * (2.0 - kEarthF);
  const double sp = std::sin(lat);
  const double cp = std::cos(lat);
  const double sl = std::sin(lon);
  const double cl = std::cos(lon);
  const double v = kEarthA / std::sqrt(1.0 - e2 * sp * sp);
  re[0] = (v + h) * cp * cl;
  re[1] = (v + h) * cp * sl;
  re[2] = (v * (1.0 - e2) + h) * sp;
}

void ecef_to_llh(const double * re, double * lat, double * lon, double * h)
{
  const double e2 = kEarthF * (2.0 - kEarthF);
  const double r2 = re[0] * re[0] + re[1] * re[1];
  double z = re[2];
  double zk = 0.0;
  double v = kEarthA;
  while (std::fabs(z - zk) >= 1e-8) {
    zk = z;
    const double sinp = z / std::sqrt(r2 + z * z);
    v = kEarthA / std::sqrt(1.0 - e2 * sinp * sinp);
    z = re[2] + v * e2 * sinp;
  }
  *lat = r2 > 1e-12 ? std::atan(z / std::sqrt(r2)) : (re[2] > 0.0 ? kPi / 2.0 : -kPi / 2.0);
  *lon = r2 > 1e-12 ? std::atan2(re[1], re[0]) : 0.0;
  *h = std::sqrt(r2 + z * z) - v;
}

double gravity_down(double lat, double h)
{
  const double e2 = kEarthF * (2.0 - kEarthF);
  const double sp = std::sin(lat);
  const double g0 =
    9.7803253359 * (1.0 + 0.001931853 * sp * sp) / std::sqrt(1.0 - e2 * sp * sp);
  return g0 * (1.0 - 2.0 * h / kEarthA);
}

/// 与 ins.cc rpy2dcm 同一组 3-2-1 角，输出机体到北东地。角的单位是弧度。
void rpy_to_cbn(double roll, double pitch, double yaw, double * c_bn)
{
  const double sp = std::sin(roll);
  const double cp = std::cos(roll);
  const double st = std::sin(pitch);
  const double ct = std::cos(pitch);
  const double sy = std::sin(yaw);
  const double cy = std::cos(yaw);
  double c_nb[9];
  c_nb[0] = ct * cy;
  c_nb[3] = ct * sy;
  c_nb[6] = -st;
  c_nb[1] = -cp * sy + sp * st * cy;
  c_nb[4] = cp * cy + sp * st * sy;
  c_nb[7] = sp * ct;
  c_nb[2] = sp * sy + cp * st * cy;
  c_nb[5] = -sp * cy + cp * st * sy;
  c_nb[8] = cp * ct;
  transpose(c_nb, c_bn);
}

void cbn_to_rpy(const double * c_bn, double * roll, double * pitch, double * yaw)
{
  double c_nb[9];
  transpose(c_bn, c_nb);
  *roll = std::atan2(c_nb[7], c_nb[8]);
  *pitch = -std::asin(std::max(-1.0, std::min(1.0, c_nb[6])));
  *yaw = std::atan2(c_nb[3], c_nb[0]);
}

void rotvec_to_quat(const double * rv, double * q)
{
  const double ang = norm3(rv);
  if (ang < 1e-8) {
    q[0] = 1.0;
    q[1] = q[2] = q[3] = 0.0;
    return;
  }
  const double s = std::sin(ang / 2.0) / ang;
  q[0] = std::cos(ang / 2.0);
  q[1] = rv[0] * s;
  q[2] = rv[1] * s;
  q[3] = rv[2] * s;
}

void quat_mul(const double * qab, const double * qca, double * qcb)
{
  qcb[0] = qab[0] * qca[0] - qab[1] * qca[1] - qab[2] * qca[2] - qab[3] * qca[3];
  qcb[1] = qab[1] * qca[0] + qab[0] * qca[1] - qab[3] * qca[2] + qab[2] * qca[3];
  qcb[2] = qab[2] * qca[0] + qab[3] * qca[1] + qab[0] * qca[2] - qab[1] * qca[3];
  qcb[3] = qab[3] * qca[0] - qab[2] * qca[1] + qab[1] * qca[2] + qab[0] * qca[3];
}

void quat_normalize(double * q)
{
  const double n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  if (n < 1e-15) {
    q[0] = 1.0;
    q[1] = q[2] = q[3] = 0.0;
    return;
  }
  q[0] /= n;
  q[1] /= n;
  q[2] /= n;
  q[3] /= n;
  if (q[0] < 0.0) {
    q[0] = -q[0];
    q[1] = -q[1];
    q[2] = -q[2];
    q[3] = -q[3];
  }
}

void quat_to_dcm(const double * q, double * c)
{
  const double a = q[0];
  const double b = q[1];
  const double d = q[2];
  const double e = q[3];
  c[0] = a * a + b * b - d * d - e * e;
  c[1] = 2.0 * (b * d + a * e);
  c[2] = 2.0 * (b * e - a * d);
  c[3] = 2.0 * (b * d - a * e);
  c[4] = a * a - b * b + d * d - e * e;
  c[5] = 2.0 * (d * e + a * b);
  c[6] = 2.0 * (b * e + a * d);
  c[7] = 2.0 * (d * e - a * b);
  c[8] = a * a - b * b - d * d + e * e;
}

void dcm_to_quat(const double * c, double * q)
{
  const double tr = c[0] + c[4] + c[8];
  const double pa = 1.0 + tr;
  const double pb = 1.0 + 2.0 * c[0] - tr;
  const double pc = 1.0 + 2.0 * c[4] - tr;
  const double pd = 1.0 + 2.0 * c[8] - tr;
  if (pa >= pb && pa >= pc && pa >= pd) {
    q[0] = 0.5 * std::sqrt(std::max(pa, 0.0));
    q[1] = (c[5] - c[7]) / (4.0 * q[0]);
    q[2] = (c[6] - c[2]) / (4.0 * q[0]);
    q[3] = (c[1] - c[3]) / (4.0 * q[0]);
  } else if (pb >= pc && pb >= pd) {
    q[1] = 0.5 * std::sqrt(std::max(pb, 0.0));
    q[2] = (c[1] + c[3]) / (4.0 * q[1]);
    q[3] = (c[6] + c[2]) / (4.0 * q[1]);
    q[0] = (c[5] - c[7]) / (4.0 * q[1]);
  } else if (pc >= pd) {
    q[2] = 0.5 * std::sqrt(std::max(pc, 0.0));
    q[3] = (c[5] + c[7]) / (4.0 * q[2]);
    q[0] = (c[6] - c[2]) / (4.0 * q[2]);
    q[1] = (c[1] + c[3]) / (4.0 * q[2]);
  } else {
    q[3] = 0.5 * std::sqrt(std::max(pd, 0.0));
    q[0] = (c[1] - c[3]) / (4.0 * q[3]);
    q[1] = (c[6] + c[2]) / (4.0 * q[3]);
    q[2] = (c[5] + c[7]) / (4.0 * q[3]);
  }
  quat_normalize(q);
}

void sculling(
  const double * w_prev, double dt_prev, const double * f_prev, const double * w, double dt,
  const double * f, double * dv, double * da)
{
  double dap[3];
  double dvp[3];
  double dak[3];
  double dvk[3];
  double dv1[3];
  double dv2[3];
  double dv3[3];
  double dv4[3];
  double tmp[3];
  for (int i = 0; i < 3; ++i) {
    dap[i] = w_prev[i] * dt_prev;
    dvp[i] = f_prev[i] * dt_prev;
    dak[i] = w[i] * dt;
    dvk[i] = f[i] * dt;
  }
  cross3(dak, dvk, dv1);
  cross3(dap, dvk, dv2);
  cross3(dvp, dak, dv3);
  const double domg = norm3(dak);
  double a1 = 0.0;
  double a2 = 0.0;
  if (std::fabs(domg) > 1e-6) {
    a1 = (1.0 - std::cos(domg)) / (domg * domg);
    a2 = (1.0 - std::sin(domg) / domg) / (domg * domg);
  } else {
    const double d2 = domg * domg;
    a1 = 0.5 - d2 / 24.0 + (d2 * d2) / 720.0;
    a2 = 1.0 / 6.0 - d2 / 120.0 + (d2 * d2) / 5040.0;
  }
  cross3(dak, dv1, dv4);
  for (int i = 0; i < 3; ++i) {
    dv[i] = a1 * dv1[i] + a2 * dv4[i] + (dv2[i] + dv3[i]) / 12.0;
  }
  cross3(dap, dak, tmp);
  for (int i = 0; i < 3; ++i) {
    da[i] = tmp[i] / 12.0;
  }
}

void ship_dcm(double psi, const double * c_ne, double * c_se)
{
  const double c = std::cos(psi);
  const double s = std::sin(psi);
  double c_sn[9] = {c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0};
  mat_mul(c_ne, c_sn, c_se);
}

}  // namespace

HullVariant HullFusion::parse(const std::string & text)
{
  if (text == "a" || text == "heading") {
    return HullVariant::HeadingOnly;
  }
  if (text == "b" || text == "ship") {
    return HullVariant::ShipNav;
  }
  return HullVariant::Auto;
}

void HullFusion::set_mode_a_origin(const ModeAOrigin & origin)
{
  origin_ = origin;
  armed_ = origin.valid;
  have_prev_ = false;
  t_sec_ = -1.0;
  wheel_sec_ = -1.0;
  nhc_sec_ = -1.0;
  ve_[0] = ve_[1] = ve_[2] = 0.0;
  for (int i = 0; i < 3; ++i) {
    ba_[i] = 0.0;
    bg_[i] = 0.0;
    f_b_[i] = 0.0;
  }
  for (int i = 0; i < 5; ++i) {
    scale_[i] = 0.0;
  }
  reset_error_covariance();
  if (!armed_) {
    nav_ = NavStatus{};
    return;
  }

  const double lat = origin.latitude_deg * kDeg2Rad;
  const double lon = origin.longitude_deg * kDeg2Rad;
  llh_to_ecef(lat, lon, origin.height_m, re_);
  double c_ne[9];
  ned_to_ecef(lat, lon, c_ne);
  double c_bn[9];
  rpy_to_cbn(
    origin.roll_deg * kDeg2Rad, origin.pitch_deg * kDeg2Rad, origin.yaw_deg * kDeg2Rad, c_bn);
  mat_mul(c_ne, c_bn, Cbe_);
  // 还没有 3 min 艏向。准船体航向先冻结为 0：艏向指北，横滚纵摇为 0。
  ship_dcm(0.0, c_ne, Cse_);
  double sp[3];
  mat_vec(Cse_, origin.p_s, sp);
  for (int i = 0; i < 3; ++i) {
    rs_[i] = re_[i] - sp[i];
  }
  project_state();
}

void HullFusion::project_state()
{
  nav_ = NavStatus{};
  nav_.backend = "hull_a";
  if (!armed_) {
    return;
  }
  nav_.filter_wired = true;
  nav_.solution_valid = true;

  double dr[3];
  for (int i = 0; i < 3; ++i) {
    dr[i] = re_[i] - rs_[i];
  }
  double cse_t[9];
  transpose(Cse_, cse_t);
  mat_vec(cse_t, dr, nav_.p_s);
  mat_vec(cse_t, ve_, nav_.v_s);

  double c_bs[9];
  mat_mul(cse_t, Cbe_, c_bs);
  double roll = 0.0;
  double pitch = 0.0;
  double yaw = 0.0;
  cbn_to_rpy(c_bs, &roll, &pitch, &yaw);
  nav_.roll_s_deg = roll * kRad2Deg;
  nav_.pitch_s_deg = pitch * kRad2Deg;
  nav_.yaw_s_deg = yaw * kRad2Deg;
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      nav_.C_bs[row * 3 + col] = c_bs[row + 3 * col];
    }
  }
  double q[4];
  dcm_to_quat(c_bs, q);
  nav_.quat_s_w = q[0];
  nav_.quat_s_x = q[1];
  nav_.quat_s_y = q[2];
  nav_.quat_s_z = q[3];

  double lat = 0.0;
  double lon = 0.0;
  double h = 0.0;
  ecef_to_llh(re_, &lat, &lon, &h);
  nav_.latitude_deg = lat * kRad2Deg;
  nav_.longitude_deg = lon * kRad2Deg;
  nav_.height_m = h;

  double c_ne[9];
  ned_to_ecef(lat, lon, c_ne);
  double c_en[9];
  transpose(c_ne, c_en);
  mat_vec(c_en, ve_, nav_.v_n);
  double c_bn[9];
  mat_mul(c_en, Cbe_, c_bn);
  cbn_to_rpy(c_bn, &roll, &pitch, &yaw);
  nav_.roll_n_deg = roll * kRad2Deg;
  nav_.pitch_n_deg = pitch * kRad2Deg;
  nav_.yaw_n_deg = yaw * kRad2Deg;
}

bool HullFusion::integrate_fog(
  double t_sec, double gyro_x, double gyro_y, double gyro_z, double acc_x, double acc_y,
  double acc_z)
{
  if (!armed_) {
    return false;
  }
  if (t_sec_ < 0.0) {
    t_sec_ = t_sec;
    return true;
  }
  const double dt = t_sec - t_sec_;
  if (!(dt > 1e-6) || dt > kMaxDt) {
    t_sec_ = t_sec;
    have_prev_ = false;
    return false;
  }

  const double fur_w[3] = {gyro_x, gyro_y, gyro_z};
  const double fur_f[3] = {acc_x, acc_y, acc_z};
  double w[3];
  double f[3];
  fur_to_frd(fur_w, w);
  fur_to_frd(fur_f, f);
  for (int i = 0; i < 3; ++i) {
    w[i] -= bg_[i];
    f[i] -= ba_[i];
  }

  double dv[3] = {0.0, 0.0, 0.0};
  double da[3] = {0.0, 0.0, 0.0};
  if (have_prev_) {
    sculling(prev_w_, prev_dt_, prev_f_, w, dt, f, dv, da);
  }

  double domg_b[3];
  for (int i = 0; i < 3; ++i) {
    domg_b[i] = w[i] * dt + da[i];
  }
  double dq_b[4];
  double dq_e[4];
  double q_old[4];
  double q_tmp[4];
  double q_new[4];
  rotvec_to_quat(domg_b, dq_b);
  const double domg_e[3] = {0.0, 0.0, -kOmegaEarth * dt};
  rotvec_to_quat(domg_e, dq_e);
  dcm_to_quat(Cbe_, q_old);
  quat_mul(q_old, dq_b, q_tmp);
  quat_mul(dq_e, q_tmp, q_new);
  quat_normalize(q_new);

  double c_old[9];
  for (int i = 0; i < 9; ++i) {
    c_old[i] = Cbe_[i];
  }
  quat_to_dcm(q_new, Cbe_);
  double d_ce[9];
  quat_to_dcm(dq_e, d_ce);

  double dvb[3];
  for (int i = 0; i < 3; ++i) {
    dvb[i] = f[i] * dt + dv[i];
  }
  double mid[3];
  double dvf[3];
  mat_vec(c_old, dvb, mid);
  mat_vec(d_ce, mid, dvf);

  double lat = 0.0;
  double lon = 0.0;
  double h = 0.0;
  ecef_to_llh(re_, &lat, &lon, &h);
  double c_ne[9];
  ned_to_ecef(lat, lon, c_ne);
  const double gn[3] = {0.0, 0.0, gravity_down(lat, h)};
  double ge[3];
  mat_vec(c_ne, gn, ge);
  const double omega_e[3] = {0.0, 0.0, kOmegaEarth};
  double cori[3];
  cross3(omega_e, ve_, cori);

  double ve_old[3];
  for (int i = 0; i < 3; ++i) {
    ve_old[i] = ve_[i];
    ve_[i] += dvf[i] + (ge[i] - 2.0 * cori[i]) * dt;
    re_[i] += 0.5 * (ve_old[i] + ve_[i]) * dt;
  }

  for (int i = 0; i < 3; ++i) {
    prev_w_[i] = w[i];
    prev_f_[i] = f[i];
  }
  prev_dt_ = dt;
  have_prev_ = true;
  t_sec_ = t_sec;
  for (int i = 0; i < 3; ++i) {
    f_b_[i] = f[i];
  }
  if (wheel_aiding_) {
    propagate_error(dt, f_b_);
  }
  project_state();
  if (wheel_aiding_) {
    update_nhc(t_sec, w);
  }
  return true;
}

bool HullFusion::level_specific_force_holds()
{
  ModeAOrigin origin;
  origin.valid = true;
  origin.latitude_deg = 30.0;
  origin.longitude_deg = 114.0;
  const double g = gravity_down(origin.latitude_deg * kDeg2Rad, 0.0);
  HullFusion fusion;
  fusion.set_mode_a_origin(origin);
  constexpr int kSteps = 200;
  constexpr double kDt = 0.01;
  for (int i = 0; i <= kSteps; ++i) {
    if (!fusion.integrate_fog(i * kDt, 0.0, 0.0, 0.0, 0.0, g, 0.0)) {
      return false;
    }
  }
  const NavStatus nav = fusion.step(false);
  const double speed = norm3(nav.v_n);
  const double drift = std::sqrt(
    nav.p_s[0] * nav.p_s[0] + nav.p_s[1] * nav.p_s[1] + nav.p_s[2] * nav.p_s[2]);
  return speed < 0.05 && drift < 0.05 && std::fabs(nav.height_m) < 0.05;
}

NavStatus HullFusion::step(bool xcar_fresh) const
{
  const bool use_ship = variant_ == HullVariant::ShipNav ||
                        (variant_ == HullVariant::Auto && xcar_fresh);
  if (use_ship || !armed_) {
    NavStatus status;
    status.backend = use_ship ? "hull_b" : "hull_a";
    return status;
  }
  return nav_;
}

}  // namespace localization_v2
}  // namespace naiad
