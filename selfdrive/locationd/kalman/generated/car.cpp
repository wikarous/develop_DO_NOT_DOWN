
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6873866917236686546) {
   out_6873866917236686546[0] = delta_x[0] + nom_x[0];
   out_6873866917236686546[1] = delta_x[1] + nom_x[1];
   out_6873866917236686546[2] = delta_x[2] + nom_x[2];
   out_6873866917236686546[3] = delta_x[3] + nom_x[3];
   out_6873866917236686546[4] = delta_x[4] + nom_x[4];
   out_6873866917236686546[5] = delta_x[5] + nom_x[5];
   out_6873866917236686546[6] = delta_x[6] + nom_x[6];
   out_6873866917236686546[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3308672848029761950) {
   out_3308672848029761950[0] = -nom_x[0] + true_x[0];
   out_3308672848029761950[1] = -nom_x[1] + true_x[1];
   out_3308672848029761950[2] = -nom_x[2] + true_x[2];
   out_3308672848029761950[3] = -nom_x[3] + true_x[3];
   out_3308672848029761950[4] = -nom_x[4] + true_x[4];
   out_3308672848029761950[5] = -nom_x[5] + true_x[5];
   out_3308672848029761950[6] = -nom_x[6] + true_x[6];
   out_3308672848029761950[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8751023711912772613) {
   out_8751023711912772613[0] = 1.0;
   out_8751023711912772613[1] = 0.0;
   out_8751023711912772613[2] = 0.0;
   out_8751023711912772613[3] = 0.0;
   out_8751023711912772613[4] = 0.0;
   out_8751023711912772613[5] = 0.0;
   out_8751023711912772613[6] = 0.0;
   out_8751023711912772613[7] = 0.0;
   out_8751023711912772613[8] = 0.0;
   out_8751023711912772613[9] = 1.0;
   out_8751023711912772613[10] = 0.0;
   out_8751023711912772613[11] = 0.0;
   out_8751023711912772613[12] = 0.0;
   out_8751023711912772613[13] = 0.0;
   out_8751023711912772613[14] = 0.0;
   out_8751023711912772613[15] = 0.0;
   out_8751023711912772613[16] = 0.0;
   out_8751023711912772613[17] = 0.0;
   out_8751023711912772613[18] = 1.0;
   out_8751023711912772613[19] = 0.0;
   out_8751023711912772613[20] = 0.0;
   out_8751023711912772613[21] = 0.0;
   out_8751023711912772613[22] = 0.0;
   out_8751023711912772613[23] = 0.0;
   out_8751023711912772613[24] = 0.0;
   out_8751023711912772613[25] = 0.0;
   out_8751023711912772613[26] = 0.0;
   out_8751023711912772613[27] = 1.0;
   out_8751023711912772613[28] = 0.0;
   out_8751023711912772613[29] = 0.0;
   out_8751023711912772613[30] = 0.0;
   out_8751023711912772613[31] = 0.0;
   out_8751023711912772613[32] = 0.0;
   out_8751023711912772613[33] = 0.0;
   out_8751023711912772613[34] = 0.0;
   out_8751023711912772613[35] = 0.0;
   out_8751023711912772613[36] = 1.0;
   out_8751023711912772613[37] = 0.0;
   out_8751023711912772613[38] = 0.0;
   out_8751023711912772613[39] = 0.0;
   out_8751023711912772613[40] = 0.0;
   out_8751023711912772613[41] = 0.0;
   out_8751023711912772613[42] = 0.0;
   out_8751023711912772613[43] = 0.0;
   out_8751023711912772613[44] = 0.0;
   out_8751023711912772613[45] = 1.0;
   out_8751023711912772613[46] = 0.0;
   out_8751023711912772613[47] = 0.0;
   out_8751023711912772613[48] = 0.0;
   out_8751023711912772613[49] = 0.0;
   out_8751023711912772613[50] = 0.0;
   out_8751023711912772613[51] = 0.0;
   out_8751023711912772613[52] = 0.0;
   out_8751023711912772613[53] = 0.0;
   out_8751023711912772613[54] = 1.0;
   out_8751023711912772613[55] = 0.0;
   out_8751023711912772613[56] = 0.0;
   out_8751023711912772613[57] = 0.0;
   out_8751023711912772613[58] = 0.0;
   out_8751023711912772613[59] = 0.0;
   out_8751023711912772613[60] = 0.0;
   out_8751023711912772613[61] = 0.0;
   out_8751023711912772613[62] = 0.0;
   out_8751023711912772613[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2421591416499194210) {
   out_2421591416499194210[0] = state[0];
   out_2421591416499194210[1] = state[1];
   out_2421591416499194210[2] = state[2];
   out_2421591416499194210[3] = state[3];
   out_2421591416499194210[4] = state[4];
   out_2421591416499194210[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2421591416499194210[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2421591416499194210[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8475765992517630234) {
   out_8475765992517630234[0] = 1;
   out_8475765992517630234[1] = 0;
   out_8475765992517630234[2] = 0;
   out_8475765992517630234[3] = 0;
   out_8475765992517630234[4] = 0;
   out_8475765992517630234[5] = 0;
   out_8475765992517630234[6] = 0;
   out_8475765992517630234[7] = 0;
   out_8475765992517630234[8] = 0;
   out_8475765992517630234[9] = 1;
   out_8475765992517630234[10] = 0;
   out_8475765992517630234[11] = 0;
   out_8475765992517630234[12] = 0;
   out_8475765992517630234[13] = 0;
   out_8475765992517630234[14] = 0;
   out_8475765992517630234[15] = 0;
   out_8475765992517630234[16] = 0;
   out_8475765992517630234[17] = 0;
   out_8475765992517630234[18] = 1;
   out_8475765992517630234[19] = 0;
   out_8475765992517630234[20] = 0;
   out_8475765992517630234[21] = 0;
   out_8475765992517630234[22] = 0;
   out_8475765992517630234[23] = 0;
   out_8475765992517630234[24] = 0;
   out_8475765992517630234[25] = 0;
   out_8475765992517630234[26] = 0;
   out_8475765992517630234[27] = 1;
   out_8475765992517630234[28] = 0;
   out_8475765992517630234[29] = 0;
   out_8475765992517630234[30] = 0;
   out_8475765992517630234[31] = 0;
   out_8475765992517630234[32] = 0;
   out_8475765992517630234[33] = 0;
   out_8475765992517630234[34] = 0;
   out_8475765992517630234[35] = 0;
   out_8475765992517630234[36] = 1;
   out_8475765992517630234[37] = 0;
   out_8475765992517630234[38] = 0;
   out_8475765992517630234[39] = 0;
   out_8475765992517630234[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8475765992517630234[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8475765992517630234[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8475765992517630234[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8475765992517630234[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8475765992517630234[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8475765992517630234[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8475765992517630234[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8475765992517630234[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8475765992517630234[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8475765992517630234[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8475765992517630234[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8475765992517630234[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8475765992517630234[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8475765992517630234[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8475765992517630234[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8475765992517630234[56] = 0;
   out_8475765992517630234[57] = 0;
   out_8475765992517630234[58] = 0;
   out_8475765992517630234[59] = 0;
   out_8475765992517630234[60] = 0;
   out_8475765992517630234[61] = 0;
   out_8475765992517630234[62] = 0;
   out_8475765992517630234[63] = 1;
}
void h_25(double *state, double *unused, double *out_5205267724522008689) {
   out_5205267724522008689[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6669771596278746555) {
   out_6669771596278746555[0] = 0;
   out_6669771596278746555[1] = 0;
   out_6669771596278746555[2] = 0;
   out_6669771596278746555[3] = 0;
   out_6669771596278746555[4] = 0;
   out_6669771596278746555[5] = 0;
   out_6669771596278746555[6] = 1;
   out_6669771596278746555[7] = 0;
}
void h_24(double *state, double *unused, double *out_2530730875656060729) {
   out_2530730875656060729[0] = state[4];
   out_2530730875656060729[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5182526297256752593) {
   out_5182526297256752593[0] = 0;
   out_5182526297256752593[1] = 0;
   out_5182526297256752593[2] = 0;
   out_5182526297256752593[3] = 0;
   out_5182526297256752593[4] = 1;
   out_5182526297256752593[5] = 0;
   out_5182526297256752593[6] = 0;
   out_5182526297256752593[7] = 0;
   out_5182526297256752593[8] = 0;
   out_5182526297256752593[9] = 0;
   out_5182526297256752593[10] = 0;
   out_5182526297256752593[11] = 0;
   out_5182526297256752593[12] = 0;
   out_5182526297256752593[13] = 1;
   out_5182526297256752593[14] = 0;
   out_5182526297256752593[15] = 0;
}
void h_30(double *state, double *unused, double *out_4519327413977263149) {
   out_4519327413977263149[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4067313812378998465) {
   out_4067313812378998465[0] = 0;
   out_4067313812378998465[1] = 0;
   out_4067313812378998465[2] = 0;
   out_4067313812378998465[3] = 0;
   out_4067313812378998465[4] = 1;
   out_4067313812378998465[5] = 0;
   out_4067313812378998465[6] = 0;
   out_4067313812378998465[7] = 0;
}
void h_26(double *state, double *unused, double *out_8225367038014836794) {
   out_8225367038014836794[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5407219630005040061) {
   out_5407219630005040061[0] = 0;
   out_5407219630005040061[1] = 0;
   out_5407219630005040061[2] = 0;
   out_5407219630005040061[3] = 0;
   out_5407219630005040061[4] = 0;
   out_5407219630005040061[5] = 0;
   out_5407219630005040061[6] = 0;
   out_5407219630005040061[7] = 1;
}
void h_27(double *state, double *unused, double *out_6019300973347613614) {
   out_6019300973347613614[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3270708381515193441) {
   out_3270708381515193441[0] = 0;
   out_3270708381515193441[1] = 0;
   out_3270708381515193441[2] = 0;
   out_3270708381515193441[3] = 1;
   out_3270708381515193441[4] = 0;
   out_3270708381515193441[5] = 0;
   out_3270708381515193441[6] = 0;
   out_3270708381515193441[7] = 0;
}
void h_29(double *state, double *unused, double *out_849401770996740336) {
   out_849401770996740336[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2844204770730783509) {
   out_2844204770730783509[0] = 0;
   out_2844204770730783509[1] = 1;
   out_2844204770730783509[2] = 0;
   out_2844204770730783509[3] = 0;
   out_2844204770730783509[4] = 0;
   out_2844204770730783509[5] = 0;
   out_2844204770730783509[6] = 0;
   out_2844204770730783509[7] = 0;
}
void h_28(double *state, double *unused, double *out_6878993397801909301) {
   out_6878993397801909301[0] = state[5];
   out_6878993397801909301[1] = state[6];
}
void H_28(double *state, double *unused, double *out_315058093637771077) {
   out_315058093637771077[0] = 0;
   out_315058093637771077[1] = 0;
   out_315058093637771077[2] = 0;
   out_315058093637771077[3] = 0;
   out_315058093637771077[4] = 0;
   out_315058093637771077[5] = 1;
   out_315058093637771077[6] = 0;
   out_315058093637771077[7] = 0;
   out_315058093637771077[8] = 0;
   out_315058093637771077[9] = 0;
   out_315058093637771077[10] = 0;
   out_315058093637771077[11] = 0;
   out_315058093637771077[12] = 0;
   out_315058093637771077[13] = 0;
   out_315058093637771077[14] = 1;
   out_315058093637771077[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
