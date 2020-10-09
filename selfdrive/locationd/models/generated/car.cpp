
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
void err_fun(double *nom_x, double *delta_x, double *out_6670387772445835062) {
   out_6670387772445835062[0] = delta_x[0] + nom_x[0];
   out_6670387772445835062[1] = delta_x[1] + nom_x[1];
   out_6670387772445835062[2] = delta_x[2] + nom_x[2];
   out_6670387772445835062[3] = delta_x[3] + nom_x[3];
   out_6670387772445835062[4] = delta_x[4] + nom_x[4];
   out_6670387772445835062[5] = delta_x[5] + nom_x[5];
   out_6670387772445835062[6] = delta_x[6] + nom_x[6];
   out_6670387772445835062[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1723616194634524998) {
   out_1723616194634524998[0] = -nom_x[0] + true_x[0];
   out_1723616194634524998[1] = -nom_x[1] + true_x[1];
   out_1723616194634524998[2] = -nom_x[2] + true_x[2];
   out_1723616194634524998[3] = -nom_x[3] + true_x[3];
   out_1723616194634524998[4] = -nom_x[4] + true_x[4];
   out_1723616194634524998[5] = -nom_x[5] + true_x[5];
   out_1723616194634524998[6] = -nom_x[6] + true_x[6];
   out_1723616194634524998[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5963083271190511387) {
   out_5963083271190511387[0] = 1.0;
   out_5963083271190511387[1] = 0.0;
   out_5963083271190511387[2] = 0.0;
   out_5963083271190511387[3] = 0.0;
   out_5963083271190511387[4] = 0.0;
   out_5963083271190511387[5] = 0.0;
   out_5963083271190511387[6] = 0.0;
   out_5963083271190511387[7] = 0.0;
   out_5963083271190511387[8] = 0.0;
   out_5963083271190511387[9] = 1.0;
   out_5963083271190511387[10] = 0.0;
   out_5963083271190511387[11] = 0.0;
   out_5963083271190511387[12] = 0.0;
   out_5963083271190511387[13] = 0.0;
   out_5963083271190511387[14] = 0.0;
   out_5963083271190511387[15] = 0.0;
   out_5963083271190511387[16] = 0.0;
   out_5963083271190511387[17] = 0.0;
   out_5963083271190511387[18] = 1.0;
   out_5963083271190511387[19] = 0.0;
   out_5963083271190511387[20] = 0.0;
   out_5963083271190511387[21] = 0.0;
   out_5963083271190511387[22] = 0.0;
   out_5963083271190511387[23] = 0.0;
   out_5963083271190511387[24] = 0.0;
   out_5963083271190511387[25] = 0.0;
   out_5963083271190511387[26] = 0.0;
   out_5963083271190511387[27] = 1.0;
   out_5963083271190511387[28] = 0.0;
   out_5963083271190511387[29] = 0.0;
   out_5963083271190511387[30] = 0.0;
   out_5963083271190511387[31] = 0.0;
   out_5963083271190511387[32] = 0.0;
   out_5963083271190511387[33] = 0.0;
   out_5963083271190511387[34] = 0.0;
   out_5963083271190511387[35] = 0.0;
   out_5963083271190511387[36] = 1.0;
   out_5963083271190511387[37] = 0.0;
   out_5963083271190511387[38] = 0.0;
   out_5963083271190511387[39] = 0.0;
   out_5963083271190511387[40] = 0.0;
   out_5963083271190511387[41] = 0.0;
   out_5963083271190511387[42] = 0.0;
   out_5963083271190511387[43] = 0.0;
   out_5963083271190511387[44] = 0.0;
   out_5963083271190511387[45] = 1.0;
   out_5963083271190511387[46] = 0.0;
   out_5963083271190511387[47] = 0.0;
   out_5963083271190511387[48] = 0.0;
   out_5963083271190511387[49] = 0.0;
   out_5963083271190511387[50] = 0.0;
   out_5963083271190511387[51] = 0.0;
   out_5963083271190511387[52] = 0.0;
   out_5963083271190511387[53] = 0.0;
   out_5963083271190511387[54] = 1.0;
   out_5963083271190511387[55] = 0.0;
   out_5963083271190511387[56] = 0.0;
   out_5963083271190511387[57] = 0.0;
   out_5963083271190511387[58] = 0.0;
   out_5963083271190511387[59] = 0.0;
   out_5963083271190511387[60] = 0.0;
   out_5963083271190511387[61] = 0.0;
   out_5963083271190511387[62] = 0.0;
   out_5963083271190511387[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_3146926422355496422) {
   out_3146926422355496422[0] = state[0];
   out_3146926422355496422[1] = state[1];
   out_3146926422355496422[2] = state[2];
   out_3146926422355496422[3] = state[3];
   out_3146926422355496422[4] = state[4];
   out_3146926422355496422[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3146926422355496422[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3146926422355496422[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2263044337683930338) {
   out_2263044337683930338[0] = 1;
   out_2263044337683930338[1] = 0;
   out_2263044337683930338[2] = 0;
   out_2263044337683930338[3] = 0;
   out_2263044337683930338[4] = 0;
   out_2263044337683930338[5] = 0;
   out_2263044337683930338[6] = 0;
   out_2263044337683930338[7] = 0;
   out_2263044337683930338[8] = 0;
   out_2263044337683930338[9] = 1;
   out_2263044337683930338[10] = 0;
   out_2263044337683930338[11] = 0;
   out_2263044337683930338[12] = 0;
   out_2263044337683930338[13] = 0;
   out_2263044337683930338[14] = 0;
   out_2263044337683930338[15] = 0;
   out_2263044337683930338[16] = 0;
   out_2263044337683930338[17] = 0;
   out_2263044337683930338[18] = 1;
   out_2263044337683930338[19] = 0;
   out_2263044337683930338[20] = 0;
   out_2263044337683930338[21] = 0;
   out_2263044337683930338[22] = 0;
   out_2263044337683930338[23] = 0;
   out_2263044337683930338[24] = 0;
   out_2263044337683930338[25] = 0;
   out_2263044337683930338[26] = 0;
   out_2263044337683930338[27] = 1;
   out_2263044337683930338[28] = 0;
   out_2263044337683930338[29] = 0;
   out_2263044337683930338[30] = 0;
   out_2263044337683930338[31] = 0;
   out_2263044337683930338[32] = 0;
   out_2263044337683930338[33] = 0;
   out_2263044337683930338[34] = 0;
   out_2263044337683930338[35] = 0;
   out_2263044337683930338[36] = 1;
   out_2263044337683930338[37] = 0;
   out_2263044337683930338[38] = 0;
   out_2263044337683930338[39] = 0;
   out_2263044337683930338[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2263044337683930338[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2263044337683930338[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2263044337683930338[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2263044337683930338[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2263044337683930338[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2263044337683930338[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2263044337683930338[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2263044337683930338[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2263044337683930338[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2263044337683930338[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2263044337683930338[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2263044337683930338[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2263044337683930338[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2263044337683930338[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2263044337683930338[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2263044337683930338[56] = 0;
   out_2263044337683930338[57] = 0;
   out_2263044337683930338[58] = 0;
   out_2263044337683930338[59] = 0;
   out_2263044337683930338[60] = 0;
   out_2263044337683930338[61] = 0;
   out_2263044337683930338[62] = 0;
   out_2263044337683930338[63] = 1;
}
void h_25(double *state, double *unused, double *out_5106299219859223393) {
   out_5106299219859223393[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4628561448638128277) {
   out_4628561448638128277[0] = 0;
   out_4628561448638128277[1] = 0;
   out_4628561448638128277[2] = 0;
   out_4628561448638128277[3] = 0;
   out_4628561448638128277[4] = 0;
   out_4628561448638128277[5] = 0;
   out_4628561448638128277[6] = 1;
   out_4628561448638128277[7] = 0;
}
void h_24(double *state, double *unused, double *out_8231067413544154697) {
   out_8231067413544154697[0] = state[4];
   out_8231067413544154697[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8259688551304485119) {
   out_8259688551304485119[0] = 0;
   out_8259688551304485119[1] = 0;
   out_8259688551304485119[2] = 0;
   out_8259688551304485119[3] = 0;
   out_8259688551304485119[4] = 1;
   out_8259688551304485119[5] = 0;
   out_8259688551304485119[6] = 0;
   out_8259688551304485119[7] = 0;
   out_8259688551304485119[8] = 0;
   out_8259688551304485119[9] = 0;
   out_8259688551304485119[10] = 0;
   out_8259688551304485119[11] = 0;
   out_8259688551304485119[12] = 0;
   out_8259688551304485119[13] = 1;
   out_8259688551304485119[14] = 0;
   out_8259688551304485119[15] = 0;
}
void h_30(double *state, double *unused, double *out_8365820242781765213) {
   out_8365820242781765213[0] = state[4];
}
void H_30(double *state, double *unused, double *out_889755044116959247) {
   out_889755044116959247[0] = 0;
   out_889755044116959247[1] = 0;
   out_889755044116959247[2] = 0;
   out_889755044116959247[3] = 0;
   out_889755044116959247[4] = 1;
   out_889755044116959247[5] = 0;
   out_889755044116959247[6] = 0;
   out_889755044116959247[7] = 0;
}
void h_26(double *state, double *unused, double *out_2122427936646130454) {
   out_2122427936646130454[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5771970396916960493) {
   out_5771970396916960493[0] = 0;
   out_5771970396916960493[1] = 0;
   out_5771970396916960493[2] = 0;
   out_5771970396916960493[3] = 0;
   out_5771970396916960493[4] = 0;
   out_5771970396916960493[5] = 0;
   out_5771970396916960493[6] = 0;
   out_5771970396916960493[7] = 1;
}
void h_27(double *state, double *unused, double *out_4499837805307058846) {
   out_4499837805307058846[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6869984384567007633) {
   out_6869984384567007633[0] = 0;
   out_6869984384567007633[1] = 0;
   out_6869984384567007633[2] = 0;
   out_6869984384567007633[3] = 1;
   out_6869984384567007633[4] = 0;
   out_6869984384567007633[5] = 0;
   out_6869984384567007633[6] = 0;
   out_6869984384567007633[7] = 0;
}
void h_29(double *state, double *unused, double *out_5068634776543713056) {
   out_5068634776543713056[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2873556906080135365) {
   out_2873556906080135365[0] = 0;
   out_2873556906080135365[1] = 1;
   out_2873556906080135365[2] = 0;
   out_2873556906080135365[3] = 0;
   out_2873556906080135365[4] = 0;
   out_2873556906080135365[5] = 0;
   out_2873556906080135365[6] = 0;
   out_2873556906080135365[7] = 0;
}
void h_28(double *state, double *unused, double *out_914175090668094469) {
   out_914175090668094469[0] = state[5];
   out_914175090668094469[1] = state[6];
}
void H_28(double *state, double *unused, double *out_4477589415406949333) {
   out_4477589415406949333[0] = 0;
   out_4477589415406949333[1] = 0;
   out_4477589415406949333[2] = 0;
   out_4477589415406949333[3] = 0;
   out_4477589415406949333[4] = 0;
   out_4477589415406949333[5] = 1;
   out_4477589415406949333[6] = 0;
   out_4477589415406949333[7] = 0;
   out_4477589415406949333[8] = 0;
   out_4477589415406949333[9] = 0;
   out_4477589415406949333[10] = 0;
   out_4477589415406949333[11] = 0;
   out_4477589415406949333[12] = 0;
   out_4477589415406949333[13] = 0;
   out_4477589415406949333[14] = 1;
   out_4477589415406949333[15] = 0;
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
