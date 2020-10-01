/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6873866917236686546);
void inv_err_fun(double *nom_x, double *true_x, double *out_3308672848029761950);
void H_mod_fun(double *state, double *out_8751023711912772613);
void f_fun(double *state, double dt, double *out_2421591416499194210);
void F_fun(double *state, double dt, double *out_8475765992517630234);
void h_25(double *state, double *unused, double *out_5205267724522008689);
void H_25(double *state, double *unused, double *out_6669771596278746555);
void h_24(double *state, double *unused, double *out_2530730875656060729);
void H_24(double *state, double *unused, double *out_5182526297256752593);
void h_30(double *state, double *unused, double *out_4519327413977263149);
void H_30(double *state, double *unused, double *out_4067313812378998465);
void h_26(double *state, double *unused, double *out_8225367038014836794);
void H_26(double *state, double *unused, double *out_5407219630005040061);
void h_27(double *state, double *unused, double *out_6019300973347613614);
void H_27(double *state, double *unused, double *out_3270708381515193441);
void h_29(double *state, double *unused, double *out_849401770996740336);
void H_29(double *state, double *unused, double *out_2844204770730783509);
void h_28(double *state, double *unused, double *out_6878993397801909301);
void H_28(double *state, double *unused, double *out_315058093637771077);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
