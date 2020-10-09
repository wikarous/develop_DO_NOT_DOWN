/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6670387772445835062);
void inv_err_fun(double *nom_x, double *true_x, double *out_1723616194634524998);
void H_mod_fun(double *state, double *out_5963083271190511387);
void f_fun(double *state, double dt, double *out_3146926422355496422);
void F_fun(double *state, double dt, double *out_2263044337683930338);
void h_25(double *state, double *unused, double *out_5106299219859223393);
void H_25(double *state, double *unused, double *out_4628561448638128277);
void h_24(double *state, double *unused, double *out_8231067413544154697);
void H_24(double *state, double *unused, double *out_8259688551304485119);
void h_30(double *state, double *unused, double *out_8365820242781765213);
void H_30(double *state, double *unused, double *out_889755044116959247);
void h_26(double *state, double *unused, double *out_2122427936646130454);
void H_26(double *state, double *unused, double *out_5771970396916960493);
void h_27(double *state, double *unused, double *out_4499837805307058846);
void H_27(double *state, double *unused, double *out_6869984384567007633);
void h_29(double *state, double *unused, double *out_5068634776543713056);
void H_29(double *state, double *unused, double *out_2873556906080135365);
void h_28(double *state, double *unused, double *out_914175090668094469);
void H_28(double *state, double *unused, double *out_4477589415406949333);
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
