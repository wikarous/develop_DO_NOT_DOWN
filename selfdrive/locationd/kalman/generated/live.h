/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_296826815739529133);
void inv_err_fun(double *nom_x, double *true_x, double *out_3208350509227892251);
void H_mod_fun(double *state, double *out_53166510160572014);
void f_fun(double *state, double dt, double *out_2514301809927179535);
void F_fun(double *state, double dt, double *out_4375158928109943950);
void h_3(double *state, double *unused, double *out_4109786953687873403);
void H_3(double *state, double *unused, double *out_3389870431374178295);
void h_4(double *state, double *unused, double *out_6870436501360339662);
void H_4(double *state, double *unused, double *out_4821798180251159592);
void h_9(double *state, double *unused, double *out_4833513399040659291);
void H_9(double *state, double *unused, double *out_6604613676500683868);
void h_10(double *state, double *unused, double *out_1079147239559160588);
void H_10(double *state, double *unused, double *out_6643670829633620854);
void h_12(double *state, double *unused, double *out_7165087673226052407);
void H_12(double *state, double *unused, double *out_3226820342086172692);
void h_13(double *state, double *unused, double *out_683024649669865070);
void H_13(double *state, double *unused, double *out_6941203397910262862);
void h_14(double *state, double *unused, double *out_4833513399040659291);
void H_14(double *state, double *unused, double *out_6604613676500683868);
void h_19(double *state, double *unused, double *out_4550226675792640235);
void H_19(double *state, double *unused, double *out_6502652099534772056);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);