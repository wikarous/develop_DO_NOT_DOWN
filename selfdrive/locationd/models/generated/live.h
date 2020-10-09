/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3278828955549405305);
void inv_err_fun(double *nom_x, double *true_x, double *out_3270294043112162685);
void H_mod_fun(double *state, double *out_1335644347630386101);
void f_fun(double *state, double dt, double *out_7863232366499697514);
void F_fun(double *state, double dt, double *out_2957341874054674020);
void h_3(double *state, double *unused, double *out_1036843898292969884);
void H_3(double *state, double *unused, double *out_2841554328525422957);
void h_4(double *state, double *unused, double *out_7800680430184608793);
void H_4(double *state, double *unused, double *out_728017061541415058);
void h_9(double *state, double *unused, double *out_5597121128106113287);
void H_9(double *state, double *unused, double *out_8496954625026249569);
void h_10(double *state, double *unused, double *out_4973311885430861647);
void H_10(double *state, double *unused, double *out_3185535861216182169);
void h_12(double *state, double *unused, double *out_5318269436441727973);
void H_12(double *state, double *unused, double *out_2843317785675348239);
void h_31(double *state, double *unused, double *out_4302711670370225146);
void H_31(double *state, double *unused, double *out_7028772316665312343);
void h_32(double *state, double *unused, double *out_3885081002420569401);
void H_32(double *state, double *unused, double *out_7783634600653251257);
void h_13(double *state, double *unused, double *out_9151854354421688967);
void H_13(double *state, double *unused, double *out_4209812053670735666);
void h_14(double *state, double *unused, double *out_5597121128106113287);
void H_14(double *state, double *unused, double *out_8496954625026249569);
void h_19(double *state, double *unused, double *out_6661139289296815191);
void H_19(double *state, double *unused, double *out_5137528230682428171);
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
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);