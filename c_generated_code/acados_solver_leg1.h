/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#ifndef ACADOS_SOLVER_leg1_H_
#define ACADOS_SOLVER_leg1_H_

#include "acados/utils/types.h"

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#define LEG1_NX     19
#define LEG1_NZ     0
#define LEG1_NU     15
#define LEG1_NP     48
#define LEG1_NBX    0
#define LEG1_NBX0   19
#define LEG1_NBU    0
#define LEG1_NSBX   0
#define LEG1_NSBU   0
#define LEG1_NSH    3
#define LEG1_NSH0   0
#define LEG1_NSG    0
#define LEG1_NSPHI  0
#define LEG1_NSHN   0
#define LEG1_NSGN   0
#define LEG1_NSPHIN 0
#define LEG1_NSPHI0 0
#define LEG1_NSBXN  0
#define LEG1_NS     3
#define LEG1_NS0    0
#define LEG1_NSN    0
#define LEG1_NG     0
#define LEG1_NBXN   0
#define LEG1_NGN    0
#define LEG1_NY0    62
#define LEG1_NY     62
#define LEG1_NYN    21
#define LEG1_N      15
#define LEG1_NH     3
#define LEG1_NHN    0
#define LEG1_NH0    0
#define LEG1_NPHI0  0
#define LEG1_NPHI   0
#define LEG1_NPHIN  0
#define LEG1_NR     0

#ifdef __cplusplus
extern "C" {
#endif


// ** capsule for solver data **
typedef struct leg1_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_out *sens_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan_t *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

    // number of expected runtime parameters
    unsigned int nlp_np;

    /* external functions */
    // dynamics

    external_function_external_param_casadi *discr_dyn_phi_fun;
    external_function_external_param_casadi *discr_dyn_phi_fun_jac_ut_xt;




    // cost

    external_function_external_param_casadi *cost_y_fun;
    external_function_external_param_casadi *cost_y_fun_jac_ut_xt;
    external_function_external_param_casadi *cost_y_hess;



    external_function_external_param_casadi cost_y_0_fun;
    external_function_external_param_casadi cost_y_0_fun_jac_ut_xt;
    external_function_external_param_casadi cost_y_0_hess;



    external_function_external_param_casadi cost_y_e_fun;
    external_function_external_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_external_param_casadi cost_y_e_hess;


    // constraints
    external_function_external_param_casadi *nl_constr_h_fun_jac;
    external_function_external_param_casadi *nl_constr_h_fun;







} leg1_solver_capsule;

ACADOS_SYMBOL_EXPORT leg1_solver_capsule * leg1_acados_create_capsule(void);
ACADOS_SYMBOL_EXPORT int leg1_acados_free_capsule(leg1_solver_capsule *capsule);

ACADOS_SYMBOL_EXPORT int leg1_acados_create(leg1_solver_capsule * capsule);

ACADOS_SYMBOL_EXPORT int leg1_acados_reset(leg1_solver_capsule* capsule, int reset_qp_solver_mem);

/**
 * Generic version of leg1_acados_create which allows to use a different number of shooting intervals than
 * the number used for code generation. If new_time_steps=NULL and n_time_steps matches the number used for code
 * generation, the time-steps from code generation is used.
 */
ACADOS_SYMBOL_EXPORT int leg1_acados_create_with_discretization(leg1_solver_capsule * capsule, int n_time_steps, double* new_time_steps);
/**
 * Update the time step vector. Number N must be identical to the currently set number of shooting nodes in the
 * nlp_solver_plan. Returns 0 if no error occurred and a otherwise a value other than 0.
 */
ACADOS_SYMBOL_EXPORT int leg1_acados_update_time_steps(leg1_solver_capsule * capsule, int N, double* new_time_steps);
/**
 * This function is used for updating an already initialized solver with a different number of qp_cond_N.
 */
ACADOS_SYMBOL_EXPORT int leg1_acados_update_qp_solver_cond_N(leg1_solver_capsule * capsule, int qp_solver_cond_N);
ACADOS_SYMBOL_EXPORT int leg1_acados_update_params(leg1_solver_capsule * capsule, int stage, double *value, int np);
ACADOS_SYMBOL_EXPORT int leg1_acados_update_params_sparse(leg1_solver_capsule * capsule, int stage, int *idx, double *p, int n_update);

ACADOS_SYMBOL_EXPORT int leg1_acados_solve(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void leg1_acados_batch_solve(leg1_solver_capsule ** capsules, int N_batch);
ACADOS_SYMBOL_EXPORT int leg1_acados_free(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void leg1_acados_print_stats(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT int leg1_acados_custom_update(leg1_solver_capsule* capsule, double* data, int data_len);


ACADOS_SYMBOL_EXPORT ocp_nlp_in *leg1_acados_get_nlp_in(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *leg1_acados_get_nlp_out(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_out *leg1_acados_get_sens_out(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_solver *leg1_acados_get_nlp_solver(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_config *leg1_acados_get_nlp_config(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT void *leg1_acados_get_nlp_opts(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_dims *leg1_acados_get_nlp_dims(leg1_solver_capsule * capsule);
ACADOS_SYMBOL_EXPORT ocp_nlp_plan_t *leg1_acados_get_nlp_plan(leg1_solver_capsule * capsule);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_leg1_H_
