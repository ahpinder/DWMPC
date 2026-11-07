#include "controllers/dwmpc/acadosInterface.hpp"

#include "acados_solver_front.h"
#include "acados_solver_back.h"
#include "acados_solver_leg0.h"
#include "acados_solver_leg1.h"
#include "acados_solver_leg2.h"
#include "acados_solver_leg3.h"

#include <algorithm>
#include <cerrno>
#include <system_error>
#include <time.h>

#define SOLVER_OPS(NAME)                                                                                   \
    {#NAME, SolverOps{                                                                                     \
        []() { return static_cast<void *>(NAME##_acados_create_capsule()); },                              \
        [](void *c) { return NAME##_acados_create(static_cast<NAME##_solver_capsule *>(c)); },             \
        [](void *c) { return NAME##_acados_free(static_cast<NAME##_solver_capsule *>(c)); },              \
        [](void *c) { return NAME##_acados_free_capsule(static_cast<NAME##_solver_capsule *>(c)); },      \
        [](void *c) { return NAME##_acados_solve(static_cast<NAME##_solver_capsule *>(c)); },             \
        [](void *c, int k, double *p, int np) {                                                            \
            return NAME##_acados_update_params(static_cast<NAME##_solver_capsule *>(c), k, p, np); },      \
        [](void *c) { return NAME##_acados_get_nlp_config(static_cast<NAME##_solver_capsule *>(c)); },    \
        [](void *c) { return NAME##_acados_get_nlp_dims(static_cast<NAME##_solver_capsule *>(c)); },      \
        [](void *c) { return NAME##_acados_get_nlp_in(static_cast<NAME##_solver_capsule *>(c)); },        \
        [](void *c) { return NAME##_acados_get_nlp_out(static_cast<NAME##_solver_capsule *>(c)); },       \
        [](void *c) { return NAME##_acados_get_nlp_solver(static_cast<NAME##_solver_capsule *>(c)); },    \
        [](void *c) { return NAME##_acados_get_nlp_opts(static_cast<NAME##_solver_capsule *>(c)); }}}

// generated solvers available by name (see generate_ocp.py)
static const std::map<std::string,SolverOps> &solverRegistry()
{
    static const std::map<std::string,SolverOps> registry{
        SOLVER_OPS(front), SOLVER_OPS(back),
        SOLVER_OPS(leg0), SOLVER_OPS(leg1), SOLVER_OPS(leg2), SOLVER_OPS(leg3)};
    return registry;
}

acadosInterface::acadosInterface()
{}
acadosInterface::~acadosInterface()
{
    clear();
}
void acadosInterface::clear()
{
    run = false;
    for (auto &entry : problems_)
    {
        entry.second->start.store(true);
        entry.second->start.notify_one();
    }
    for (auto &entry : problems_)
    {
        Problem &p = *entry.second;
        if (p.worker.joinable())
            p.worker.join();
        p.ops->free(p.capsule);
        p.ops->free_capsule(p.capsule);
    }
    problems_.clear();
}
void acadosInterface::init(const int N, const std::vector<std::string> &problems)
{
    clear();
    run = true;
    N_ = N;
    for (const auto &name : problems)
    {
        if (name == "wb")
            continue;
        auto it = solverRegistry().find(name);
        if (it == solverRegistry().end())
            throw std::runtime_error("acadosInterface: no generated solver named '" + name + "'");

        auto p = std::make_unique<Problem>();
        p->ops = &it->second;
        p->capsule = p->ops->create_capsule();
        if (p->ops->create(p->capsule) != 0)
            throw std::runtime_error("acadosInterface: could not create the solver '" + name + "'");
        p->config = p->ops->get_config(p->capsule);
        p->dims = p->ops->get_dims(p->capsule);
        p->in = p->ops->get_in(p->capsule);
        p->out = p->ops->get_out(p->capsule);
        p->solver = p->ops->get_solver(p->capsule);
        p->opts = p->ops->get_opts(p->capsule);
        p->nx = *p->dims->nx;
        p->nu = *p->dims->nu;
        Problem &ref = *p;
        problems_[name] = std::move(p);
        ref.worker = std::thread(&acadosInterface::work, this, std::ref(ref));
    }
}
/**
 * Set the initial condition for the optimization problem.
 *
 * @param initial_condition The initial condition vector.
 * @param problem The name of the problem (e.g. "front" or "leg0").
 */
void acadosInterface::setInitialCondition(std::vector<double> &initial_condition,const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    std::vector<int> idxbx0(p.nx);

    for (int counter{0};counter<p.nx;counter++)
    {
        idxbx0[counter] = counter;
    }

    // Set the initial condition constraints for the optimization problem
    ocp_nlp_constraints_model_set(p.config, p.dims, p.in, 0, "idxbx", idxbx0.data());
    ocp_nlp_constraints_model_set(p.config, p.dims, p.in, 0, "lbx", initial_condition.data());
    ocp_nlp_constraints_model_set(p.config, p.dims, p.in, 0, "ubx", initial_condition.data());
}
/**
 * Set the parameter values for the optimization problem.
 *
 * @param parameter The parameter values for each time step size [np X N].
 * @param np The number of parameters per stage.
 * @param problem The name of the problem.
 */
void acadosInterface::setParameter(std::vector<std::vector<double>> &parameter, const int np, const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    for (int k = 0; k <= N_; k++)
    {
        p.ops->update_params(p.capsule, k, parameter[k].data(), np);
    }
}
void acadosInterface::setConstraints( std::vector<std::vector<double>> &constraint, const int nlh0, const int nlh, const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    for (int k = 0; k <= N_; k++)
    { 
        ocp_nlp_constraints_model_set(p.config, p.dims, p.in, k, "lh", constraint[k].data());
    }
}
/**
 * Set the weight matrix for the cost function of the optimization problem.
 *
 * @param weight The weight matrix for each time step size [nx+nu X nx+nu] X N [nx X nx] X 1 for the last stage.
 * @param problem The name of the problem.
 */
void acadosInterface::setWeight(std::vector<std::vector<double>> &weight,const int &n_w0,const int &n_w, const int &n_we,const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    for(int k{0}; k < N_+1;k++)
    {
        ocp_nlp_cost_model_set(p.config, p.dims, p.in, k, "W", weight[k].data());
    }
}

/**
 * Set the reference values for the optimization problem.
 *
 * @param reference The reference values for each time step size [nx_ref+nu_ref X N] and [nx_ref X 1] for the last stage.
 * @param problem The name of the problem.
 */
void acadosInterface::setReference(std::vector<std::vector<double>> &reference, const int n_ref0, const int n_ref, const int n_ref_e, const std::string &problem)
{   
    Problem &p = *problems_.at(problem);
    for (int k = 0; k <= N_; k++)
    {
        ocp_nlp_cost_model_set(p.config, p.dims, p.in, k, "yref", reference[k].data());
    }
}

/**
 * Run one RTI phase of every problem, each one in its own thread, and wait for all of them.
 */
void acadosInterface::runAll(int rti_phase)
{
    const auto start = timing_enabled_ ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    for (auto &entry : problems_)
    {
        Problem &p = *entry.second;
        ocp_nlp_solver_opts_set(p.config, p.opts, "rti_phase", &rti_phase);
    }
    for (auto &entry : problems_)
    {
        entry.second->start.store(true);
        entry.second->start.notify_one();
    }
    waitWorkers();
    if (timing_enabled_)
    {
        last_timing_ = {};
        last_timing_.wall_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        for (const auto &entry : problems_)
        {
            const Problem &p = *entry.second;
            if (p.clock_error)
                throw std::system_error(p.clock_error, std::generic_category(), "worker CPU clock");
            last_timing_.worker_wall_sum_ms += p.wall_ms;
            last_timing_.worker_wall_max_ms = std::max(last_timing_.worker_wall_max_ms, p.wall_ms);
            last_timing_.worker_cpu_sum_ms += p.cpu_ms;
        }
    }
}
/**
 * Solve (feedback phase) all the problems.
 */
void acadosInterface::solve()
{
    runAll(2);
}
/**
 * Prepare (preparation phase) all the problems.
 */
void acadosInterface::prepare()
{
    runAll(1);
}
void acadosInterface::work(Problem &p)
{
    while(run)
    {
        p.start.wait(false);
        if(!run)
            break;
        if (timing_enabled_)
        {
            timespec cpu_start{}, cpu_end{};
            p.clock_error = 0;
            if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpu_start) != 0)
                p.clock_error = errno;
            const auto start = std::chrono::steady_clock::now();
            p.status = p.ops->solve(p.capsule);
            const auto end = std::chrono::steady_clock::now();
            if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &cpu_end) != 0 && !p.clock_error)
                p.clock_error = errno;
            p.wall_ms = std::chrono::duration<double, std::milli>(end - start).count();
            p.cpu_ms = (cpu_end.tv_sec - cpu_start.tv_sec) * 1000.0 +
                       (cpu_end.tv_nsec - cpu_start.tv_nsec) / 1e6;
        }
        else
            p.status = p.ops->solve(p.capsule);
        p.start.store(false);
        p.start.notify_all(); // wake the thread sleeping in waitWorkers
    }
}

/**
 * Update the prediction vectors for states and inputs.
 *
 * This function retrieves the predicted states and inputs from the solver output
 * and updates the provided vectors accordingly.
 *
 * @param x The vector of predicted states [nx X N+1].
 * @param u The vector of predicted inputs [nu X N].
 * @param problem The name of the problem.
 */
void acadosInterface::UpdatePrediction(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &u, const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    std::vector<double> xtraj(p.nx);
    std::vector<double> utraj(p.nu);

    for (int k{0}; k < N_; k++)
    {
        ocp_nlp_out_get(p.config, p.dims, p.out, k, "x", xtraj.data());
        x.push_back(xtraj);
        ocp_nlp_out_get(p.config, p.dims, p.out, k, "u", utraj.data());
        u.push_back(utraj);
    }

    // Get the predicted state at the last time step N
    ocp_nlp_out_get(p.config, p.dims, p.out, N_, "x", xtraj.data());
    x.push_back(xtraj);
}
void acadosInterface::warmstart(const std::string &problem)
{
    Problem &p = *problems_.at(problem);
    std::vector<double> xtraj(p.nx);
    std::vector<double> utraj(p.nu);

    for (int k{0}; k < N_; k++)
    {
        ocp_nlp_out_get(p.config, p.dims, p.out, k, "x", xtraj.data());
        ocp_nlp_out_set(p.config, p.dims, p.out, k, "x", xtraj.data());
        ocp_nlp_out_get(p.config, p.dims, p.out, k, "u", utraj.data());
        ocp_nlp_out_set(p.config, p.dims, p.out, k, "u", utraj.data());
    }

    ocp_nlp_out_get(p.config, p.dims, p.out, N_, "x", xtraj.data());
    ocp_nlp_out_set(p.config, p.dims, p.out, N_, "x", xtraj.data());
}
void acadosInterface::waitWorkers(){
    // spin briefly (the solves are short) and then sleep, so that on machines with fewer cores
    // than threads the waiting thread does not steal the time of the workers
    static const int spin_iterations = [](){
        const char *value = std::getenv("DWMPC_SPIN");
        return value ? std::atoi(value) : 2000;
    }();
    for (auto &entry : problems_)
    {
        for (int i{0}; i < spin_iterations && entry.second->start.load(); i++) {}
        entry.second->start.wait(true);
    }
}
