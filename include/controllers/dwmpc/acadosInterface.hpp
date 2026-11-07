#include "acados_c/ocp_nlp_interface.h"

#include <vector>
#include <map>
#include <memory>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <stdlib.h>
#include <csignal>
#include <iostream>
#include <thread>
#include <atomic>
#include <string>
#include <stdexcept>

// Entry points of one generated acados solver, so that any generated problem
// (front, back, leg0, ...) can be handled by name.
struct SolverOps
{
    void *(*create_capsule)();
    int (*create)(void *capsule);
    int (*free)(void *capsule);
    int (*free_capsule)(void *capsule);
    int (*solve)(void *capsule);
    int (*update_params)(void *capsule, int stage, double *value, int np);
    ocp_nlp_config *(*get_config)(void *capsule);
    ocp_nlp_dims *(*get_dims)(void *capsule);
    ocp_nlp_in *(*get_in)(void *capsule);
    ocp_nlp_out *(*get_out)(void *capsule);
    ocp_nlp_solver *(*get_solver)(void *capsule);
    void *(*get_opts)(void *capsule);
};

struct SolverTiming
{
    double wall_ms{0};
    double worker_wall_sum_ms{0};
    double worker_wall_max_ms{0};
    double worker_cpu_sum_ms{0};
};

class acadosInterface {
    
    public:

        acadosInterface();
        ~acadosInterface();
        // create the solvers of the given problems ("wb" is not a solver and is ignored)
        void init(const int N, const std::vector<std::string> &problems);
        // setter
        void setInitialCondition(std::vector<double> &initial_condition,const std::string &problem);
        void warmstart(const std::string &problem);
        void setParameter(std::vector<std::vector<double>> &parameter,const  int np, const std::string &problem);
        void setConstraints(std::vector<std::vector<double>> &constraint, const int nlh0, const int nlh, const std::string &problem);
        void setWeight(std::vector<std::vector<double>> &weight,const int &n_w0,const int &n_w, const int &n_we,const std::string &problem);
        void setReference(std::vector<std::vector<double>> &reference,const int n_ref0, const int n_ref, const int n_ref_e, const std::string &problem);
        void solve();
        void prepare();
        void enableTiming(bool enabled) { timing_enabled_ = enabled; }
        const SolverTiming &lastTiming() const { return last_timing_; }
        // getter
        void UpdatePrediction(std::vector<std::vector<double>> &x,std::vector<std::vector<double>> &u, const std::string &problem);   

    private:

        struct Problem
        {
            const SolverOps *ops{nullptr};
            void *capsule{nullptr};
            ocp_nlp_config *config{nullptr};
            ocp_nlp_dims *dims{nullptr};
            ocp_nlp_in *in{nullptr};
            ocp_nlp_out *out{nullptr};
            ocp_nlp_solver *solver{nullptr};
            void *opts{nullptr};
            int nx{0};
            int nu{0};
            int status{1};
            double wall_ms{0};
            double cpu_ms{0};
            int clock_error{0};
            std::atomic_bool start{false};
            std::thread worker;
        };

        std::map<std::string,std::unique_ptr<Problem>> problems_;

        int N_{0}; // prediction horizon
        bool timing_enabled_{false};
        SolverTiming last_timing_{};

        void work(Problem &problem);
        void runAll(int rti_phase);
        void waitWorkers();
        void clear();
        std::atomic_bool run{true};
};
