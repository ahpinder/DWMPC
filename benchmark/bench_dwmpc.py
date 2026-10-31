"""Headless benchmark of the DWMPC controller: times mpc.run() and mpc.prepare().

Usage: python3 benchmark/bench_dwmpc.py [--duration 20] [--warmup 50] [--vel 0.3] [--out results.json]
"""
import argparse
import json
import os
import platform
import subprocess
import time

import numpy as np
import pydwmpc
from gym_quadruped.quadruped_env import QuadrupedEnv

parser = argparse.ArgumentParser()
parser.add_argument("--robot", default="aliengo")
parser.add_argument("--scene", default="flat")
parser.add_argument("--duration", type=float, default=20.0, help="simulated seconds")
parser.add_argument("--warmup", type=int, default=50, help="MPC calls excluded from timing stats")
parser.add_argument("--vel", type=float, default=0.3, help="forward velocity command [m/s]")
parser.add_argument("--mpc-freq", type=float, default=100.0)
parser.add_argument("--sim-freq", type=float, default=200.0)
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--out", default=None, help="optional JSON output path")
args = parser.parse_args()

np.random.seed(args.seed)
decimation = int(round(args.sim_freq / args.mpc_freq))

env = QuadrupedEnv(robot=args.robot,
                   scene=args.scene,
                   sim_dt=1.0 / args.sim_freq,
                   ref_base_lin_vel=args.vel,
                   ground_friction_coeff=1.5,
                   base_vel_command_type="forward",
                   state_obs_names=tuple(QuadrupedEnv.ALL_OBS))
env.reset(random=False)

mpc = pydwmpc.Dwmpc()
mpc.init()
mpc.startWalking()
mpc.setStepHeight(0.1)

tau = pydwmpc.DoubleVector([0.0] * 12)
des_q = pydwmpc.DoubleVector([0.0] * 12)
des_dq = pydwmpc.DoubleVector([0.0] * 12)
contact = pydwmpc.DoubleVector([0.0] * 4)
Kp, Kd = 20, 3

run_ns, prep_ns = [], []
start_x = float(env.mjData.qpos[0])
fell = False
steps_done = 0

for step in range(int(args.duration * args.sim_freq)):
    qpos = env.mjData.qpos
    qvel = env.mjData.qvel

    if step % decimation == 0:
        feet = env.feet_pos('world')
        foot_op = np.array([feet.FL, feet.FR, feet.RL, feet.RR], order="F")
        cs = env.feet_contact_state()[0]
        contact_op = np.array([float(cs.FL), float(cs.FR), float(cs.RL), float(cs.RR)])
        quat = np.array([qpos[4], qpos[5], qpos[6], qpos[3]])
        ref_lin, ref_ang = env.target_base_vel()
        omega = env.base_configuration[:3, :3] @ qvel[3:6]

        t0 = time.perf_counter_ns()
        mpc.run(qpos[:3].copy(), quat, qpos[7:].copy(), qvel[:3].copy(), omega, qvel[6:].copy(),
                1 / args.mpc_freq, contact_op, foot_op,
                env.heading_orientation_SO3.transpose() @ ref_lin, ref_ang,
                np.array([0.0, 0.0, 0.0, 1.0]), contact, tau, des_q, des_dq)
        t1 = time.perf_counter_ns()
        mpc.prepare()
        t2 = time.perf_counter_ns()
        run_ns.append(t1 - t0)
        prep_ns.append(t2 - t1)

    action_torque = tau + Kp * (des_q.getList() - qpos[7:]) + Kd * (des_dq.getList() - qvel[6:])
    action = np.zeros(env.mjModel.nu)
    action[env.legs_tau_idx.FL] = action_torque[:3]
    action[env.legs_tau_idx.FR] = action_torque[3:6]
    action[env.legs_tau_idx.RL] = action_torque[6:9]
    action[env.legs_tau_idx.RR] = action_torque[9:]
    _, _, terminated, _, _ = env.step(action=action)
    steps_done = step + 1
    if terminated:
        fell = True
        break

distance_x = float(env.mjData.qpos[0]) - start_x
env.close()


def stats(ns, period_ms):
    ms = np.array(ns[args.warmup:], dtype=float) / 1e6
    return {
        "n": int(ms.size),
        "mean_ms": float(ms.mean()), "median_ms": float(np.median(ms)),
        "p95_ms": float(np.percentile(ms, 95)), "p99_ms": float(np.percentile(ms, 99)),
        "max_ms": float(ms.max()), "std_ms": float(ms.std()),
        "frac_over_period": float((ms > period_ms).mean()),
    }


if len(run_ns) <= args.warmup:
    raise SystemExit("Too few MPC calls recorded (robot likely fell early); no stats.")

period_ms = 1000.0 / args.mpc_freq
try:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
except Exception:
    commit = os.environ.get("GIT_COMMIT", "unknown")
cpu = "unknown"
try:
    cpu = next(l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name"))
except Exception:
    pass

result = {
    "config": vars(args),
    "commit": commit,
    "host": {"machine": platform.machine(), "cpu_count": os.cpu_count(), "cpu": cpu},
    "mpc_period_ms": period_ms,
    "first_call_ms": (run_ns[0] + prep_ns[0]) / 1e6,
    "run": stats(run_ns, period_ms),
    "prepare": stats(prep_ns, period_ms),
    "run_plus_prepare": stats([a + b for a, b in zip(run_ns, prep_ns)], period_ms),
    "quality": {"fell": fell, "sim_time_s": steps_done / args.sim_freq, "distance_x_m": distance_x},
}

print(f"first (initialisation) call: {result['first_call_ms']:.1f} ms")
print(f"MPC calls timed: {result['run']['n']} (period {period_ms:.1f} ms), fell={fell}, "
      f"distance_x={distance_x:.2f} m over {result['quality']['sim_time_s']:.1f} s")
for key in ("run", "prepare", "run_plus_prepare"):
    s = result[key]
    print(f"{key:>17}: mean {s['mean_ms']:.2f}  median {s['median_ms']:.2f}  p95 {s['p95_ms']:.2f}  "
          f"p99 {s['p99_ms']:.2f}  max {s['max_ms']:.2f} ms  over-period {100 * s['frac_over_period']:.1f}%")
if args.out:
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
