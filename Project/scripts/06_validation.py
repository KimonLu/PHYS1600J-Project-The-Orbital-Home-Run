#!/usr/bin/env python3
from __future__ import annotations
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from orbital_home_run import C, propagate, acceleration_central, horizontal_surface_orbit
from plotting import setup, save, OUT


def main() -> None:
    setup()
    u = 1.2
    analytic = horizontal_surface_orbit(u)
    state0 = np.array([C.radius, 0.0, 0.0, 0.0, u*C.circular_speed, 0.0])
    rows = []
    steps = [60, 30, 15, 8, 4, 2]
    for step in steps:
        tr = propagate(state0, analytic["period"], acceleration_central, max_step=step, rtol=2e-11, atol=2e-6, samples=2001)
        s = tr["state"]
        pos_err = np.linalg.norm(s[-1,:3]-state0[:3])
        vel_err = np.linalg.norm(s[-1,3:]-state0[3:])
        r = np.linalg.norm(s[:,:3],axis=1)
        v2 = np.sum(s[:,3:]**2,axis=1)
        energy = 0.5*v2-C.mu/r
        h = np.linalg.norm(np.cross(s[:,:3],s[:,3:]),axis=1)
        rows.append([step,pos_err,vel_err,np.max(np.abs((energy-energy[0])/energy[0])),np.max(np.abs((h-h[0])/h[0])),r.min()-analytic["rp"],r.max()-analytic["ra"]])
    df = pd.DataFrame(rows,columns=["max_step_s","position_closure_error_m","velocity_closure_error_m_s","max_relative_energy_error","max_relative_angular_momentum_error","periapsis_error_m","apoapsis_error_m"])
    df.to_csv(OUT/"numerical_convergence.csv",index=False)

    fig,ax=plt.subplots(figsize=(3.45,2.8))
    ax.loglog(df["max_step_s"],df["position_closure_error_m"],marker="o")
    ax.invert_xaxis()
    ax.set_xlabel("Maximum integrator step (s)")
    ax.set_ylabel("One-period closure error (m)")
    ax.set_title("Numerical convergence against the analytic orbit")
    save(fig,"fig12_numerical_convergence")


if __name__=="__main__":
    main()
