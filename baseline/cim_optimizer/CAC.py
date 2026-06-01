import torch
import numpy as np
import math
import time
torch.backends.cudnn.benchmark = True


def CIM_CAC_GPU(T_time, J, batch_size=1, time_step=0.05, r=None, alpha=3.0,
                beta=0.25, gamma=0.00011, delta=10, mu=1, rho=3, tau=1000,
                noise=0, H0=None, stop_when_solved=False, num_sol=10,
                custom_fb_schedule=None, custom_pump_schedule=None,
                cac_nonlinearity=torch.tanh, device=torch.device('cpu')):
    """CIM solver with chaotic amplitude control; no external h field.

    Reference:
    Leleu, T., Khoyratee, F., Levi, T. et al. Scaling advantage of chaotic
    amplitude control for high-performance combinatorial optimization.
    Commun Phys 4, 266 (2021).
    """
    MAX_FLOAT = 500
    EMAX_FLOAT = 32
    J = torch.from_numpy(J).float().to(device)
    N = J.size()[1]
    if r is None:
        r = 0.8 - (N / 220) ** 2

    xi = torch.zeros(batch_size).to(device)
    t_c = torch.zeros(batch_size).to(device)
    H = torch.zeros(batch_size).to(device)
    spin_amplitude_trajectory = torch.zeros(batch_size, N, T_time).to(device)
    error_var_data = torch.zeros(batch_size, N, T_time).to(device)
    energy_plot_data = torch.zeros(batch_size, T_time).to(device)
    t_opt = torch.zeros(batch_size).to(device)

    x = 0.001 * torch.rand(batch_size, N).to(device) - 0.0005
    error_var = torch.ones(batch_size, N).to(device).float()
    effective_tau = tau / time_step

    if custom_fb_schedule is None:
        beta = torch.ones(T_time).to(device) * beta
    else:
        beta = custom_fb_schedule(torch.arange(0, T_time).to(device))
    if custom_pump_schedule is None:
        r = torch.ones(T_time).to(device) * r
    else:
        r = custom_pump_schedule(torch.arange(0, T_time).to(device))

    sig = ((2 * (x > 0) - 1).float()).to(device)
    H = (-1/2 * (torch.bmm(
        sig.view(batch_size, 1, N),
        (sig @ J).view(batch_size, N, 1)))[:, :, 0]).view(batch_size)

    H_opt = H
    sig_opt = sig
    a = alpha * torch.ones(batch_size).to(device)

    for t in range(T_time):
        x_ = x
        spin_amplitude_trajectory[:, :, t] = x_
        error_var_data[:, :, t] = error_var
        sig = ((2 * (x > 0) - 1).float())
        H = (-1/2 * (torch.bmm(
            sig.view(batch_size, 1, N),
            (sig @ J).view(batch_size, N, 1)))[:, :, 0]).view(batch_size)

        x_squared = x ** 2
        MVM = x @ J
        x += time_step * (x * ((r[t] - 1) - mu * x_squared))
        x += time_step * beta[t] * (MVM * error_var)
        x += beta[t] * noise * (torch.rand(batch_size, N, device=device) - 0.5)

        energy_plot_data[:, t] = H

        delta_e = -xi[:, None] * (x_ ** 2 - a[:, None]) * error_var
        error_var += delta_e * time_step
        error_var[error_var > EMAX_FLOAT] = EMAX_FLOAT

        xi += gamma * time_step
        dH = H - H_opt
        a = alpha + rho * cac_nonlinearity(delta * dH)

        t_c[t_c < t - effective_tau] = t
        xi[t_c < t - effective_tau] = 0
        t_opt[H < H_opt] = t
        t_c[H < H_opt] = t
        H_opt = torch.minimum(H_opt, H)

    spin_amplitude_trajectory = spin_amplitude_trajectory.cpu()
    spin_plot_data = 2 * (spin_amplitude_trajectory > 0) - 1
    energy_plot_data = energy_plot_data.cpu()

    for k in range(batch_size):
        sig_opt[k, :] = spin_plot_data[k, :, t_opt.long()[k]]
    sig_opt = sig_opt.cpu()
    error_var_data = error_var_data.cpu()
    return (sig_opt.numpy(), spin_amplitude_trajectory.numpy(), t,
            energy_plot_data.numpy(), error_var_data.numpy())
