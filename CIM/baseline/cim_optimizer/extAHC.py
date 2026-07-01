import torch
import numpy as np
import copy
torch.backends.cudnn.benchmark = True


def CIM_ext_f_AHC_GPU(time_stop, J, h, batch_size=1, nsub=1, dt=0.015625,
                       F_h=2.0, alpha=1.0, beta=0.05, delta=0.25, eps=0.333,
                       lambd=0.001, pi=-0.225, rho=1.0, tau=100, noise=0,
                       ahc_ext_nonlinearity=None, custom_fb_schedule=None,
                       custom_pump_schedule=None, device=torch.device('cpu')):
    """CIM solver with amplitude-heterogeneity correction and external-field.

    Reference:
    T. Leleu, Y. Yamamoto, P. L. McMahon, and K. Aihara, Destabilization
    of local minima in analog spin systems by correction of amplitude
    heterogeneity. Phys Rev Lett 122, 040607 (2019).
    """
    J = torch.from_numpy(J).float().to(device)
    h = torch.from_numpy(h).float().to(device)
    N = J.size()[1]
    if h is None:
        h = torch.zeros(batch_size, N).float().to(device)

    energy_plot_data = torch.zeros(batch_size, time_stop).to(device)
    error_var_data = torch.zeros(batch_size, N, time_stop).to(device)
    spin_amplitude_trajectory = torch.zeros(batch_size, N, time_stop).to(device)
    MAX_FLOAT = 500
    EMAX_FLOAT = 32

    theta0 = (2 * (torch.rand(batch_size, N) - 0.5)).to(device)
    x = (0.14 * (torch.cos(2 * torch.pi * theta0))).to(device)

    sig = ((2 * (x > 0) - 1).float()).to(device)
    mu = x @ J
    h_ = sig @ J

    H = (-1/2 * (torch.bmm(
        sig.view(batch_size, 1, N),
        (sig @ J).view(batch_size, N, 1)))[:, :, 0]).view(batch_size) - sig @ h
    e = (0.9 * torch.ones(batch_size, N) + 0.6 * torch.rand(batch_size, N)).to(device)
    g = (0.9 * torch.ones(batch_size, N) + 0.6 * torch.rand(batch_size, N)).to(device)

    ic = torch.ones(batch_size).to(device)
    beta = torch.zeros(batch_size).to(device)
    p = pi * torch.ones(batch_size).to(device)
    a = alpha * torch.ones(batch_size).to(device)

    if custom_fb_schedule is None:
        eps = (torch.ones(time_stop) * eps).to(device)
    else:
        eps = custom_fb_schedule((torch.arange(0, time_stop, dt))).to(device)
    if custom_pump_schedule is None:
        pi = (torch.ones(time_stop) * pi).to(device)
    else:
        pi = custom_pump_schedule(torch.arange(0, time_stop)).to(device)

    rho = 0
    if ahc_ext_nonlinearity is None:
        ahc_ext_nonlinearity = lambda c: torch.pow(c, 3)

    H_opt = H
    sig_opt = sig
    t_opt = torch.zeros(batch_size).float().to(device)
    Jm = ((abs(J) > 0).int()).to(device)
    norm_ = torch.ones(N).to(device)
    norm = (norm_.float() @ Jm.float()).to(device)

    dt_sub = dt / nsub
    effective_tau = tau / dt

    for t in range(time_stop):
        spin_amplitude_trajectory[:, :, t] = x
        error_var_data[:, :, t] = e
        x_prev = copy.deepcopy(x)
        e_prev = copy.deepcopy(e)
        g_prev = copy.deepcopy(g)

        for z in range(nsub):
            x_squared = x ** 2
            dxdt = (p - 1)[:, None] * x - ahc_ext_nonlinearity(x)
            dxdt += noise * (torch.rand(batch_size, N, device=device) - 0.5)
            dedt = -beta[:, None] * (x_squared - a[:, None]) * e_prev
            g_ = 5 * g_prev
            xmax = torch.sqrt(a / 2)
            theta = torch.log(torch.cosh(g_ * x_prev)) / g_
            theta += torch.abs(xmax)[:, None]
            theta += -torch.log(torch.cosh(g_ * xmax[:, None])) / g_

            mx = F_h * (theta @ Jm.float()) / norm
            dtheta = x_prev ** 2 - a[:, None] / 2
            dgdt = beta[:, None] * (dtheta * g_prev)

            x = x + dxdt * dt_sub
            e = e + dedt * dt_sub
            g = g + dgdt * dt_sub

            g[g > EMAX_FLOAT] = EMAX_FLOAT
            e[e > EMAX_FLOAT] = EMAX_FLOAT
            e = torch.abs(e)
            g = torch.abs(g)

        x = x + (eps[t] * dt) * e_prev * (x_prev @ J)
        x = x + (eps[t] * dt) * e_prev * mx * h / 2
        sig = ((2 * (x > 0) - 1).float())
        h_ = sig @ J

        H = (-1/2 * (torch.bmm(
            sig.view(batch_size, 1, N),
            (h_).view(batch_size, N, 1)))[:, :, 0]).view(batch_size) - sig @ h
        energy_plot_data[:, t] = H

        dH = (H - H_opt)
        a = alpha + rho * torch.tanh(delta * dH)
        p = pi[t] - rho * torch.tanh(delta * dH)
        beta = beta + lambd * dt

        ic[ic < t - effective_tau] = t
        beta[ic < t - effective_tau] = 0
        t_opt[H < H_opt] = t
        ic[H < H_opt] = t
        H_opt = torch.minimum(H_opt, H)

    H_opt = H_opt.cpu()
    spin_amplitude_trajectory = spin_amplitude_trajectory.cpu()
    t = time_stop
    spin_plot_data = 2 * (spin_amplitude_trajectory > 0) - 1
    for k in range(batch_size):
        sig_opt[k, :] = spin_plot_data[k, :, t_opt.long()[k]]
    energy_plot_data = energy_plot_data.cpu()
    sig_opt = sig_opt.cpu()
    error_var_data = error_var_data.cpu()
    return (sig_opt.numpy(), spin_amplitude_trajectory.numpy(), t,
            energy_plot_data.numpy(), error_var_data.numpy())
