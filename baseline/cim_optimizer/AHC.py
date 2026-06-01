import torch
import numpy as np
torch.backends.cudnn.benchmark = True


def CIM_AHC_GPU(T_time, J, batch_size=1, time_step=0.05, r=0.2, beta=0.05,
                eps=0.07, mu=1, noise=0, custom_fb_schedule=None,
                custom_pump_schedule=None, random_number_function=None,
                ahc_nonlinearity=None, device=torch.device('cpu')):
    """CIM solver with amplitude-heterogeneity correction; no external h field.

    Reference:
    T. Leleu, Y. Yamamoto, P. L. McMahon, and K. Aihara, Destabilization
    of local minima in analog spin systems by correction of amplitude
    heterogeneity. Phys Rev Lett 122, 040607 (2019).
    """
    J = torch.from_numpy(J).float().to(device)
    N = J.size()[1]

    end_ising_energy = (1e20 * torch.ones(batch_size)).to(device)
    target_a_baseline = 0.2
    target_a = (target_a_baseline * torch.ones(batch_size)).to(device)

    ticks = int(T_time / time_step)
    spin_amplitude_trajectory = torch.zeros(batch_size, N, ticks).to(device)
    error_var_data = torch.zeros(batch_size, N, ticks).to(device)
    energy_plot_data = torch.zeros(batch_size, ticks).to(device)
    t_opt = torch.zeros(batch_size).to(device)
    EMAX_FLOAT = 32.0

    x = 0.001 * torch.rand(batch_size, N).to(device) - 0.0005
    error_var = torch.ones(batch_size, N).to(device).float()
    etc_flag = torch.ones(batch_size, N).to(device)
    sig_ = ((2 * (x > 0) - 1).float()).to(device)
    sig_opt = sig_

    if random_number_function is None:
        random_number_function = lambda c: torch.rand(c, 3)

    if custom_fb_schedule is None:
        eps = (torch.ones(ticks) * eps).to(device)
    else:
        eps = custom_fb_schedule((torch.arange(0, ticks, time_step))).to(device)
    if custom_pump_schedule is None:
        r = (torch.ones(ticks) * r).to(device)
    else:
        r = (custom_pump_schedule(np.arange(0, ticks, time_step))).to(device)

    if ahc_nonlinearity is None:
        ahc_nonlinearity = lambda c: torch.pow(c, 3)

    for t in range(ticks):
        sig = ((2 * (x > 0) - 1).float())
        sig_ = sig

        curr_ising_energy = (-1/2 * (torch.bmm(
            sig.view(batch_size, 1, N),
            (sig @ J).view(batch_size, N, 1)))[:, :, 0]).view(batch_size)
        energy_plot_data[:, t] = curr_ising_energy

        spin_amplitude_trajectory[:, :, t] = x
        error_var_data[:, :, t] = error_var
        x_squared = x ** 2
        MVM = x @ J
        x += time_step * (x * ((r[t] - 1) - mu * x_squared))
        x += time_step * eps[t] * (MVM * error_var)
        x += eps[t] * noise * (torch.rand(N, device=device) - 0.5)

        delta_a = eps[t] * torch.mean((sig @ J) * sig * etc_flag, 1)
        target_a = target_a_baseline + delta_a
        x_squared = x ** 2

        error_var += time_step * (-beta * ((x_squared) - target_a[:, None]) * error_var)
        error_var[error_var > EMAX_FLOAT] = EMAX_FLOAT

        comparison = torch.any(sig_ != sig, 1)
        etc_flag[comparison, :] = error_var[comparison, :]
        t_opt[curr_ising_energy < end_ising_energy] = t
        end_ising_energy = torch.minimum(end_ising_energy, curr_ising_energy)

    sig = ((2 * (x > 0) - 1).float())
    spin_amplitude_trajectory = spin_amplitude_trajectory.cpu()
    spin_plot_data = 2 * (spin_amplitude_trajectory > 0) - 1
    energy_plot_data = energy_plot_data.cpu()
    error_var_data = error_var_data.cpu()
    for k in range(batch_size):
        sig_opt[k, :] = spin_plot_data[k, :, t_opt.long()[k]]
    sig_opt = sig_opt.cpu()

    return (sig_opt.numpy(), spin_amplitude_trajectory.numpy(), t,
            energy_plot_data.numpy(), error_var_data.numpy())
