import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


"""
generating gaussian approximation to total density that's defined in psp file
"""


def g_integration(A, b):
    temp = np.sqrt(np.pi / b)
    return A * temp * temp * temp


def print_get_integration(x0):
    integration = 0
    n = int(len(x0) / 2)
    for i in range(n):
        A = x0[i + n]
        b = 1 / (2 * x0[i] * x0[i])
        temp = g_integration(A, b)
        print("{}\t{}\t{}".format(A, b, temp))
        integration += temp
    print("integration")
    print(integration)
    return integration


def get_integration(x0):
    integration = 0
    n = int(len(x0) / 2)
    for i in range(n):
        A = x0[i + n]
        b = 1 / (2 * x0[i] * x0[i])
        temp = g_integration(A, b)
        integration += temp
    #     print("integration")
    #     print(integration)
    return integration


def get_integration_comp(x0):
    integration = 0
    int_abs = 0
    n = int(len(x0) / 2)
    for i in range(n):
        A = x0[i + n]
        b = 1 / (2 * x0[i] * x0[i])
        temp = g_integration(A, b)
        int_abs += abs(temp)
        integration += temp
    #     print("integration")
    #     print(integration)
    return integration, int_abs


def normalize(x0, num_elec):
    integration = get_integration(x0)
    n = int(len(x0) / 2)
    factor = num_elec / integration
    for i in range(n):
        x0[i + n] = x0[i + n] * factor
    return x0


def mse_function(x0, r, ref, num_elec):
    x0 = normalize(x0, num_elec)

    n = int(len(x0) / 2)
    result = np.zeros(len(r))
    for i in range(n):
        result += np.exp(-(1 / (2 * x0[i] * x0[i])) * r * r) * x0[i + n]
    return np.mean(np.square(result - ref))


def mae_function(x0, r, ref, num_elec):
    # x0 = normalize(x0, num_elec)

    n = int(len(x0) / 2)
    result = np.zeros(len(r))
    for i in range(n):
        result += np.exp(-(1 / (2 * x0[i] * x0[i])) * r * r) * x0[i + n]
    return np.mean(np.abs(result - ref))


def mae_function_regularized(x0, r, ref, num_elec):
    # x0 = normalize(x0, num_elec)

    n = int(len(x0) / 2)
    result = np.zeros(len(r))
    for i in range(n):
        result += np.exp(-(1 / (2 * x0[i] * x0[i])) * r * r) * x0[i + n]

    integral, integral_abs = get_integration_comp(x0)
    regularize = (integral - num_elec) ** 2 + 0.03 * integral_abs
    #     regularize = abs(get_integration(x0) - num_elec)
    return np.mean(np.abs(result - ref)) + 0.01 * regularize


def get_result(x0, r):
    n = int(len(x0) / 2)
    result = np.zeros(len(r))
    for i in range(n):
        result += np.exp(-(1 / (2 * x0[i] * x0[i])) * r * r) * x0[i + n]
    return result


def plot_comparison(r, ref, predict, plot_filename="test.png"):
    f, ax = plt.subplots(figsize=(20, 5))
    ax.plot(r, ref)
    ax.plot(r, predict)
    plt.savefig(plot_filename)
    plt.show()
    print(np.mean(np.square(predict - ref)))
    return


def predict_and_plot(x0, r, ref, plot_filename):
    predict = get_result(x0, r)
    plot_comparison(r, ref, predict, plot_filename)


def optimize_coeff(x0, r, ref, num_elec, show=True, plot_filename="test.png"):
    print("\n=========================start=================================\n")
    #     if show:
    #         predict_and_plot(x0, r, ref)
    #     print(x0)
    res = minimize(
        mae_function_regularized,
        x0,
        (r, ref, num_elec),
        method="nelder-mead",
        options={"maxiter": 100000, "disp": True, "xatol": 1e-11, "adaptive": True},
    )

    #     result = normalize(res.x, num_elec)
    result = res.x
    #     print(num_elec)
    integration = print_get_integration(result)

    if show:
        predict_and_plot(result, r, ref, plot_filename=plot_filename)
    print("\n==========================================================\n")
    return result, integration


def get_optimized_gaussian(r, ref, num_gaussian, num_elec, plot_filename="test.png"):
    x0 = np.concatenate((np.logspace(-1, 0.0, num=num_gaussian), np.ones(num_gaussian)))
    #     print(num_gaussian)
    optimized, integration = optimize_coeff(
        x0, r, ref, num_elec, plot_filename=plot_filename
    )
    error = mae_function(optimized, r, ref, num_elec)
    return optimized, error, integration


def save_gaussian(r, x0, atom, density_type, num_gaussian):
    #     log_filename = "info.log"
    result_filename = "{}_{}_{}.g".format(atom, density_type, num_gaussian)

    f = open(result_filename, "w")
    n = int(len(x0) / 2)
    for i in range(n):
        f.write("{}\t{}\n".format(x0[i + n], (1 / (2 * x0[i] * x0[i]))))
    f.close()

    return


def log(log_filename, message):
    f = open(log_filename, "a")
    f.write(message)
    f.close()
    return


def normalize_ref_data(r, density, num_elec):
    integration = 0
    for i in range(len(r) - 1):
        r_low = r[i]
        r_high = r[i + 1]
        rho_ave = (density[i] + density[i + 1]) / 2
        integration += (4 / 3) * np.pi * (r_high**3 - r_low**3) * rho_ave
    print("ref data integration:{}".format(integration))
    factor = num_elec / integration
    #     print("{}\t{}".format(integration, factor))
    density_normalized = density * factor
    return density_normalized


def optimize_atom_and_save(atom, num_elec, num_gaussian, log_filename="info.log"):
    csv_filename = "{}_pot.tsv".format(atom)
    data = []
    with open(csv_filename, "r") as dest_f:
        Lines = dest_f.readlines()
        for line in Lines:
            temp = line.strip().split()
            data.append(np.asarray(temp, dtype=np.float))
    data_array = np.array(data)

    r = data_array[:, 1] * 0.529177249
    pseudo_density = data_array[:, 2] / (4 * np.pi)
    core_density = data_array[:, 4] / (4 * np.pi)
    total_density = pseudo_density + core_density
    total_density = normalize_ref_data(r, total_density, num_elec)

    plot_filename = "{}_{}_totaldensity_comparison.png".format(atom, num_gaussian)
    optimized_gaussian, error, integration = get_optimized_gaussian(
        r, total_density, num_gaussian, num_elec, plot_filename=plot_filename
    )
    save_gaussian(
        r,
        optimized_gaussian,
        atom,
        density_type="totaldensity",
        num_gaussian=num_gaussian,
    )
    msg = "atom:{}\ttype:{}\tngaussian:{}\tmse:{}\tintegration:{}\n".format(
        atom, "totaldensity", num_gaussian, error, integration
    )
    print(msg)
    log(log_filename, msg)

    return


if __name__ == "__main__":
    for atom, num_elec in [
        ("H", 1),
        ("C", 6),
        ("N", 7),
        ("O", 8),
        ("F", 9),
        ("Pt", 78),
        ("Cu", 29),
        ("Au", 79),
        ("Rh", 45),
        ("Si", 14),
        ("Zr", 40),
        ("Zn", 30),
    ]:
        for num_gaussian in range(2, 9):
            optimize_atom_and_save(
                atom, num_elec, num_gaussian, log_filename="info.log"
            )
