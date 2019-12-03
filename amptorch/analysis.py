import copy
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parity_plot(calc, images, label, data="energy"):
    # TODO Parity plots without the calc.get_property(). Batch computation via
    # model forward pass
    images = copy.deepcopy(images)
    """Constructs a parity plot"""
    fig = plt.figure(figsize=(7.0, 7.0))
    ax = fig.add_subplot(111)
    targets = []
    preds = []
    if not os.path.exists("results/plots/parity"):
        os.makedirs("results/plots/parity")
    if data == "energy":
        for image in images:
            targets.append(image.get_potential_energy())
            image.set_calculator(calc)
            preds.append(image.get_potential_energy())
        targets = np.array(targets).reshape(-1, 1)
        preds = np.array(preds).reshape(-1, 1)
        energy_min = min(targets)
        energy_max = max(targets)
        ax.plot(targets, preds, "bo", markersize=3)
        ax.plot([energy_min, energy_max], [energy_min, energy_max], "r-", lw=0.3)
        ax.set_xlabel("ab initio energy, eV")
        ax.set_ylabel("PyTorch energy, eV")
        ax.set_title("Energies")
        fig.savefig("results/plots/parity/parity_E_"+label+".pdf")
    if data == "forces":
        for image in images:
            targets.append(image.get_forces().reshape(-1, 1))
            image.set_calculator(calc)
            preds.append(image.get_forces().reshape(-1, 1))
        targets = np.array(targets).reshape(-1, )
        preds = np.array(preds).reshape(-1, )
        force_min = min(targets)
        force_max = max(targets)
        ax.plot(targets, preds, "bo", markersize=3)
        ax.plot([force_min, force_max], [force_min, force_max], "r-", lw=0.3)
        ax.set_xlabel("ab initio forces, eV/A")
        ax.set_ylabel("PyTorch forces, eV/A")
        ax.set_title("Forces")
        fig.savefig("results/plots/parity/parity_F_"+label+".pdf")
