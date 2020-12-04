import hashlib
from ..util import list_symbols_to_indices
import numpy as np


class GaussianDescriptorSet:
    def __init__(
        self,
        elements,
        cutoff=6.5,
        cutoff_params={"cutoff_func": "cosine"},
        default_interactions=False,
    ):
        self.elements = list(elements)
        self.element_indices = list_symbols_to_indices(elements)
        self.cutoff = cutoff
        cutoff_params = {
            (k.lower() if isinstance(k, str) else k): (
                v.lower() if isinstance(v, str) else v
            )
            for k, v in cutoff_params.items()
        }
        self.cutoff_params = cutoff_params
        self.all_interactions = set()
        self.interactions = {
            element: {"G2": set(), "G4": set(), "G5": set()} for element in elements
        }
        self.descriptor_setup = None
        self.descriptor_setup_hash = None

    def batch_add_descriptors(
        self, number, param1s, param2s, param3s, cutoff=None, update=True
    ):
        for element_i in self.elements:
            for j, element_j in enumerate(self.elements):
                if number == 2:
                    for eta, rs in zip(param1s, param2s):
                        self.add_g2(element_i, element_j, eta, rs, cutoff, False)
                else:
                    for element_k in self.elements[j:]:
                        for eta, zeta, gamma in zip(param1s, param2s, param3s):
                            if number == 4:
                                self.add_g4(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    False,
                                )
                            else:
                                self.add_g5(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    False,
                                )
        if update:
            self.update()

    def add_g2(self, element_i, element_j, eta=3.0, rs=0.0, cutoff=None, update=True):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        g2_params = (
            2,
            self.element_indices[self.elements.index(element_j)],
            0,
            cutoff or self.cutoff,
            eta / (cutoff or self.cutoff) ** 2.0,
            rs,
            0.0,
        )
        # print(element_i, g2_params)
        self.interactions[element_i]["G2"].add(g2_params)
        if update:
            self.update()
        return self.interactions[element_i]["G2"]

    def add_g4(
        self,
        element_i,
        element_j,
        element_k,
        eta=0.005,
        zeta=1.0,
        gamma=1.0,
        cutoff=None,
        update=True,
    ):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        assert element_k in self.elements, f"{element_k} is not in {self.elements}"
        # element_j, element_k = sorted([element_j, element_k])
        g4_params = (
            4,
            self.element_indices[self.elements.index(element_j)],
            self.element_indices[self.elements.index(element_k)],
            cutoff or self.cutoff,
            eta / (cutoff or self.cutoff) ** 2.0,
            zeta,
            gamma,
        )
        # print(element_i, g4_params)
        self.interactions[element_i]["G4"].add(g4_params)
        if update:
            self.update()
        return self.interactions[element_i]["G4"]

    def add_g5(
        self,
        element_i,
        element_j,
        element_k,
        eta=0.005,
        zeta=1.0,
        gamma=1.0,
        cutoff=None,
        update=True,
    ):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        assert element_k in self.elements, f"{element_k} is not in {self.elements}"
        # element_j, element_k = sorted([element_j, element_k])
        g5_params = (
            5,
            self.element_indices[self.elements.index(element_j)],
            self.element_indices[self.elements.index(element_k)],
            cutoff or self.cutoff,
            eta,
            zeta,
            gamma,
        )
        # print(element_i, g5_params)
        self.interactions[element_i]["G5"].add(g5_params)
        self.descriptor_setup = None
        self.descriptor_setup_hash = None
        if update:
            self.update()
        return self.interactions[element_i]["G5"]

    def update(self):
        self.descriptor_setup = self._get_descriptor_setup()
        self.descriptor_setup_hash = self._get_descriptor_setup_hash()

    def process_combinatorial_Gs(self, Gs):
        for element in self.elements:
            if element in Gs:
                # self.interactions[element] =
                self._process_element_combinatorial_params(element, Gs[element])
            elif "default" in Gs:
                # self.interactions[element] =
                self._process_element_combinatorial_params(element, Gs["default"])
            else:
                raise NotImplementedError(
                    "Symmetry function parameters not defined properly - element {} passed but not present in {}".format(
                        element, self.elements
                    )
                )
        self.update()

    def _process_element_combinatorial_params(self, element_i, element_Gs):
        cutoff = element_Gs["cutoff"]
        if "G2" in element_Gs:
            for eta in np.array(element_Gs["G2"]["etas"]):
                for rs in element_Gs["G2"]["rs_s"]:
                    for element_j in self.elements:
                        self.add_g2(element_i, element_j, eta, rs, cutoff, update=False)

        if "G4" in element_Gs:
            for eta in np.array(element_Gs["G4"]["etas"]):
                for zeta in element_Gs["G4"]["zetas"]:
                    for gamma in element_Gs["G4"]["gammas"]:
                        for j, element_j in enumerate(self.elements):
                            for element_k in self.elements[j:]:
                                self.add_g4(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    update=False,
                                )

        if "G5" in element_Gs:
            for eta in element_Gs["G5"]["etas"]:
                for zeta in element_Gs["G5"]["zetas"]:
                    for gamma in element_Gs["G5"]["gammas"]:
                        for j, element_j in enumerate(self.elements):
                            for element_k in self.elements[j:]:
                                self.add_g5(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    update=False,
                                )

    def _get_descriptor_setup_hash(self):
        string = (
            "cosine"
            if self.cutoff_params["cutoff_func"] == "cosine"
            else "polynomial%.15f" % self.cutoff_params["gamma"]
        )
        descriptor_setup = self.descriptor_setup
        for element in descriptor_setup.keys():
            string += element
            for desc in descriptor_setup[element]:
                for num in desc:
                    string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        return hash_result

    def _get_descriptor_setup(self):
        descriptor_setup = {element: None for element in self.elements}
        for element, descriptors in self.interactions.items():
            g2s, g4s, g5s = descriptors["G2"], descriptors["G4"], descriptors["G5"]
            g2s = [list(params) for params in sorted(g2s)]
            g4s = [list(params) for params in sorted(g4s)]
            g5s = [list(params) for params in sorted(g5s)]
            descriptor_setup[element] = np.array(g2s + g4s + g5s)
        return descriptor_setup

    def __eq__(self, other):
        return self.descriptor_setup_hash == other.descriptor_setup_hash

    def __hash__(self):
        return self.descriptor_setup_hash

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        interaction_counts = []
        for element, intxns in self.interactions.items():
            interaction_counts.append(
                "%s: {#G2: %d, #G4: %d, #G5: %d}"
                % (element, len(intxns["G2"]), len(intxns["G4"]), len(intxns["G5"]))
            )
        return "GaussianDescriptorSet(%s)" % "".join(interaction_counts)
