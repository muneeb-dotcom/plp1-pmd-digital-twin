import json
import pandas as pd
import numpy as np

class PMDTwin:
    """Mechanism-to-phenotype model for PLP1 variants.

    NOT a mechanistic simulation. Reports measured cell-state module
    values (log2 expression) relative to a healthy mature-OL baseline
    built from the same units. Extrapolates to unmeasured mechanism
    classes only via the explicit gate functions below, never by
    fabricating a data point.
    """

    def __init__(self, state_modules_path, disease_states_path):
        with open(state_modules_path) as f:
            ref = json.load(f)
        self.module_genes = ref["genes"]
        self.healthy_baseline = ref["healthy_baseline"]
        self.disease_states = pd.read_csv(disease_states_path)

    def predict_state(self, mechanism):
        subset = self.disease_states[self.disease_states["mechanism"] == mechanism]
        if subset.empty:
            return {"status": "no_measured_data", "mechanism": mechanism}
        result = {"status": "measured", "mechanism": mechanism}
        for _, row in subset.iterrows():
            result[row["module"]] = {
                "wt_mean_log2fpkm": row["wt_mean_log2fpkm"],
                "disease_mean_log2fpkm": row["jimpy_mean_log2fpkm"],
                "delta_log2fpkm": row["delta_log2fpkm"],
            }
        return result

    @staticmethod
    def _sigmoid(x):
        return 1 / (1 + np.exp(-x))

    def naive_myelination(self, maturation_z, apoptosis_z):
        return self._sigmoid(maturation_z) * self._sigmoid(-apoptosis_z)

    def gated_myelination(self, maturation_z, apoptosis_z, myelin_z, er_stress_z):
        survival = self.naive_myelination(maturation_z, apoptosis_z)
        program = self._sigmoid(myelin_z) * self._sigmoid(-er_stress_z)
        return survival * program

    def explain(self, mechanism=None):
        return self.module_genes

    def available_mechanisms(self):
        return self.disease_states["mechanism"].unique().tolist()


if __name__ == "__main__":
    twin = PMDTwin("results/state_modules.json", "results/pmd_state_vectors.csv")
    print("Available:", twin.available_mechanisms())
    print()
    for k, v in twin.predict_state("misfolding_confirmed").items():
        print(f"  {k}: {v}")
    print()
    print("duplication (no data):", twin.predict_state("duplication"))