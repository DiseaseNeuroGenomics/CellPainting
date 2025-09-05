from typing import List

import os
import copy
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.cross_decomposition import PLSRegression

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


class CreateData:

    def __init__(
        self,
        plate_schema_fn: str,
        drug_info_fn: str,
        plate_directories: List[str],
        include_phagocytosis: bool = True,
    ):

        # get the mapping between the drug code and well row and column
        self._extract_drug_info_from_plate_schema(plate_schema_fn)

        # get the mapping between the drug code and the drug name and concentration
        self._extract_drug_info_from_drug_conditions(drug_info_fn)

        # load Cell Painting, and possibly phagocytosis, data
        self._load_plate_directories(plate_directories, include_phagocytosis=include_phagocytosis)

        # given the mapping between drug codes and names/concentrations, add names and conentrations
        self._add_drug_info()

        # list of features to be added to adata.obs
        self.observation_features = [
            "Plate", "Drug", "Concentration", "All Cells - Number of Objects",
            "Non border cells - Number of Objects",
            "Non border cells - Total Spot Area - Mean per Well",
            "Non border cells - Number of Spots - Mean per Well",
            "Spots Substrate - Spot Area [px²] - Mean per Well",
            "Spots Substrate - Relative Spot Intensity - Mean per Well",
        ]

        # list of Phagocytosis features we wish to normalize
        self.features_to_normalize = [
                "All Cells - Number of Objects",
                "Non border cells - Number of Objects",
                "Non border cells - Total Spot Area - Mean per Well",
                "Non border cells - Number of Spots - Mean per Well",
                "Spots Substrate - Spot Area [px²] - Mean per Well",
                "Spots Substrate - Relative Spot Intensity - Mean per Well",
                "Phagocytosis per cell",
        ]


    def _extract_drug_info_from_plate_schema(self, plate_schema_fn):

        """
        Example plate_schema_fn = "data/Plate_setup_SchemaA_11Dx3x7conc.xlsx"
        Excel file that links the drug codes (e.g. Drug8_c5) with the well row and column
        """

        df = pd.read_excel(plate_schema_fn, sheet_name='destination_plate_schema')
        _, n_cols = df.shape
        n_rows = 18 # extra info provided lower down we want to avoid
        self.drug_map = {i: {} for i in range(1, n_rows)}
        for i in range(n_rows):
            for j in range(n_cols):
                code = df.iloc[i, j]
                if isinstance(code, str) and ("Drug" in code or "c" in code):
                    self.drug_map[i + 1][j] = {}
                    if code.startswith("Drug"):
                        code = code.split("_")
                        self.drug_map[i + 1][j] = {"Drug": code[0], "Concentration Name": code[1]}
                    elif code.startswith("c"):
                        self.drug_map[i + 1][j] = {"Drug": code, "Concentration Name": "na"}


    def _extract_drug_info_from_drug_conditions(self, drug_info_fn):

        """ Contains the mapping between the drug code and each drug name and concentration"""
        self.drug_info = pd.read_excel(drug_info_fn, sheet_name='Sheet1')

    @staticmethod
    def _column_name_correction(df, old_wave="568", new_wave="555"):

        for k in df.columns:
            if old_wave in k:
                k_new = k.replace(old_wave, new_wave)
                df.rename(columns={k: k_new}, inplace=True)

        return df

    def _load_plate_directories(self, plate_directories, include_phagocytosis=True):

        """
        Load the data for each plate.
        Each directory contains two .txt files: one for cell painting (contains CP in the name) and one
        for phagocytosis (contains Phago in the name).
        Returns a pandas dataframe containing the cell painting and phagocytosis results for each plate, row and column
        """

        dataframes = []
        for n, path in enumerate(plate_directories):
            fns = os.listdir(path)
            for fn in fns:
                if "CP" in fn: # short for Cell Painting
                    fn = os.path.join(path, fn)
                    df_cp = self._extract_info_from_txt(fn)
                    df_cp = self._column_name_correction(df_cp)

                elif include_phagocytosis and "Phago" in fn:
                    fn = os.path.join(path, fn)
                    df_phago = self._extract_info_from_txt(fn)

            # add plate number
            print(path)
            for n in range(99, 0, -1):
                if f"plate {n}" in path.lower():
                    df_cp["Plate"] = len(df_cp) * [f"Plate {n}"]
                    break

            if include_phagocytosis:
                dataframes.append(
                    pd.merge(df_cp, df_phago, on=['Row', 'Column'], how='inner')
                )
            else:
                dataframes.append(df_cp)

        self.df = pd.concat(dataframes, ignore_index=True)

    @staticmethod
    def _extract_info_from_txt(fn):

        with open(fn, "r") as f:
            for n in range(10000):
                data = f.readline()
                if n == 8:
                    columns = data.split("\t")
                    x = {k: [] for k in columns}
                elif n > 8:
                    data = data.split("\t")
                    if len(data) < 2:
                        continue
                    for k, d in zip(columns, data):
                        x[k].append(d)

        return pd.DataFrame(x)

    def _add_drug_info(self):

        """Given the mapping between the drug code and the name and concentration, add the drug name and concentraion
        for each plate, row and column"""

        drug_name = []
        drug_conc = []

        for r, c, plate in zip(self.df["Row"], self.df["Column"], self.df["Plate"]):
            val = self.drug_map[int(r)][int(c)]
            drug_code = val["Drug"]
            conc_name = val["Concentration Name"]

            if drug_code.startswith("Drug"):
                idx = (
                    (self.drug_info["Plate"] == plate) *
                    (self.drug_info["Drug"] == drug_code) *
                    (self.drug_info["Concentration Name"] == conc_name)
                )
                # print(r, c, val)
                drug_name.append(self.drug_info["Drug Name"].values[idx][0])
                drug_conc.append(self.drug_info["Concentration_uM"].values[idx][0])
            else:
                drug_name.append(drug_code)
                drug_conc.append(0.0)

        self.df["Drug"] = drug_name
        self.df["Concentration"] = drug_conc

    def _add_normalized_data(self, baseline_drug = "c_aB_DM"):

        n_samples = self.adata.shape[0]
        for k in self.features_to_normalize:
            vals = np.zeros(n_samples)

            for plate in self.adata.obs.Plate.unique():

                idx = self.adata.obs.Plate == plate
                idx_baseline = (self.adata.obs.Plate == plate) * (self.adata.obs.Drug == baseline_drug)
                baseline = np.mean(self.adata.obs[k].values[idx_baseline])
                vals[idx] = self.adata.obs[k].values[idx] - baseline

            new_k = f"{k} normalized"
            self.adata.obs[new_k] = vals

        X = copy.deepcopy(self.adata.X)
        for plate in self.adata.obs.Plate.unique():
            idx = self.adata.obs.Plate == plate
            idx_baseline = (self.adata.obs.Plate == plate) * (self.adata.obs.Drug == baseline_drug)
            baseline = np.mean(self.adata.X[idx_baseline, :])
            X[idx, :] = self.adata.X[idx, :] - baseline

        self.adata.layers = {"X_norm": X}

    def create_anndata(self, eps = 1e-4, max_val = 5.0):

        idx = []
        features = []
        for n, k in enumerate(self.df.columns):
            if "Nuclei Selected" in k:
                idx.append(n)
                features.append(k)

        x = self.df.to_numpy()
        x = x[:, idx]
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                if isinstance(x[i, j], float):
                    continue
                elif len(x[i, j]) > 0:
                    x[i, j] = float(x[i, j])
                else:
                    x[i, j] = np.nan

        # subset rows
        x = np.float32(x)
        idx_rows = ~np.all(np.isnan(x), axis=1)
        x = x[idx_rows, :]

        # subset columns
        idx = ~np.all(np.isnan(x), axis=0)
        x = x[:, idx]
        features = np.array(features)[idx]

        x -= np.nanmean(x, axis=0, keepdims=True)
        x /= (eps + np.nanstd(x, axis=0, keepdims=True))

        x = np.clip(x, -max_val, max_val)

        df_var = {"feature": features}
        df_obs = {}
        for k in self.observation_features:
            df_obs[k] = np.array(self.df[k].values)[idx_rows]
            if not k in  ["Plate", "Drug"]:
                df_obs[k] = df_obs[k].astype(np.float32)

        df_obs["Phagocytosis per cell"] = df_obs["Non border cells - Number of Spots - Mean per Well"] / (
            1e-3 + df_obs["Non border cells - Number of Objects"]
        )

        self.adata = ad.AnnData(X = x, obs = df_obs, var = df_var)
        self._add_normalized_data()

class Analysis:

    def __init__(self, adata):

        self.adata = adata.copy()
        self._define_features()

    def pls_regression(self, target_name, n_components=10):

        # k = "Non border cells - Number of Spots - Mean per Well"
        pls = PLSRegression(n_components = n_components)
        pls.fit(self.adata.X, self.adata.obs[target_name])
        self.x_pls = pls.transform(self.adata.X)
        self.adata.obsm = {"X_pls": self.x_pls}
        self.adata.uns = {"pls_coef": pls.coef_}

    def plot_all_dose_responses(self, normalized=True, idx_sort = None):

        suffix = " normalized" if normalized else ""
        obs = [
            f"Non border cells - Number of Spots - Mean per Well{suffix}",
            f"Non border cells - Number of Objects{suffix}",
            f"Phagocytosis per cell{suffix}"
        ]

        f, ax = plt.subplots(1, 3, figsize=(10, 7), sharey=True)
        for i, k in enumerate(obs):

            drugs = []
            responses = []
            for drug in self.adata.obs.Drug.unique():
                if isinstance(drug, float) or drug.startswith("c"):
                    continue
                a = self.adata[self.adata.obs.Drug == drug]
                concentrations = np.sort(a.obs.Concentration.unique())

                x = []
                y = []
                for c in concentrations:
                    a0 = a[a.obs.Concentration == c]
                    vals = a0.obs[k].values
                    x.append(1000 * c)
                    y.append(np.mean(vals))

                if len(y) == 8:
                    drugs.append(drug)
                    responses.append(y)

            responses = np.stack(responses, axis=0)
            max_vals = np.max(responses, axis=1)
            if i == 0 and idx_sort is None:
                idx_sort = np.argsort(max_vals)[::-1]

            responses = responses[idx_sort, :]
            drugs = np.array(drugs)[idx_sort]

            ax[i].imshow(responses)
            ax[i].set_yticks(np.arange(len(drugs)), drugs, fontsize=8)
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes('right', size='5%', pad=0.03)
            norm = matplotlib.colors.Normalize(vmin=np.min(responses), vmax=np.max(responses))
            f.colorbar(cm.ScalarMappable(norm=norm), cax=cax, orientation='vertical')
            ax[i].set_xlabel("Concentration")
            ax[i].set_xticks([])
            if "-" in k:
                idx = k.find("-")
                k1 = f"{k[:idx+2]} \n {k[idx+2:]}"
            else:
                k1= k
            ax[i].set_title(k1, fontsize=10)
        plt.show()

        return idx_sort


    def plot_dose_response(self, drug, k = "Non border cells - Number of Spots - Mean per Well"):

        a0 = self.adata[self.adata.obs.Drug == "c_aB_DM"]
        a1 = self.adata[self.adata.obs.Drug == "c_aB"]
        a = self.adata[self.adata.obs.Drug == drug]

        y0 = a0.obs[k].values
        y1 = a1.obs[k].values
        y = a.obs[k].values
        x = 1000 * a.obs["Concentration"].values

        plt.semilogx(len(y0) * [5], y0, 'r.', label="c_aB_DM")
        plt.semilogx(len(y1) * [3], y1, 'g.', label="c_aB")
        plt.semilogx( [5], np.mean(y0), 'ro', markersize=10)
        plt.semilogx([3], np.mean(y1), 'go', markersize=10)
        plt.semilogx(x, y, 'k.', label=drug)
        plt.legend()
        plt.ylabel(k)
        plt.xlabel("Concentration (nM)")
        #plt.xticks(np.arange(), minor=True)  # Minor ticks
        #plt.xticks([10, 100, 1000], minor=False)  # Major ticks

        plt.show()

    def _define_features(self):

        self.features = [
            ["Cell", "Membrane", "Cytoplasm", "Ring", "Nucleus"],
            ["488", "512", "555", "641", "33342"],
            ["Radial Mean", "Length", "Deviation", "Compactness", "Symmetry"],
            ["Dark", "Hole", "Edge", "Filter", "Bright", "Valley", "Saddle", "Ridge", "Spot", "Profile"]
        ]
        self.feature_size = [len(f) for f in self.features]


    def count_features(self, feature_list):

        counts = np.zeros(self.feature_size, dtype=np.int16)
        for feature in feature_list:
            i = [None, None, None, None]
            for j in range(4):
                for n, f in enumerate(self.features[j]):
                    if f in feature:
                        i[j] = n
                        break

            if np.all(np.array([not j is None for j in i])):
                counts[i[0], i[1], i[2], i[3]] += 1
            else:
                pass
                # print(feature)

        return counts

    def plot_feature_counts(self, feature_list):

        counts = self.count_features(feature_list)
        f, ax = plt.subplots(1, 4, figsize=(11, 5))
        axes = [(0, 1), (0, 2), (1, 2), (1, 3)]
        axes_remain = [(2, 3), (1, 3), (0, 3), (0, 2)]
        for n in range(4):
            c = np.sum(counts, axis=axes[n])
            ax[n].imshow(c)
            ax[n].grid(False)
            n0 = len(self.features[axes_remain[n][0]])
            n1 = len(self.features[axes_remain[n][1]])
            ax[n].set_yticks(np.arange(n0), self.features[axes_remain[n][0]], fontsize=8)
            ax[n].set_xticks(np.arange(n1), self.features[axes_remain[n][1]], fontsize=8, rotation=-45, ha="left")
            divider = make_axes_locatable(ax[n])
            cax = divider.append_axes('right', size='3%', pad=0.02)
            norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(c))
            f.colorbar(cm.ScalarMappable(norm=norm), cax=cax, orientation='vertical')

        plt.tight_layout()
        plt.show()

