import itertools
import multiprocessing as mp

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import wasserstein_distance


def _emd_donorpair_component(
    adata: AnnData,
    donorID_pair,
    obs_donorID: str = "donorID",
    index_comp: int = 1,
    basis: str = "X_diffmap",
):
    """
    Calculates the Earth mover's distance (Wasserstein distance) between two donors for a single component.

    Returns a single EMD.
    """
    donorID_1 = donorID_pair[0]
    donorID_2 = donorID_pair[1]

    import warnings

    if (basis == "X_diffmap") & (index_comp == 0):
        warnings.warn("Diffusion component 0 is not meaningful and should be ignored.", UserWarning, stacklevel=2)

    df_1 = adata.obsm[basis][adata.obs[obs_donorID] == donorID_1, index_comp]
    df_2 = adata.obsm[basis][adata.obs[obs_donorID] == donorID_2, index_comp]

    dist_EMD = wasserstein_distance(df_1, df_2)

    return dist_EMD


def _emd_donorpair_pmt(
    adata: AnnData,
    donorID_pair,
    obs_donorID: str = "donorID",
    basis: str = "X_diffmap",
    n_comp: int = 5,
    n_pmt_pval: int = 1000,
):
    """
    Calculates the Earth mover's distance (Wasserstein distance) between two donors for multiple permutations of the donorID assignment and components.

    Returns an array of EMDs.
    """
    filter_sub = [pID in donorID_pair for pID in adata.obs[obs_donorID]]
    adata_sub = adata[
        filter_sub, :
    ].copy()  # TODO do not copy in the future for more performance - needs different vector passing

    if basis == "X_diffmap":
        index_comps = range(1, n_comp + 1)
    else:
        index_comps = range(n_comp)

    array_donor_EMD = np.zeros([n_pmt_pval + 1, n_comp], dtype=np.float64)
    print(donorID_pair)
    # permuted EMD
    for seed in range(n_pmt_pval):
        rng = np.random.default_rng(seed=seed)
        adata_sub.obs["donorArray_pmt"] = rng.permutation(adata_sub.obs[obs_donorID])
        for j in range(n_comp):
            array_donor_EMD[seed, j] = _emd_donorpair_component(
                adata_sub, donorID_pair, obs_donorID="donorArray_pmt", index_comp=index_comps[j], basis="X_diffmap"
            )
    # actual EMD
    for j in range(n_comp):
        array_donor_EMD[seed + 1, j] = _emd_donorpair_component(
            adata_sub, donorID_pair, obs_donorID=obs_donorID, index_comp=index_comps[j], basis="X_diffmap"
        )

    del adata_sub
    return array_donor_EMD


def _emd_donorpair_pmt_wrapper(x):
    array_donor_EMD = _emd_donorpair_pmt(x[0], x[1], n_comp=x[2], n_pmt_pval=x[3])
    return array_donor_EMD


def emd_pval(
    adata: AnnData,
    obs_donorID: str = "donorID",
    basis: str = "X_diffmap",
    n_comp: int = 5,
    n_pmt_pval: int = 1000,
    n_cpu: int = 4,
):
    """
    Calculates the Earth mover's distance (Wasserstein distance) and calculates pvalues using a permutation test.

    Returns a pandas DataFrame.
    """
    # donor pairs
    donorIDs = adata.obs[obs_donorID].unique()
    paired_donorIDs = list(itertools.combinations(donorIDs, 2))

    pool = mp.Pool(n_cpu, maxtasksperchild=1)
    results = pool.map(
        _emd_donorpair_pmt_wrapper, [(adata, paired_donorID, n_comp, n_pmt_pval) for paired_donorID in paired_donorIDs]
    )
    pool.close()

    array_emd = np.array(results)
    array_sim_emd = array_emd[:, :-1, :]
    array_act_emd = array_emd[:, -1, :]

    pval_pmt_emd = np.zeros([len(paired_donorIDs), n_comp], dtype=np.float64)
    for i in range(len(paired_donorIDs)):
        for j in range(n_comp):
            pval_pmt_emd[i, j] = sum(array_sim_emd[i, :, j] > array_act_emd[i, j]) / array_sim_emd.shape[1]

    paired_donorIDs_1, paired_donorIDs_2 = zip(*paired_donorIDs, strict=False)
    dict_emd = {
        "donorID_1": np.tile(np.asarray(paired_donorIDs_1), n_comp),
        "donorID_2": np.tile(np.asarray(paired_donorIDs_2), n_comp),
        "basis": np.tile(basis, n_comp * len(paired_donorIDs)),
        "component": np.tile(np.arange(n_comp), len(paired_donorIDs)),
        "emd": array_act_emd.flatten(order="F"),
        "pval": pval_pmt_emd.flatten(order="F"),
    }
    if basis == "X_diffmap":
        # adding 1 as the first diffusion component is ignored
        dict_emd["component"] += 1
    # NOTE: consider moving over to polars in the future
    df_emd = pd.DataFrame(data=dict_emd)

    return df_emd


# def emd(
#         adata: AnnData,
#         obs_donorID: str = 'donorID',
#         basis: str = 'X_diffmap',
#         n_comp: int=5,
#         n_cpu: int=1,
#         ):
#     """
#     Calculates the Earth mover's distance (Wasserstein distance).
#     """

#     df_emd = emd_pval(adata,
#                       obs_donorID = obs_donorID,
#                       basis = basis,
#                       n_comp = n_comp,
#                       n_pmt_pval = 0,
#                       n_cpu = n_cpu,
#         )

#     return df_emd
