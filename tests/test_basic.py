import anndata as ad
import numpy as np
import pandas as pd
import pytest

import donorpy as dn
from donorpy.core import _emd_donorpair_pmt_wrapper


def test_package_has_version():
    assert dn.__version__ is not None


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example():
    assert 1 == 0  # This test is designed to fail.


def test_emd_identical():
    """Tests that EMD is 0 when comparing two identical donors."""

    n_cells_donor = 50
    n_comp = 3
    metadata = pd.DataFrame(
        data={"donorID": np.repeat(["donor_1", "donor_2"], n_cells_donor)},
        index=[str(x) for x in range(n_cells_donor * 2)],
    )
    rng = np.random.default_rng()
    dist_1 = rng.random((n_cells_donor, n_comp))
    dist_2 = dist_1
    X_diffmap = np.concatenate((dist_1, dist_2), axis=0)

    adata = ad.AnnData(obs=metadata, obsm={"X_diffmap": X_diffmap})

    df_emd = dn.emd_pval(
        adata,
        obs_donorID="donorID",
        basis="X_diffmap",
        n_comp=2,
        n_pmt_pval=5,
        n_cpu=1,
    )
    assert all(df_emd.emd == 0.0)
    assert all(df_emd.pval == 1.0)


def test_emd_donorpair_pmt():
    """Tests that EMD is 0 when comparing two identical donors. Tests function that is parallelised."""

    n_cells_donor = 50
    n_comp = 3
    metadata = pd.DataFrame(
        data={"donorID": np.repeat(["donor_1", "donor_2"], n_cells_donor)},
        index=[str(x) for x in range(n_cells_donor * 2)],
    )
    rng = np.random.default_rng()
    dist_1 = rng.random((n_cells_donor, n_comp))
    dist_2 = dist_1
    X_diffmap = np.concatenate((dist_1, dist_2), axis=0)

    adata = ad.AnnData(obs=metadata, obsm={"X_diffmap": X_diffmap})

    array_donor_EMD = _emd_donorpair_pmt_wrapper(
        [
            adata,
            ["donor_1", "donor_2"],
            2,
            5,
        ]
    )

    assert all(array_donor_EMD[-1, :] == 0.0)

    # Test EMD >= 1.0 for non overlapping donors
    dist_2 = dist_1 + 2.0
    X_diffmap = np.concatenate((dist_1, dist_2), axis=0)

    adata = ad.AnnData(obs=metadata, obsm={"X_diffmap": X_diffmap})

    array_donor_EMD = _emd_donorpair_pmt_wrapper(
        [
            adata,
            ["donor_1", "donor_2"],
            2,
            5,
        ]
    )

    assert all(array_donor_EMD[-1, :] >= 1.0)
