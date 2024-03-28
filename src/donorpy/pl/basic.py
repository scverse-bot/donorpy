import pandas as pd
import seaborn as sns
from anndata import AnnData


def basic_plot(adata: AnnData) -> int:
    """Generate a basic plot for an AnnData object.

    Parameters
    ----------
    adata
        The AnnData object to preprocess.

    Returns
    -------
    Some integer value.
    """
    print("Import matplotlib and implement a plotting function here.")
    return 0


class BasicClass:
    """A basic class.

    Parameters
    ----------
    adata
        The AnnData object to preprocess.
    """

    my_attribute: str = "Some attribute."
    my_other_attribute: int = 0

    def __init__(self, adata: AnnData):
        print("Implement a class here.")

    def my_method(self, param: int) -> int:
        """A basic method.

        Parameters
        ----------
        param
            A parameter.

        Returns
        -------
        Some integer value.
        """
        print("Implement a method here.")
        return 0

    def my_other_method(self, param: str) -> str:
        """Another basic method.

        Parameters
        ----------
        param
            A parameter.

        Returns
        -------
        Some integer value.
        """
        print("Implement a method here.")
        return ""


def donor_emd_heatmap(
    df_emd: pd.core.frame.DataFrame,
    compID: int = 1,
    pval_cut: float = None,
):
    """
    Returns a heatmap of donor emd values for a single component

    Parameters
    ----------
    df_emd
        The pandas DataFrame containing the EMD

    Returns
    -------
    Seaborn plot.
    """
    df_emd_comp = df_emd.query("component == @compID")
    if pval_cut:
        # set uncertain distances to 0
        df_emd_comp.loc[df_emd_comp["pval"] > 0.1, "emd"] = 0

    df_emd_comp_flip = df_emd_comp.rename(columns={"donorID_1": "donorID_2", "donorID_2": "donorID_1"})
    df_emd_comp_long = pd.concat([df_emd_comp, df_emd_comp_flip], ignore_index=True)
    df_emd_comp_wide = df_emd_comp_long.pivot(index="donorID_1", columns="donorID_2", values="emd").fillna(0)

    g = sns.clustermap(df_emd_comp_wide, metric="euclidean", cmap="mako")
    g.fig.suptitle(df_emd.basis[0] + " " + str(compID))
    ax = g.ax_heatmap
    # ax.set_title(df_emd.basis[0] + " " + str(compID))
    ax.set(xlabel=None)
    ax.set(ylabel=None)
    ax.tick_params(axis="x", rotation=90)
    ax.tick_params(axis="y", rotation=0)

    return g
