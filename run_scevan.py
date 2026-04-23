# running scevan
# contains EnsDb.Hsapiens.v86 which correlates to hg 38

# imports
# ----------------------------------------------------------------------------------------------------------------------
import rpy2.robjects as robjects
from pathlib import Path
import rdata
import shutil
import pandas as pd
from pyomics.utils import benchmark_method
import itertools
# ----------------------------------------------------------------------------------------------------------------------


def grid_by_dict(pars_dict: dict) -> list:
    keys=pars_dict.keys()
    combinations = itertools.product(*pars_dict.values())
    list_of_kwargs = [dict(zip(keys, cc)) for cc in combinations]
    return list_of_kwargs


def val_build_project() -> (Path, Path):
    cwd_path = Path.cwd()
    print(f"Current working directory of running script {Path(__file__).name}: {cwd_path}")
    path_out = cwd_path / "app" / "out"
    path_in = cwd_path / "data_input"

    if not path_in.exists():
        raise ValueError(f"Data dir '{str(path_in)}' does not exist!")

    if not path_out.exists():
        path_out.mkdir(parents=True, exist_ok=True)
        print(f"Data out-dir '{str(path_out)}' has been created...")

    return path_in, path_out


def get_hg_38_file_paths(target_path: Path) -> list:
    return list(target_path.rglob("*__hg_38__RCM.csv"))


def get_hg_38_desc_paths(target_path: Path) -> dict:
    return {p.stem: p for p in target_path.rglob("*__hg_38.txt")}


# add decorator for performance
def run_scevan(path_target: Path,
                path_out_data: Path,
                n_cores: int = 10,
                n_genes_chr: int = 5,
                perc_genes: int = 10,
                beta_vega: float = 0.5,
                cell_pre_label: bool = False):
    list_paths_target_csvs = get_hg_38_file_paths(path_target)
    if cell_pre_label:
        dict_paths_target_txts = get_hg_38_desc_paths(path_target)
        list_paths_target_csvs = [p for p in list_paths_target_csvs if p.stem.split("__RCM")[0] in dict_paths_target_txts.keys()]
    for p in list_paths_target_csvs:
        name_tag = f"{p.stem}__n,{n_cores};n,{n_genes_chr};p,{perc_genes};b,{beta_vega};c,{cell_pre_label}__scevan_"
        path_out_target = path_out_data / f"out__{name_tag}"
        path_out_target.mkdir(parents=True, exist_ok=True)

        norm_cell_vector = robjects.NULL
        if cell_pre_label:
            path_txt = dict_paths_target_txts[p.stem.split("__RCM")[0]]   # normal cells split by \n
            with open(path_txt, "r") as f:
                list_norm_cells = list(map(lambda x: x.replace("\n", ""), f.readlines()))
                norm_cell_vector = robjects.vectors.StrVector(list_norm_cells)

        @benchmark_method(str(path_out_target))
        def run_rscript(p,
                        name_tag,
                        n_cores,
                        norm_cell_vector,
                        n_genes_chr,
                        perc_genes,
                        beta_vega):
            r = robjects.r
            r.source("c_scevanR.R")
            r.r_run_scevan(str(p), name_tag, n_cores, norm_cell_vector, n_genes_chr, perc_genes, beta_vega)

        run_rscript(p,
                    name_tag,
                    n_cores,
                    norm_cell_vector,
                    n_genes_chr,
                    perc_genes,
                    beta_vega)

        # generate a csv-matrix for the CNA-output
        path_cnv_rdata = [p for p in Path("./output").glob("*__scevan__CNAmtx.RData")][0]
        df_cnv = rdata.read_rds(path_cnv_rdata)["CNA_mtx_relat"].to_pandas()
        path_annot_rdata = [p for p in Path("./output").glob("*__scevan__count_mtx_annot.RData")][0]
        df_pos = rdata.read_rds(path_annot_rdata)['count_mtx_annot'].set_index("gene_name").drop("gene_id", axis=1).rename({"seqnames":"CHR", "start":"START", "end":"END"}, axis=1).astype(int)
        df_concat = pd.concat([df_pos, df_cnv], axis=1).set_index("CHR")
        df_concat.to_csv(Path.cwd() / f"{name_tag}_GBC.csv")

        list_all_nametag_items = [p for p in Path.cwd().glob(f"*{name_tag}*")]
        for items in list_all_nametag_items:
            shutil.move(items, path_out_target / items.name)
        shutil.move("output", path_out_target / "output")


if __name__ == "__main__":
    # matrix of possible scevan hyperparameter kwargs
    # 0 for beta_vega not allowed (comparable to copykat's KS.cut)
    kwargs_gridsearch = {
        "n_cores": [5, 10, 20],
        "n_genes_chr": [1, 5, 10, 100],
        "perc_genes": [0, 5, 10, 20, 30],
        "beta_vega": [0.1, 0.5, 1, 2, 3, 4],
        "cell_pre_label": [True, False]
    }

    path_in, path_out = val_build_project()
    run_scevan(path_in, path_out, n_cores=1, cell_pre_label=False)
    run_scevan(path_in, path_out, n_cores=1, cell_pre_label=True)  # standard parameters; one core
    list_kwargs = grid_by_dict(kwargs_gridsearch)

    for kwarg_opt in list_kwargs:
        print(f"SCEVAN running with hyperparameters: {kwarg_opt}")
        run_scevan(path_in, path_out, **kwarg_opt)



