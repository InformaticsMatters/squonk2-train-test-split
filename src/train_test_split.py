#!/usr/bin/env python
"""Template echo utility."""

import argparse
import logging
import os
from collections import defaultdict

import numpy as np
import pandas as pd

# Import Data Manager DmLog utility.
# Messages emitted using this result in Task Events.
from dm_job_utilities.dm_log import DmLog
from dm_job_utilities.utils import read_delimiter
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import rdkit_utils

# import deepchem as dc
# from deepchem.data import NumpyDataset
# from deepchem.splits import ScaffoldSplitter


DEFAULT_SPLIT_RATIOS = [0.5, 0.25, 0.25]
DEFAULT_SPLIT_NAMES = ["training", "test", "validation"]

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_scaffold(mol: Chem.rdchem.Mol) -> str:
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def group_by_scaffold(mol_list):
    scaffold_to_indices = defaultdict(list)

    for idx, mol in enumerate(mol_list):
        scaffold = get_scaffold(mol)
        if scaffold is not None:
            scaffold_to_indices[scaffold].append(idx)

    return scaffold_to_indices


def weighted_scaffold_split(mol_list, split_fractions, seed=42):
    rng = np.random.default_rng(seed)

    scaffold_groups = group_by_scaffold(mol_list)

    # sort scaffold groups largest first
    groups = list(scaffold_groups.values())
    rng.shuffle(groups)
    groups.sort(key=len, reverse=True)

    split_names = list(split_fractions.keys())
    n_total = len(mol_list)

    target_sizes = {
        name: int(round(frac * n_total)) for name, frac in split_fractions.items()
    }

    splits = {name: [] for name in split_names}
    current_sizes = {name: 0 for name in split_names}

    for group in groups:
        deficits = {
            name: target_sizes[name] - current_sizes[name] for name in split_names
        }

        viable = [name for name in split_names if deficits[name] >= len(group)]

        if viable:
            max_deficit = max(deficits[n] for n in viable)
            candidates = [n for n in viable if deficits[n] == max_deficit]
        else:
            min_size = min(current_sizes.values())
            candidates = [n for n in split_names if current_sizes[n] == min_size]

        target = rng.choice(candidates)

        splits[target].extend(group)
        current_sizes[target] += len(group)

    # TODO: error handling
    check_scaffold_leakage(mol_list, splits)

    return splits


def weighted_random_split(smiles_list, split_fractions, seed=42):
    """
    Returns indices for each split.
    """

    rng = np.random.default_rng(seed)
    n_samples = len(smiles_list)

    indices = np.arange(n_samples)
    rng.shuffle(indices)

    splits = {}
    start = 0

    for name, frac in split_fractions.items():
        size = int(round(frac * n_samples))
        splits[name] = indices[start : start + size].tolist()
        start += size

    # assign leftovers (rounding)
    leftovers = indices[start:]
    for i, idx in enumerate(leftovers):
        splits[list(splits.keys())[i % len(splits)]].append(idx)

    return splits


def check_scaffold_leakage(mols, splits):
    split_scaffolds = []

    for idxs in splits:
        s = set(get_scaffold(mols[i]) for i in idxs)
        split_scaffolds.append(s)

    for i in range(len(splits)):
        for j in range(i + 1, len(splits)):
            overlap = split_scaffolds[i] & split_scaffolds[j]
            assert not overlap, f"Leakage between split {i} and {j}"


def run(
    filename,
    delimiter=None,
    id_column=None,
    mol_column=None,
    y_column=None,
    omit_fields=False,
    read_header=False,
    write_header=False,
    fragment_method="hac",
    # missing_val=None,
    split_method="random",
    n_splits=None,
    split_ratios=(),
    split_names=(),
):

    DmLog.emit_event("Splitter job started")

    SPLIT_METHODS = {
        "random": weighted_random_split,
        "scaffold": weighted_scaffold_split,
    }

    df = pd.read_csv(
        filename,
        delimiter=delimiter,
        header=0 if read_header else None,
    )

    # convert to column names
    if isinstance(id_column, int):
        id_column = df.columns[id_column]

    if isinstance(mol_column, int):
        mol_column = df.columns[mol_column]

    if isinstance(y_column, int):
        y_column = df.columns[y_column]

    columns = list(set([id_column, mol_column, y_column]))

    # need mol in several places
    df["rdkit_mol"] = df[mol_column].apply(
        lambda smiles: Chem.MolFromSmiles(smiles),  # pylint: disable=unnecessary-lambda
    )
    df["fragment"] = df["rdkit_mol"].apply(rdkit_utils.fragment, args=(fragment_method,))

    seed = 42

    if len(split_names) != n_splits:
        width = len(str(n_splits))
        split_names = [f"group_{i:0{width}d}" for i in range(1, n_splits + 1)]

    split_groups = {split_names[i]: k for i, k in enumerate(split_ratios)}

    method = SPLIT_METHODS.get(split_method, weighted_random_split)

    splits = method(df["rdkit_mol"], split_fractions=split_groups, seed=seed)
    df = df.drop(["rdkit_mol", "fragment"], axis=1)

    # write groups to files
    DmLog.emit_event("Split finished, writing out set files")
    for k, v in splits.items():
        fname = f"{k}.smi"

        if omit_fields:
            df = df.loc[:, columns]

        df.iloc[v].to_csv(fname, encoding="utf-8", index=False, header=write_header)
        os.chmod(fname, 0o664)


def list_of_strings(arg):
    return arg.split(",")


def list_of_floats(arg):
    l = list(map(float, arg.split(",")))
    if not abs(sum(l) - 1.0) < 1e-6:
        DmLog.emit_event("The sum of splits must be equal to 1")
        raise argparse.ArgumentError(arg, "The sum of splits must be equal to 1")
    return l



def main():
    parser = argparse.ArgumentParser(description="Split dataset")
    # to pass tab as the delimiter specify it as $'\t' or use one of
    # the symbolic names 'comma', 'tab', 'space' or 'pipe'
    rdkit_utils.add_common_molecule_io_args(parser, include_y_column=True)

    rdkit_generic_group = parser.add_argument_group("General RDKit options")
    rdkit_generic_group.add_argument(
        "--fragment-method",
        choices=["hac", "mw", "none"],
        default="hac",
        help="Strategy for picking largest fragment (mw or hac or none",
    )
    split_group = parser.add_argument_group("Splitting options")
    split_group.add_argument(
        "--split-method",
        choices=["random", "scaffold"],
        default="random",
    )
    split_number_group = parser.add_mutually_exclusive_group()
    split_number_group.add_argument(
        "--n-splits",
        type=int,
        default=argparse.SUPPRESS,
        help="Split the input set to number of sets of equal size",
    )
    split_number_group.add_argument(
        "--split-ratios",
        type=list_of_floats,
        default=argparse.SUPPRESS,
        help="Split input to the number of sets using the given ratios",
    )
    split_group.add_argument(
        "--split-names",
        type=list_of_strings,
        default=argparse.SUPPRESS,
    )

    args = parser.parse_args()
    delimiter = read_delimiter(args.delimiter)

    # split_ratios_provided = hasattr(args, "split_ratios")
    # n_splits_provided = hasattr(args, "n_splits")
    # split_names_provided = hasattr(args, "split_names")

    # in current setting, allow only 3 sets, training, test and validation
    split_ratios_provided = False
    n_splits_provided = False
    split_names_provided = False

    # if given, override n-splits
    if split_ratios_provided:
        split_ratios = args.split_ratios
        n_splits = len(split_ratios)
        if split_names_provided and len(args.split_names) == n_splits:
            split_names = args.split_names

    elif n_splits_provided:
        n_splits = args.n_splits
        split_ratios = [1.0 / n_splits] * n_splits
        if split_names_provided and len(args.split_names) == n_splits:
            split_names = args.split_names

    else:
        split_ratios = DEFAULT_SPLIT_RATIOS
        split_names = DEFAULT_SPLIT_NAMES
        n_splits = len(split_ratios)

    # if split_names_provided and len(args.split_names) == n_splits:
    #     split_names = args.split_names
    # else:
    #     split_names = []

    print(args)

    run(
        args.input,
        omit_fields=args.omit_fields,
        delimiter=delimiter,
        id_column=args.id_column,
        mol_column=args.mol_column,
        y_column=args.y_column,
        read_header=args.read_header,
        write_header=args.write_header,
        fragment_method=args.fragment_method,
        split_method=args.split_method,
        n_splits=n_splits,
        split_ratios=split_ratios,
        split_names=split_names,
    )


if __name__ == "__main__":
    # python src/train_test_split.py --input=data/caco2_mordred_filtered_scaled.smi --delimiter=comma --id-column=0 --mol-column=Drug --y-column=Y --read-header --write-header --split-method=random

    main()
