"""
Helper module to read herm_matrix_ness Green's functions stored in HDF5 format
(by the NESS code base).

Author: Martin Eckstein 2025
"""

import h5py
import numpy as np

# ----------------------------------------------------------------------

class Dummy:
    """Simple object to store attributes via dot notation."""
    pass

# ----------------------------------------------------------------------

def isgroup(obj):
    return isinstance(obj, h5py.Group)

def isdataset(obj):
    return isinstance(obj, h5py.Dataset)

# ----------------------------------------------------------------------

def read_herm_matrix_ness_h5group(group):
    """
    Reads a Green's function object (herm_matrix_ness) from a given HDF5 group.
    Expects the group to contain:
        - 'Nft' and 'size1' as scalar datasets
        - 'ret' and 'les' groups with 'time' and 'freq' datasets
    Returns a Dummy object with fields:
        - ret_t, ret_w, les_t, les_w (arrays)
        - Nft, size1 (ints)
    """
    d = Dummy()
    d.Nft = int(group["Nft"][()].item())
    d.size1 = int(group["size1"][()].item())

    for component in ["ret", "les"]:
        subgroup = group[component]
        time = subgroup["time"][:]
        freq = subgroup["freq"][:]

        # Reshape depending on matrix size
        if d.size1 == 1:
            time = time[:, 0, 0]
            freq = freq[:, 0, 0]

        setattr(d, f"{component}_t", time)
        setattr(d, f"{component}_w", freq)

    return d

# ----------------------------------------------------------------------

def read_herm_matrix_ness_h5file(filename, key=None):
    """
    Reads a Green's function from a full HDF5 file (not just a group).
    If `key` is provided, reads from that group inside the file.
    Returns a Dummy object with fields like 'ret_t', 'les_t', etc.
    """
    with h5py.File(filename, "r") as f:
        if key is not None:
            return read_herm_matrix_ness_h5group(f[key])
        else:
            if len(f.keys()) != 1:
                raise RuntimeError("File has multiple top-level groups, specify key=")
            root_key = list(f.keys())[0]
            return read_herm_matrix_ness_h5group(f[root_key])

import numpy as np


def read_herm_matrix_ness_text(filename):
    """
    Reads a herm_matrix_ness from a text file (produced by `print_to_file`)
    and returns a Dummy object with attributes:
      - ret_t
      - les_t
      - ret_w
      - les_w
      - Nft
      - size1
    """

    with open(filename, 'r') as f:
        header = f.readline()
        if not header.startswith("#"):
            raise ValueError("Invalid file format. Missing header.")
        
        nft, size1 = map(int, header[1:].strip().split())

        d = Dummy()
        d.Nft = nft
        d.size1 = size1

        def allocate():
            if size1 == 1:
                return np.zeros(nft, dtype=np.complex128)
            else:
                return np.zeros((nft, size1, size1), dtype=np.complex128)

        for label in ["ret_t", "les_t", "ret_w", "les_w"]:
            setattr(d, label, allocate())

        for line in f:
            if ':' not in line:
                continue
            label_part, *data = line.strip().split()
            label = label_part.rstrip(':')
            t_index = int(data[0])
            values = list(map(float, data[1:]))

            if len(values) != 2 * size1 * size1:
                raise ValueError(f"Unexpected number of values for {label} at t={t_index}")

            complex_values = np.array(values[::2]) + 1j * np.array(values[1::2])
            if size1 == 1:
                getattr(d, label)[t_index] = complex_values[0]
            else:
                getattr(d, label)[t_index] = complex_values.reshape((size1, size1))

    return d
