import numpy as np

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    valid_bases = {"A", "C", "G", "T"}

    compared_columns = 0
    mismatches = 0

    for ref_base, other_base in zip(reference_sequence.upper(), distant_sequence.upper()):
        if ref_base not in valid_bases or other_base not in valid_bases:
            continue

        compared_columns += 1

        if ref_base != other_base:
            mismatches += 1

    if compared_columns == 0:
        return 0.0

    prop = mismatches / compared_columns

    inside_log = 1 - (4 / 3) *prop 
    if inside_log <= 0:
        return float("inf")

    d_per_site = -0.75 * np.log(inside_log)

    return float(d_per_site * compared_columns)


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
