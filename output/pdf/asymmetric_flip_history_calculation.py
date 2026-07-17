#!/usr/bin/env python3
"""Scratch calculation for repeated asymmetric edge-flip history transfer.

This script deliberately separates the calculation from the explanatory
document.  It constructs a deformed four-node square, computes the affine
deformation gradient of each triangle before and after the diagonal flip,
iterates the two history rules, and evaluates the final reconstructed state
with the Conti square potential used by MTS2D.
"""

from __future__ import annotations

import numpy as np


EPSILON = 0.05
REPETITIONS = 5
BETA = -0.25
BULK_MODULUS = 4.0
NOISE = 1.0


def edge_matrix(points: dict[str, np.ndarray], triangle: tuple[str, str, str]) -> np.ndarray:
    x0, x1, x2 = (points[name] for name in triangle)
    return np.column_stack((x1 - x0, x2 - x0))


def deformation_gradient(
    current: dict[str, np.ndarray],
    reference: dict[str, np.ndarray],
    triangle: tuple[str, str, str],
) -> np.ndarray:
    return edge_matrix(current, triangle) @ np.linalg.inv(edge_matrix(reference, triangle))


def elastic_domain_margin(F: np.ndarray) -> float:
    """Positive means the square metric is inside the unreduced elastic domain."""
    C = F.T @ F
    return 0.5 * min(C[0, 0], C[1, 1]) - abs(C[0, 1])


def conti_square_energy(F: np.ndarray) -> float:
    """MTS2D Conti square energy density, before subtracting the ground state."""
    C = F.T @ F
    C11, C22, C12 = C[0, 0], C[1, 1], C[0, 1]
    x0 = C11 * C22 - C12**2
    x1 = NOISE * x0
    x2 = x0 ** (-0.5)
    x3 = C11 * x2
    x4 = C22 * x2
    x5 = (x3 - x4) ** 2
    x6 = C12 * x2
    x7 = x3 + x4 - 4.0 * x6
    x8 = 0.25 * x5 + (1.0 / 12.0) * x7**2
    x9 = x8**3
    x10 = x5 * x7 - (1.0 / 9.0) * x7**3
    x11 = x10**2
    x12 = (x3 + x4 - x6) / 3.0
    x13 = x10 * x8
    return (
        BULK_MODULUS * (x1 - np.log(x1))
        + BETA
        * (
            (1.0 / 1056.0) * x11
            + x12**4 * x8
            + (7.0 / 198.0) * x13 * (x3 + x4 - x6)
            - (41.0 / 99.0) * x9
        )
        + x10 * x12**3
        + (17.0 / 528.0) * x11
        - (8.0 / 33.0) * x13 * (x3 + x4 - x6)
        + (4.0 / 11.0) * x9
    )


def first_piola_finite_difference(F: np.ndarray, step: float = 1.0e-6) -> np.ndarray:
    P = np.empty((2, 2), dtype=float)
    for row in range(2):
        for col in range(2):
            perturbation = np.zeros((2, 2), dtype=float)
            perturbation[row, col] = step
            P[row, col] = (
                conti_square_energy(F + perturbation)
                - conti_square_energy(F - perturbation)
            ) / (2.0 * step)
    return P


def format_matrix(A: np.ndarray) -> str:
    return np.array2string(A, precision=9, suppress_small=False)


def main() -> None:
    e = EPSILON
    reference = {
        "i": np.array([0.0, 0.0]),
        "j": np.array([1.0, 0.0]),
        "k": np.array([1.0, 1.0]),
        "b": np.array([0.0, 1.0]),
    }
    current = {
        "i": np.array([0.0, 0.0]),
        "j": np.array([1.0 + e, 0.0]),
        "k": np.array([1.0, 1.0]),
        "b": np.array([0.0, 1.0 + e]),
    }

    # Counter-clockwise triangles.  Each descendant keeps the exterior edge
    # it shared with its parent.
    old_large_parent = ("i", "j", "k")
    old_small_parent = ("i", "k", "b")
    new_large = ("i", "j", "b")
    new_small = ("j", "k", "b")

    F_old_large = deformation_gradient(current, reference, old_large_parent)
    F_old_small = deformation_gradient(current, reference, old_small_parent)
    F_new_large = deformation_gradient(current, reference, new_large)
    F_new_small = deformation_gradient(current, reference, new_small)

    transfer_large = np.linalg.inv(F_new_large) @ F_old_large
    transfer_small = np.linalg.inv(F_new_small) @ F_old_small

    H_mathcal_large = np.linalg.matrix_power(transfer_large, REPETITIONS)
    H_mathcal_small = np.linalg.matrix_power(transfer_small, REPETITIONS)

    # All four local deformation gradients remain in the identity lattice
    # well, so F_P(old) = F_P(new) = I and the plastic-history update is I.
    H_plastic_large = np.eye(2)
    H_plastic_small = np.eye(2)

    ground_energy = conti_square_energy(np.eye(2))

    print(f"epsilon={EPSILON:.6f}, repetitions={REPETITIONS}")
    print("\nDeformation gradients before and after one flip")
    for name, F in (
        ("F_old_large_parent", F_old_large),
        ("F_old_small_parent", F_old_small),
        ("F_new_large", F_new_large),
        ("F_new_small", F_new_small),
    ):
        print(f"{name}=\n{format_matrix(F)}")
        print(
            f"  det={np.linalg.det(F):.9f}, "
            f"elastic-domain margin={elastic_domain_margin(F):.9f}"
        )

    print("\nCurrent triangle areas")
    for name, triangle in (
        ("old_large_parent", old_large_parent),
        ("old_small_parent", old_small_parent),
        ("new_large", new_large),
        ("new_small", new_small),
    ):
        area = 0.5 * abs(np.linalg.det(edge_matrix(current, triangle)))
        print(f"{name}: {area:.9f}")

    print("\nOne-flip direct-history transfer matrices")
    print(f"A_large=\n{format_matrix(transfer_large)}")
    print(f"A_small=\n{format_matrix(transfer_small)}")

    print(f"\nHistories after {REPETITIONS} identical asymmetric transfers")
    print(f"H_mathcal_large=\n{format_matrix(H_mathcal_large)}")
    print(f"H_mathcal_small=\n{format_matrix(H_mathcal_small)}")
    print(f"H_plastic_large=\n{format_matrix(H_plastic_large)}")
    print(f"H_plastic_small=\n{format_matrix(H_plastic_small)}")

    print("\nReconstructed state after the current element returns to F=I")
    for name, reconstructed in (
        ("mathcal_large", H_mathcal_large),
        ("mathcal_small", H_mathcal_small),
        ("T_total_large", H_plastic_large),
        ("T_total_small", H_plastic_small),
    ):
        C = reconstructed.T @ reconstructed
        excess_energy = conti_square_energy(reconstructed) - ground_energy
        P = first_piola_finite_difference(reconstructed)
        print(f"{name}=\n{format_matrix(reconstructed)}")
        print(f"  C=\n{format_matrix(C)}")
        print(
            f"  det={np.linalg.det(reconstructed):.9f}, "
            f"elastic-domain margin={elastic_domain_margin(reconstructed):.9f}, "
            f"excess energy={excess_energy:.9f}, ||P||_F={np.linalg.norm(P):.9f}"
        )

    # Exact reverse flips are a control: inv(A)^N A^N = I.
    reverse_control_large = np.linalg.matrix_power(
        np.linalg.inv(transfer_large), REPETITIONS
    ) @ H_mathcal_large
    reverse_control_small = np.linalg.matrix_power(
        np.linalg.inv(transfer_small), REPETITIONS
    ) @ H_mathcal_small
    print("\nExact reverse-flip control")
    print(f"large residual norm={np.linalg.norm(reverse_control_large - np.eye(2)):.3e}")
    print(f"small residual norm={np.linalg.norm(reverse_control_small - np.eye(2)):.3e}")


if __name__ == "__main__":
    main()
