import numpy as np

# ============================================================
#  Final Unitarity Coefficients K for all dimension-8 operators
#  Based exactly on Eboli et al. Table I, II, III.
#
#  UI = (K / |f|)^(1/4)
#
#  All K values below are the FINAL coefficients after all
#  eigenvalue normalization, symmetry factors and partial-wave
#  definitions in the paper.
# ============================================================

K = {
    # -------------------------
    # FS (Scalar) – Table I
    # -------------------------
    "FS0": 32 * np.pi,
    "FS1": (96/7) * np.pi,
    "FS2": (96/5) * np.pi,

    # -------------------------
    # FM (Mixed) – Table II
    # -------------------------
    "FM0": (32/np.sqrt(6)) * np.pi,
    "FM1": (128/np.sqrt(6)) * np.pi,
    "FM2": (16/np.sqrt(2)) * np.pi,
    "FM3": (64/np.sqrt(2)) * np.pi,
    "FM4": 32 * np.pi,
    "FM5": 64 * np.pi,
    "FM7": (256/np.sqrt(6)) * np.pi,

    # -------------------------
    # FT (Tensor) – Table III
    # -------------------------
    "FT0": (12/5) * np.pi,
    "FT1": (24/5) * np.pi,
    "FT2": (96/13) * np.pi,
    "FT3": (32/3) * np.pi,
    "FT4": 16 * np.pi,
    "FT5": (8/np.sqrt(3)) * np.pi,
    "FT6": (48/7) * np.pi,
    "FT7": (32/np.sqrt(3)) * np.pi,
    "FT8": (3/2) * np.pi,
    "FT9": (24/7) * np.pi,
}


# ============================================================
#  Compute UI
# ============================================================
def compute_UI(operator, f):
    """
    operator : string (e.g. "FT8")
    f        : Wilson coefficient value (TeV^-4)
    """
    if operator not in K:
        raise ValueError(f"Unknown operator: {operator}")

    if f == 0:
        return np.inf

    return (K[operator] / abs(f)) ** 0.25


# ============================================================
#  Example: Your input limits (2σ high)
# ============================================================
limits = {
    'FS0': 36.77,
    'FS1': 28.32,
    'FS2': 36.76,

    'FM0': 6.63,
    'FM1': 22.82,
    'FM2': 7.84,
    'FM3': 26.56,
    'FM4': 18.24,
    'FM5': 26.43,
    'FM7': 40.78,

    'FT0': 1.45,
    'FT1': 1.59,
    'FT2': 3.59,
    'FT3': 3.71,
    'FT4': 8.28,
    'FT5': 3.14,
    'FT6': 3.91,
    'FT7': 10.09,
    'FT8': 2.05,
    'FT9': 4.26,
}

########exp#######
# limits = {
#     'FS0': 40.34,
#     'FS1': 31.60,
#     'FS2': 40.33,

#     'FM0': 6.91,
#     'FM1': 23.98,
#     'FM2': 8.39,
#     'FM3': 28.62,
#     'FM4': 19.53,
#     'FM5': 28.22,
#     'FM7': 43.99,

#     'FT0': 1.49,
#     'FT1': 1.68,
#     'FT2': 3.70,
#     'FT3': 3.81,
#     'FT4': 8.61,
#     'FT5': 3.25,
#     'FT6': 4.17,
#     'FT7': 10.48,
#     'FT8': 2.18,
#     'FT9': 4.56,
# }

# limits = {
#     'FT0': 0.24,
#     'FT1': 0.31,
#     'FT2': 0.63,
#     'FT8': 0.43,
#     'FT9': 0.92,
# }
print(f"{'Operator':<10} | {'f (TeV^-4)':<12} | {'UI (TeV)':<10}")
print("-" * 40)

for op, f in limits.items():
    ui = compute_UI(op, f)
    print(f"{op:<10} | {f:<12.3f} | {ui:<10.3f}")
