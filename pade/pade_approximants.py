#!/usr/bin/env python3
"""
Improved Padé Approximations for kh_numeric over k₀h ∈ [0, 2π] with Domain Rescaling
and Automatic Generation of a Routines File.

This program does the following:
  • Computes Padé approximants for the function defined by
         k₀h = kh * tanh(kh)
    via Newton–Raphson, and with the transformation u = √(k₀h).
    
  • Approximates g(u)=kh_numeric(u²) by a rational function of the form:
         R(u) = u * (p₀ + p₁·v² + ... + p_M·v^(2M)) / (1 + q₁·v² + ... + q_M·v^(2M))
    where v = u/√(2π) and used_degree = 2M+1 (if an even degree is requested, it is reduced by one).
    
  • Writes a detailed report to "pade_output.txt" and plots the approximants.
  
  • Also outputs a file "pade_routines.txt" containing a function
         def pade2025(k0h, formula):
    which has an if/elif chain; each branch hard-codes one approximant as a Python expression.
  
See the code below.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# -------------------------------
# Exact dispersion function (kh_numeric)
# -------------------------------
def kh_numeric(k0h, tol=1e-20, max_iter=1000):
    """
    Compute the 'exact' nondimensional wavenumber, kh, by solving
         k₀h = kh * tanh(kh)
    using Newton–Raphson iteration.
    """
    if k0h == 0:
        return 0.0
    kh = k0h / math.tanh((6/5)**k0h * math.sqrt(k0h))  # initial guess
    for _ in range(max_iter):
        f = k0h - kh * math.tanh(kh)
        df = -math.tanh(kh) - kh / math.cosh(kh) ** 2
        dkh = f / df
        kh_new = kh - dkh
        if abs(dkh / kh) < tol:
            return kh_new
        kh = kh_new
    return kh

# -------------------------------
# Compute Padé approximant coefficients (with rescaled domain)
# -------------------------------
def compute_pade_coeffs(requested_degree, num_points=300):
    """
    Fit a Padé approximant of the form:
       R(u) = u · (p₀ + p₁·v² + … + p_M·v^(2M)) / (1 + q₁·v² + … + q_M·v^(2M))
    to approximate g(u) = kh_numeric(u²) for u ∈ [0, uₘₐₓ],
    where uₘₐₓ = √(2π) and v = u/uₘₐₓ ∈ [0,1].

    If an even degree is requested, it is reduced by one.
    Here, used_degree = 2M+1 and we solve for x = [p₀, …, p_M, q₁, …, q_M].

    Returns:
      requested_degree, used_degree, p (numpy array, length M+1), q (numpy array, length M)
    """
    u_max = np.sqrt(2 * math.pi)
    used_degree = requested_degree if requested_degree % 2 == 1 else requested_degree - 1
    M = (used_degree - 1) // 2  # so that used_degree = 2M + 1

    # Use Chebyshev nodes for v in [0,1]:
    theta = np.linspace(0, np.pi, num_points)
    v_vals = 0.5 * (1 - np.cos(theta))  # Chebyshev nodes mapped to [0,1]
    u_vals = u_max * v_vals             # u = u_max * v, so u ∈ [0, u_max]
    g_vals = np.array([kh_numeric(u ** 2) for u in u_vals])
    
    # Build the design matrix A and right-hand side b.
    A = np.zeros((num_points, 2 * M + 1))
    for i, v in enumerate(v_vals):
        # Numerator part: columns 0 to M (for p₀ ... p_M)
        for j in range(M + 1):
            A[i, j] = u_max * (v ** (2 * j + 1))
        # Denominator part: columns M + 1 to 2 * M (for q₁ ... q_M)
        for k in range(1, M + 1):
            A[i, M + k] = -g_vals[i] * (v ** (2 * k))
    b = g_vals.copy()

    # --- Column scaling to improve conditioning ---
    col_scales = np.max(np.abs(A), axis=0)
    col_scales[col_scales == 0] = 1.0
    A_scaled = A / col_scales
    x_scaled, residuals, rank, s = np.linalg.lstsq(A_scaled, b, rcond=None)
    x = x_scaled / col_scales

    p = x[:M + 1]
    q = x[M + 1:]
    return requested_degree, used_degree, p, q

# -------------------------------
# Evaluate the Padé approximant with rescaled domain
# -------------------------------
def pade_approximation(k0h, p, q):
    """
    Evaluate the Padé approximant for a given k0h.
    Given:
       u = √(k0h) and v = u/√(2π),
    the approximant is:
       R(u) = u · (p₀ + p₁·v² + ... + p_M·v^(2M)) / (1 + q₁·v² + ... + q_M·v^(2M))
    """
    u_max = math.sqrt(2 * math.pi)
    u = math.sqrt(k0h)
    v = u / u_max
    M = len(p) - 1
    num = 0.0
    for j in range(M + 1):
        num += p[j] * (v ** (2 * j))
    num *= u_max * v  # recall u_max*v equals u
    denom = 1.0
    for k in range(1, M + 1):
        denom += q[k - 1] * (v ** (2 * k))
    return num / denom

# -------------------------------
# Main driver
# -------------------------------
def main():
    # Define full domain: k₀h in [0, 2π]
    k0h_vals = np.linspace(0, 2 * math.pi, 1000)
    exact_vals = np.array([kh_numeric(val) for val in k0h_vals])

    # Compute approximants for requested degrees.
    unique_results = {}
    requested_degrees = list(range(5, 31))
    for deg in requested_degrees:
        req_deg, used_deg, p, q = compute_pade_coeffs(deg, num_points=300)
        if used_deg not in unique_results:
            approx_vals = np.array([pade_approximation(val, p, q) for val in k0h_vals])
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_errors = np.where(np.abs(exact_vals) > 1e-12,
                                      np.abs(approx_vals - exact_vals) / np.abs(exact_vals),
                                      0.0)
            avg_err = np.mean(rel_errors)
            max_err = np.max(rel_errors)
            unique_results[used_deg] = (req_deg, used_deg, p, q, avg_err, max_err, approx_vals)

    # Convert dictionary to a list for later processing.
    series_results = list(unique_results.values())

    # Write a detailed report to pade_output.txt
    with open("pade_output.txt", "w", encoding="utf-8") as f:
        f.write("=== Summary of Rescaled Padé Approximations for kh_numeric ===\n")
        f.write("Fitting is performed for k₀h in [0, 2π] with u = √(k₀h) and v = u/√(2π) ∈ [0,1]\n")
        f.write("The Padé approximant is of the form:\n")
        f.write("   R(u) = u · (p₀ + p₁·v² + ... + p_M·v^(2M)) / (1 + q₁·v² + ... + q_M·v^(2M))\n\n")
        for req_deg, used_deg, p, q, avg_err, max_err, _ in series_results:
            f.write(f"Requested Degree: {req_deg}   (Used Degree: {used_deg})\n")
            f.write("Numerator coefficients (for u_max*v^(2j+1) terms):\n")
            for j, coeff in enumerate(p):
                f.write(f"  Coefficient for v^{2*j} : {coeff:.12e}\n")
            f.write("Denominator coefficients (for v^(2k) terms, k>=1):\n")
            for k, coeff in enumerate(q, start=1):
                f.write(f"  Coefficient for v^{2*k} : {coeff:.12e}\n")
            f.write(f"Average relative error: {avg_err*100:.6e} %\n")
            f.write(f"Maximum relative error: {max_err*100:.6e} %\n")
            f.write("-" * 60 + "\n")

        # Detailed table with 50 sample points
        f.write("\n=== Detailed Sample Evaluation (50 sample points) ===\n")
        header = "k0h".ljust(10) + "Exact kh".ljust(15)
        for req_deg, used_deg, _, _, _, _, _ in series_results:
            header += f"Pade(Deg {req_deg})".ljust(20) + f"RelErr({req_deg})".ljust(15)
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")

        sample_indices = np.linspace(0, len(k0h_vals) - 1, 50, dtype=int)
        for idx in sample_indices:
            line = f"{k0h_vals[idx]:<10.4f}" + f"{exact_vals[idx]:<15.6f}"
            for req_deg, used_deg, p, q, _, _, _ in series_results:
                approx_val = pade_approximation(k0h_vals[idx], p, q)
                exact_val = kh_numeric(k0h_vals[idx])
                rel_err = (abs(approx_val - exact_val) / abs(exact_val) * 100 
                           if abs(exact_val) > 1e-12 else 0.0)
                line += f"{approx_val:<20.6f}" + f"{rel_err:<15.6f}"
            f.write(line + "\n")

    # Output pade_routines.txt containing a single function pade2025()
    with open("pade_routines.txt", "w", encoding="utf-8") as f:
        f.write("import math\n\n")
        f.write("def pade2025(k0h, formula):\n")
        f.write('    """\n')
        f.write("    Parameters:\n")
        f.write("      k0h    : Nondimensional deep-water parameter.\n")
        f.write(f"      formula: Integer (1 to {len(series_results)}) specifying which formula to use.\n")
        f.write("    Returns:\n")
        f.write("      Approximated nondimensional wavenumber, kh.\n")
        f.write('    """\n')

        # Write the if/elif chain using the unique series_results.
        for idx, result in enumerate(series_results, start=1):
            req_deg, used_deg, p, q, avg_err, max_err, _ = result
            if idx == 1:
                f.write("    if formula == 1:\n")
            else:
                f.write(f"    elif formula == {idx}:\n")

            M = len(p) - 1
            
            # Create sympy symbols
            k0h_sym = sp.Symbol('k0h')
            u_max = math.sqrt(2 * math.pi)
            v_sym = k0h_sym**0.5 / u_max  # v = sqrt(k0h) / sqrt(2π)

            # Build pure polynomial numerators and denominators
            num_expr = sum(coef * (v_sym**(2*j)) for j, coef in enumerate(p))
            denom_expr = 1 + sum(coef * (v_sym**(2*k)) for k, coef in enumerate(q, start=1))

            # Convert to plain polynomial expressions in terms of k0h
            num_expr = sp.simplify(num_expr * (k0h_sym**0.5))  # Include u as part of the numerator
            denom_expr = sp.simplify(denom_expr)

            # Format for output
            formatted_num_expr = sp.expand(num_expr)
            formatted_denom_expr = sp.expand(denom_expr)

            f.write(f"        return ( ({formatted_num_expr}) ) / ( ({formatted_denom_expr}) )\n")

        f.write("    else:\n")
        f.write("        return -1\n")

    # -------------------------------
    # Plotting the approximants against the exact solution.
    # -------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(k0h_vals, exact_vals, 'k-', linewidth=2, label='kh_numeric (exact)')
    colors = plt.cm.plasma(np.linspace(0, 1, len(series_results)))
    for (req_deg, used_deg, p, q, avg_err, max_err, approx_vals), col in zip(series_results, colors):
        label_str = f"Req Deg {req_deg} (used {used_deg}): avg err {avg_err*100:.2e}%, max err {max_err*100:.2e}%"
        plt.plot(k0h_vals, approx_vals, '-', color=col, label=label_str)
    plt.xlabel("k0h")
    plt.ylabel("kh")
    plt.title("kh_numeric vs. Padé Approximants (Fitted over [0, 2π])")
    plt.legend(fontsize=8, loc="upper left")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('pade_plot.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()