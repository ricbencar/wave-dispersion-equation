/*
 * ============================================================
 * WAVE DISPERSION EQUATION – Linear Gravity Waves (Airy Theory)
 * ============================================================
 *
 * This module provides a comprehensive suite of solutions and approximations for analyzing wave dispersion,
 * essential for wave prediction, oceanographic calculations, and coastal engineering design. It includes a
 * reference "exact" solution using the Newton-Raphson method, classical and contemporary explicit approximations,
 * and high-order Padé approximants for high precision.
 *
 * **Background:**
 *
 * The **linear wave dispersion equation** relates wave frequency (or period) to wavenumber and water depth
 * for gravity waves:
 *
 *     ω² = g · k · tanh(k · h)
 *
 * where:
 *     - ω (omega) is the angular frequency (ω = 2π/T, T = wave period),
 *     - g is the gravitational acceleration,
 *     - k is the wavenumber (k = 2π/L, L = wavelength),
 *     - h is the water depth.
 *
 * The transcendental nature of this equation prevents closed-form solutions for *k*. Therefore, it is
 * nondimensionalized to:
 *
 *     k₀h = kh · tanh(kh)
 *
 * This equation is solved iteratively for the dimensionless wavenumber *kh* (k₀ = ω²/g). An accurate
 * *kh* evaluation is vital for computing wave phase speed, group velocity, and understanding various
 * nearshore processes. Explicit approximations bypass the need for iteration but must be chosen
 * carefully based on accuracy requirements.
 *
 * **Module Contents:**
 *
 *   - Reference "Exact" Solution: kh_numeric() implements the Newton-Raphson iteration method for a
 *     highly precise solution of wave dispersion, acting as a benchmark for other techniques.
 *
 *   - Classical Approximations: Established methods from researchers like Hunt, Eckart, Nielsen, and Gilbert.
 *
 *   - Contemporary Approximations: Recent techniques from researchers such as Guo, Beji, Vatankhah &
 *     Aghashariatmadari, Simarro & Orfila, Yu, Fenton & McKee, Guan & Ju, and Iwagaki.
 *
 *   - High-Order Padé Approximations: Carvalho's 2025 high-order Padé approximants deliver exceptional
 *     precision, addressing increasing complexity in wave calculations as a robust alternative to simpler
 *     methods.
 *
 * ## Compilation
 *
 * To compile the program, use the following command:
 *
 * ```sh
 * g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic
 *     -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -lm
 *     -o wave-disp-equation wave-disp-equation.cpp
 * ```
 *
 * **Explanation of compile options:**
 *   - `-O3`                   : Enables high-level optimizations for maximum performance.
 *   - `-fopenmp`              : Enables OpenMP support for multi-threading.
 *   - `-march=native`         : Optimizes the code for the architecture of the compiling machine.
 *   - `-std=c++17`            : Uses the C++17 standard.
 *   - `-Wall -Wextra -pedantic`: Activates a broad set of compiler warnings to ensure code quality.
 *   - `-Wconversion`          : Warns about implicit type conversions.
 *   - `-Wsign-conversion`     : Warns about implicit sign conversions.
 *   - `-static, -static-libgcc, -static-libstdc++`: Links libraries statically, enhancing portability.
 *   - `-lm`                   : Explicitly link the math library (sometimes needed).
 *
 * ## Usage
 *
 * After compilation, run the program with:
 *
 * ```sh
 * ./wave-disp-equation.exe
 * ```
 *
 * The program will calculate and display error statistics for various wave dispersion approximation methods
 * compared against the reference Newton-Raphson solution.
 *
 * **References:**
 *
 *  1. Wikipedia. "Airy wave theory". [Wikipedia](https://en.wikipedia.org/wiki/Airy_wave_theory).
 *  2. Wikipedia. "Dispersion (water waves)". [Wikipedia](https://en.wikipedia.org/wiki/Dispersion_(water_waves)).
 *  3. Yu, J. (2014). "A Note on Approximations of the Dispersion Relationship of Water Waves",
 *     *Journal of Engineering Mechanics (ASCE)*, 140(1), 233–237.
 *  4. Simarro, G. & Orfila, A. (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves:
 *     Another discussion", *Coastal Engineering*, 80, 15–16.
 *  5. Vatankhah, A.R. & Aghashariatmadari, Z. (2013). "Improved explicit approximation of linear dispersion relationship
 *     for gravity waves: a discussion", *Coastal Engineering*, 78, 21–22.
 *  6. Beji, S. (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves",
 *     *Coastal Engineering*, 73, 11–12.
 *  7. **You, Z.J.** (2008). "A close approximation of wave dispersion relation for direct calculations",
 *     *Applied Ocean Research*, 30(2), 141–143.
 *  8. Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
 *     *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
 *  9. Yamaguchi, M. and H. Nonaka: Comparative study of explicit solutions to wave dispersion equation,
 *     *Annu. Jour. Eng.*, Ehime Univ., Vol. 6, 2007 in CD-ROM.
 * 10. Carvalho, R. (2006). Unpublished work based on gene expression programming for wave dispersion equations.
 * 11. You, Z.J. "Discussion of 'Simple and explicit solution to the wave dispersion equation'",
 *     [Coastal Engineering 45 (2002) 71-74], Coastal Eng., Vol. 48, pp.133-135, 2003.
 * 12. Guo, J. (2002). "Simple and explicit solution of the wave dispersion equation",
 *     *Coastal Engineering*, 45, 71–74.
 * 13. Fenton, J.D. & McKee, W.D. (1990). "On calculating the lengths of water waves",
 *     *Coastal Engineering*, 14, 499–513.
 * 14. Fenton, J.D. "The numerical solution of steady water wave problems", *Computers & Geosciences*,
 *     Vol. 4, No. 3, pp.357-368, 1988.
 * 15. Fenton, J.D. (1972). "A ninth-order solution for the solitary wave",
 *     *Journal of Fluid Mechanics*, 53, 257–271.
 * 16. Wu, C. S. and E. B. Thornton. "Wave numbers of linear progressive waves",
 *     *Journal of Waterway, Port, Coastal and Ocean Engineering*, ASCE, Vol. 112, No. 4, pp.536-540, 1986.
 * 17. Nielsen, P. "Explicit solutions to practical wave problems", Proc. 19th ICCE, Vol. 1, pp.968-982, 1984.
 * 18. Nielsen, P. "Explicit formulae for practical wave calculations", Coastal Eng., No. 6, pp.389-398, 1982.
 * 19. Hunt, J.N. (1979). "A simple approximation for the dispersion relation of water waves",
 *     *Journal of Waterway, Port, Coastal and Ocean Engineering*, 105(4), 457–459.
 * 20. Hunt, J. N. "Direct solution of wave dispersion equation", *J. Waterway, Port, Coastal and Ocean Div.*,
 *     Proc. ASCE, Vol. 105, No. WW4, pp.457-459, 1979.
 * 21. Eckart, C. (1951). "The propagation of gravity waves from deep to shallow water",
 *     U.S. Department of Commerce, National Bureau of Standards Circular 521.
 *
 * These references have been arranged by both historical and topical relevance. They span classical
 * methods (e.g., Eckart, 1951; Hunt, 1979), modern explicit approximations (e.g., Guo, 2002; Beji, 2013;
 * Vatankhah & Aghashariatmadari, 2013; Simarro & Orfila, 2013; Yu, 2014), as well as pivotal contributions
 * from Fenton and colleagues (Fenton & McKee, 1990; Fenton, 1972) and innovative computational approaches
 * (Carvalho, 2006 & 2025) to improve the dispersion relation accuracy.
 */

// Standard library includes
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <limits>
#include <functional>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cstdlib>

// Use namespaces to simplify code
using namespace std;

// Constants
const double PI = 3.14159265358979323846;

// =============================================================================
// EXACT SOLUTION (NEWTON–RAPHSON) - Reference Implementation
// =============================================================================

/**
 * Compute the 'exact' nondimensional wavenumber *kh* by numerically solving the dispersion
 * relation using the Newton–Raphson method.
 *
 * **Equation Solved:**
 * The nondimensional dispersion relation (derived from Airy wave theory) is:
 *     f(kh) = k0h - kh * tanh(kh) = 0
 * where k0h = k₀·h (with k₀ = ω²/g) and kh = k·h.
 *
 * **Method:**
 *   - Newton–Raphson iteration:
 *         kh_new = kh - f(kh) / f'(kh)
 *   - f'(kh) = -tanh(kh) - kh * sech²(kh)
 *   - Initial guess: kh₀ ≈ k0h / tanh((6/5)^k0h * sqrt(k0h)) (Carvalho, 2006 style).  This initialization provides
 *     a reasonable starting point and promotes quicker convergence.
 *
 * **Convergence Criteria:**
 *   Iteratively adjust *kh* until the relative change |Δkh/kh| is below the tolerance `tol`.
 *
 * **References:**
 *   - Fenton & McKee (1990); Yamaguchi & Nonaka (2007); Press et al. (1992).
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter (k₀·h).  Must be non-negative.
 *   tol (float): Relative convergence tolerance (default: 1e-15). Adjust for desired precision.
 *   max_iter (int): Maximum number of iterations (default: 100). Increase if convergence fails.
 *
 * **Returns:**
 *   float: Computed nondimensional wavenumber *kh*. Returns 0.0 if `k0h` is 0.
 */
double kh_numeric(double k0h, double tol = 1e-15, int max_iter = 100)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double kh = k0h / tanh(pow(6.0 / 5.0, k0h) * sqrt(k0h));
    for (int i = 0; i < max_iter; ++i)
    {
        double f = k0h - kh * tanh(kh);
        double df = -tanh(kh) - kh / pow(cosh(kh), 2);
        double dkh = f / df;
        double kh_new = kh - dkh;
        if (abs(dkh / kh) < tol)
        {
            return kh_new;
        }
        kh = kh_new;
    }
    return kh;
}

// =============================================================================
// PADE APPROXIMANT - a ratio of two power series
// =============================================================================

/**
 * Approximations using Padé approximants for the nondimensional wavenumber.
 *
 * **Description:**
 * Padé approximants are rational functions that approximate a given function by matching its Taylor
 * series expansion up to a specified order. They use ratios of polynomials instead of solely
 * polynomial expansions as in Taylor series, offering better accuracy, particularly around
 * singularities and for functions exhibiting poles. This suite of approximations corresponds to
 * approximations derived in Carvalho (2025) using gene expression programming, to fit a ratio
 * of power series for efficient wave dispersion estimates.
 *
 * **Advantages of Padé Approximants:**
 * - Improved Accuracy: Approximations are often superior to Taylor series, mainly in regions where behavior
 *   cannot be fully explained through polynomial formulations.
 * - Convergence Characteristics: Demonstrate accelerated convergence with cases that may be divergent for Taylor series.
 * - Singularity Representations: Allows effective characterization close to points involving singular values or poles.
 * - Analytic Continuations: It makes provision to broaden approximation to areas beyond standard convergence boundaries with conventional Taylor functions.
 *
 * **Error Characteristics:**
 *   Potential errors can arise from numerical instability in the polynomial evaluation for `k0h > 2π.
 *   Performance assessments must be conducted for choosing ideal formulation to the given input.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter (k₀·h). Should be >=0 and <= 2π.
 *   formula (int): An integer (1 to 13) indicating which formula to compute and use.
 *
 * **Returns:**
 *   float: An approximation to nondimensional wavenumber, kh.
 *   returns -1.0 when 'formula' parameters are inappropriate / non compliant of input regulations.
 *
 * **References:**
 *     - R. Carvalho (2025). Work published on GitHub, actually the code you're reading right now.
 */
double pade2025(double k0h, int formula)
{
    if (k0h < 0 || k0h > 2 * PI)
    {
        return -1.0;
    }
    if (formula < 1 || formula > 13)
    {
        return -1.0;
    }

    // Padé approximant formulas with coefficients derived from gene expression programming
    if (formula == 1)
    {
        return ((1.00649052194019 * pow(k0h, 0.5) + 0.423646282789217 * pow(k0h, 1.5) + 0.175406661440005 * pow(k0h, 2.5))) /
               ((0.306955955676234 * pow(k0h, 1.0) + 0.0328975279727171 * pow(k0h, 2.0) + 1));
    }
    else if (formula == 2)
    {
        return ((0.998980252114366 * pow(k0h, 0.5) + 0.0240176797055886 * pow(k0h, 1.5) + 0.102524886754552 * pow(k0h, 2.5) + 0.0317327085938995 * pow(k0h, 3.5))) /
               ((-0.150350405960952 * pow(k0h, 1.0) + 0.112157962910113 * pow(k0h, 2.0) + 0.00294483072586115 * pow(k0h, 3.0) + 1));
    }
    else if (formula == 3)
    {
        return ((1.00006668638419 * pow(k0h, 0.5) + 0.322645945302282 * pow(k0h, 1.5) + 0.0860384450810725 * pow(k0h, 2.5) + 0.051143347041175 * pow(k0h, 3.5) + 0.0153420957423937 * pow(k0h, 4.5))) /
               ((0.157166943736625 * pow(k0h, 1.0) + 0.0245168267924732 * pow(k0h, 2.0) + 0.0462567432956417 * pow(k0h, 3.0) + 0.00175392506101448 * pow(k0h, 4.0) + 1));
    }
    else if (formula == 4)
    {
        return ((0.999996682596798 * pow(k0h, 0.5) - 0.0889915717930786 * pow(k0h, 1.5) + 0.147076211695128 * pow(k0h, 2.5) + 0.0123471280480147 * pow(k0h, 3.5) + 0.00866458140843225 * pow(k0h, 4.5) + 0.00204463718201973 * pow(k0h, 5.5))) /
               ((-0.255723982020183 * pow(k0h, 1.0) + 0.159493904911975 * pow(k0h, 2.0) - 0.0106101311382749 * pow(k0h, 3.0) + 0.00784491418150148 * pow(k0h, 4.0) + 0.000184273251439305 * pow(k0h, 5.0) + 1));
    }
    else if (formula == 5)
    {
        return ((0.999998218345888 * pow(k0h, 0.5) - 0.424362176674708 * pow(k0h, 1.5) + 0.171875463304611 * pow(k0h, 2.5) - 0.0357487982640122 * pow(k0h, 3.5) + 0.00410625374333464 * pow(k0h, 4.5) - 0.000978753693904127 * pow(k0h, 5.5) - 0.000636955605769902 * pow(k0h, 6.5))) /
               ((-0.591069429462395 * pow(k0h, 1.0) + 0.240083348894323 * pow(k0h, 2.0) - 0.0617593442909405 * pow(k0h, 3.0) + 0.0104920694265126 * pow(k0h, 4.0) - 0.00231970889331938 * pow(k0h, 5.0) - 5.65924775627923e-5 * pow(k0h, 6.0) + 1));
    }
    else if (formula == 6)
    {
        return ((1.0000000012405 * pow(k0h, 0.5) - 0.350251200743747 * pow(k0h, 1.5) + 0.229153326540668 * pow(k0h, 2.5) - 0.0205204312544928 * pow(k0h, 3.5) + 0.0133231478358294 * pow(k0h, 4.5) + 0.0010401274983046 * pow(k0h, 5.5) + 0.00048671850792775 * pow(k0h, 6.5) + 8.40088474488992e-5 * pow(k0h, 7.5))) /
               ((-0.516917882097882 * pow(k0h, 1.0) + 0.284751410622371 * pow(k0h, 2.0) - 0.0555622365621819 * pow(k0h, 3.0) + 0.0161071584333013 * pow(k0h, 4.0) - 0.000808341017586247 * pow(k0h, 5.0) + 0.000381511960690599 * pow(k0h, 6.0) + 6.40735447518177e-6 * pow(k0h, 7.0) + 1));
    }
    else if (formula == 7)
    {
        return ((0.999999995257458 * pow(k0h, 0.5) - 0.543811114837314 * pow(k0h, 1.5) + 0.297774393256421 * pow(k0h, 2.5) - 0.0648661921727468 * pow(k0h, 3.5) + 0.0174768559302056 * pow(k0h, 4.5) - 0.00151793039097231 * pow(k0h, 5.5) + 0.000295750461715408 * pow(k0h, 6.5) - 4.1567697098083e-6 * pow(k0h, 7.5) - 1.62498860684328e-5 * pow(k0h, 8.5))) /
               ((-0.710477969437912 * pow(k0h, 1.0) + 0.385633880896532 * pow(k0h, 2.0) - 0.110812474410256 * pow(k0h, 3.0) + 0.0270500418196004 * pow(k0h, 4.0) - 0.00394602590941703 * pow(k0h, 5.0) + 0.000553921270262525 * pow(k0h, 6.0) - 6.65533846723705e-5 * pow(k0h, 7.0) - 1.24656671282763e-6 * pow(k0h, 8.0) + 1));
    }
    else if (formula == 8)
    {
        return ((1.00000000020126 * pow(k0h, 0.5) - 0.388439115555858 * pow(k0h, 1.5) + 0.310223332529737 * pow(k0h, 2.5) - 0.0496321949056331 * pow(k0h, 3.5) + 0.0293825301580729 * pow(k0h, 4.5) - 0.000149900084396432 * pow(k0h, 5.5) + 0.00139739652490532 * pow(k0h, 6.5) + 0.000171592322622253 * pow(k0h, 7.5) + 4.94349555930422e-5 * pow(k0h, 8.5) + 5.36329658499187e-6 * pow(k0h, 9.5))) /
               ((-0.555105771884377 * pow(k0h, 1.0) + 0.372185262898829 * pow(k0h, 2.0) - 0.0980736559689507 * pow(k0h, 3.0) + 0.036689986163131 * pow(k0h, 4.0) - 0.00440840853967443 * pow(k0h, 5.0) + 0.00139722437126806 * pow(k0h, 6.0) - 3.64106868131082e-6 * pow(k0h, 7.0) + 2.72315576473091e-5 * pow(k0h, 8.0) + 3.80576928768955e-7 * pow(k0h, 9.0) + 1));
    }
    else if (formula == 9)
    {
        return ((1.00000000054683 * pow(k0h, 0.5) - 0.302517970258141 * pow(k0h, 1.5) + 0.216194173804703 * pow(k0h, 2.5) - 0.00911867675112112 * pow(k0h, 3.5) + 0.0131923444114312 * pow(k0h, 4.5) + 0.00197230469011907 * pow(k0h, 5.5) + 0.000583685774943952 * pow(k0h, 6.5) + 0.000117061426950459 * pow(k0h, 7.5) + 2.16577683514269e-5 * pow(k0h, 8.5) - 4.06436904033371e-6 * pow(k0h, 9.5) - 1.35655741732812e-7 * pow(k0h, 10.5))) /
               ((-0.469184609052338 * pow(k0h, 1.0) + 0.263835664763448 * pow(k0h, 2.0) - 0.0421256860069538 * pow(k0h, 3.0) + 0.0141905686336732 * pow(k0h, 4.0) + 0.000171234247175628 * pow(k0h, 5.0) + 0.000284037180073895 * pow(k0h, 6.0) + 7.23728845968762e-5 * pow(k0h, 7.0) - 3.33925175755566e-6 * pow(k0h, 8.0) - 1.1196042731312e-6 * pow(k0h, 9.0) - 6.10361897619335e-9 * pow(k0h, 10.0) + 1));
    }
    else if (formula == 10)
    {
        return ((1.00000000069543 * pow(k0h, 0.5) - 0.334954381847524 * pow(k0h, 1.5) + 0.229997336203832 * pow(k0h, 2.5) - 0.0170932674073579 * pow(k0h, 3.5) + 0.0142286437913579 * pow(k0h, 4.5) + 0.00157766182461221 * pow(k0h, 5.5) + 0.000555414435408452 * pow(k0h, 6.5) + 0.000109241385068583 * pow(k0h, 7.5) + 2.04181414663277e-5 * pow(k0h, 8.5) - 4.47936721320148e-6 * pow(k0h, 9.5) + 7.20245847805242e-8 * pow(k0h, 10.5) + 2.27591359161482e-9 * pow(k0h, 11.5))) /
               ((-0.501621014271701 * pow(k0h, 1.0) + 0.283044819037831 * pow(k0h, 2.0) - 0.0523102821262626 * pow(k0h, 3.0) + 0.0164455383396673 * pow(k0h, 4.0) - 0.000365622576684072 * pow(k0h, 5.0) + 0.000304882740985283 * pow(k0h, 6.0) + 7.15596032735346e-5 * pow(k0h, 7.0) - 5.64260090297976e-6 * pow(k0h, 8.0) - 7.22847827757772e-7 * pow(k0h, 9.0) + 3.16954475149092e-8 * pow(k0h, 10.0) - 4.18523172095159e-11 * pow(k0h, 11.0) + 1));
    }
    else if (formula == 11)
    {
        return ((1.00000000021134 * pow(k0h, 0.5) - 0.439538511010958 * pow(k0h, 1.5) + 0.262071966075091 * pow(k0h, 2.5) - 0.0352260757847662 * pow(k0h, 3.5) + 0.0136888861362354 * pow(k0h, 4.5) + 0.00119983590894612 * pow(k0h, 5.5) + 0.000306132196963262 * pow(k0h, 6.5) + 9.3344593984067e-5 * pow(k0h, 7.5) + 1.13485236045952e-5 * pow(k0h, 8.5) - 2.19033671564094e-6 * pow(k0h, 9.5) - 1.52303393432862e-7 * pow(k0h, 10.5) + 3.29680588537e-8 * pow(k0h, 11.5) - 3.59648926971857e-9 * pow(k0h, 12.5))) /
               ((-0.606205166733527 * pow(k0h, 1.0) + 0.332550449569891 * pow(k0h, 2.0) - 0.0755002415059252 * pow(k0h, 3.0) + 0.0186169006337531 * pow(k0h, 4.0) - 0.000624409536143469 * pow(k0h, 5.0) + 0.000110698204767146 * pow(k0h, 6.0) + 7.31839474973857e-5 * pow(k0h, 7.0) - 4.17221912299122e-6 * pow(k0h, 8.0) - 1.00985652913037e-6 * pow(k0h, 9.0) + 1.15600808535017e-7 * pow(k0h, 10.0) - 7.31093937337407e-9 * pow(k0h, 11.0) - 3.576922099426e-10 * pow(k0h, 12.0) + 1));
    }
    else if (formula == 12)
    {
        return ((1.00000000034658 * pow(k0h, 0.5) - 0.35765423836608 * pow(k0h, 1.5) + 0.220474157851537 * pow(k0h, 2.5) - 0.0143101981637556 * pow(k0h, 3.5) + 0.0106882664277003 * pow(k0h, 4.5) + 0.00177809831685665 * pow(k0h, 5.5) + 0.000425551341246499 * pow(k0h, 6.5) + 8.39454668343144e-5 * pow(k0h, 7.5) + 1.68918586510378e-5 * pow(k0h, 8.5) - 3.47736023040874e-6 * pow(k0h, 9.5) + 1.1564308170262e-7 * pow(k0h, 10.5) - 5.38824115168048e-9 * pow(k0h, 11.5) + 2.55400838802905e-9 * pow(k0h, 12.5) - 5.77246087402385e-10 * pow(k0h, 13.5))) /
               ((-0.524320887192481 * pow(k0h, 1.0) + 0.277305164506462 * pow(k0h, 2.0) - 0.0478782214771677 * pow(k0h, 3.0) + 0.0124223110237922 * pow(k0h, 4.0) + 0.00037129622072912 * pow(k0h, 5.0) + 0.000158246646088704 * pow(k0h, 6.0) + 5.24876551697954e-5 * pow(k0h, 7.0) - 1.38557325635002e-6 * pow(k0h, 8.0) - 7.96232034767084e-7 * pow(k0h, 9.0) + 8.35464703802028e-9 * pow(k0h, 10.0) + 1.1834526426828e-8 * pow(k0h, 11.0) - 1.48873774124982e-9 * pow(k0h, 12.0) - 5.59846037648501e-11 * pow(k0h, 13.0) + 1));
    }
    else if (formula == 13)
    {
        return ((1.00000000043044 * pow(k0h, 0.5) - 0.341214787680155 * pow(k0h, 1.5) + 0.216029315236116 * pow(k0h, 2.5) - 0.0116443516054976 * pow(k0h, 3.5) + 0.0108812744435703 * pow(k0h, 4.5) + 0.00184089184824533 * pow(k0h, 5.5) + 0.000466873470691525 * pow(k0h, 6.5) + 8.83178779786945e-5 * pow(k0h, 7.5) + 1.82969284264395e-5 * pow(k0h, 8.5) - 3.70984687913519e-6 * pow(k0h, 9.5) + 1.15252770559743e-7 * pow(k0h, 10.5) + 4.05866020101675e-10 * pow(k0h, 11.5) + 1.00682972256747e-9 * pow(k0h, 12.5) - 1.03386709606535e-10 * pow(k0h, 13.5) + 2.54165346162395e-11 * pow(k0h, 14.5))) /
               ((-0.507881432557414 * pow(k0h, 1.0) + 0.27012036160695 * pow(k0h, 2.0) - 0.0445169260696995 * pow(k0h, 3.0) + 0.0122182321499095 * pow(k0h, 4.0) + 0.00039913233839911 * pow(k0h, 5.0) + 0.000190210573213751 * pow(k0h, 6.0) + 5.3696821111172e-5 * pow(k0h, 7.0) - 1.69377956737975e-6 * pow(k0h, 8.0) - 7.45889615784955e-7 * pow(k0h, 9.0) + 1.5632250659559e-9 * pow(k0h, 10.0) + 1.19811786331838e-8 * pow(k0h, 11.0) - 1.3440816876962e-9 * pow(k0h, 12.0) + 1.13151579925971e-10 * pow(k0h, 13.0) + 1.46210486272321e-12 * pow(k0h, 14.0) + 1));
    }
    else
    {
        return -1.0; // Should never get here with the checks at the top
    }
}

// =============================================================================
// CARVALHO (2025) GEP-based approximations
// =============================================================================

/**
 * Approximations using Carvalho's (2025) Gene Expression Programming (GEP) solutions
 * for estimating the nondimensional wavenumber *kh*.
 *
 * **Description:**
 * These approximations, derived via Gene Expression Programming (GEP) (Carvalho, 2006 & 2025),
 * represent one of the early attempts to utilize GEP to fit a closed-form,
 * explicit expression for the dispersion relation over the whole range of nondimensional water depths.
 *
 * **Method:**
 * This function employs a suite of pre-computed GEP formulas. Each formula provides a different
 * algebraic expression estimating *kh* based on *k0h*. The *formula* parameter selects the
 * specific GEP expression to use. Many expressions incorporate combinations of hyperbolic functions,
 * power functions, and exponential terms designed to mimic the behavior of the exact dispersion relation.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter (k₀·h). Must be non-negative.
 *   formula (int): An integer ranging from 1 to 20, inclusive, selecting which GEP formula to use.
 *
 * References:
 *   - Carvalho, R. (2006). Unpublished work; see also Yamaguchi & Nonaka (2007) for discussion.
 *   - Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
 *     *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
 *   - Ferreira, C. (2006). *Gene Expression Programming: Mathematical Modeling by an Artificial Intelligence*.
 *     2nd Edition. Springer-Verlag, Germany. This book provides a comprehensive introduction to Gene Expression
 *     Programming (GEP), a novel genetic algorithm that uses the principles of evolution to solve complex problems
 *     in mathematics and modeling. Available at: [Springer Link](https://link.springer.com/book/10.1007/3-540-32849-1).
 *   - Ferreira, C. (2004). "Gene Expression Programming and the Evolution of Computer Programs".
 *     In Leandro N. de Castro and Fernando J. Von Zuben (Eds.), *Recent Developments in Biologically Inspired Computing*,
 *     pp. 82-103. Idea Group Publishing.
 *   - Ferreira, C. (2002). "Gene Expression Programming in Problem Solving". In R. Roy, M. Káppen,
 *     S. Ovaska, T. Furuhashi, and F. Hoffmann (Eds.), *Soft Computing and Industry: Recent Applications*,
 *     pp. 635-654. Springer-Verlag.
 *   - Ferreira, C. (2002). *Gene Expression Programming: Mathematical Modeling by an Artificial Intelligence*.
 *     Angra do Heroismo, Portugal. Online version.
 *   - Ferreira, C. (2001). "Gene Expression Programming: A New Adaptive Algorithm for Solving Problems",
 *     *Complex Systems*, 13(2), 87-129.
 *
 * Returns:
 *   float: Approximated nondimensional wavenumber, kh.
 */
double carvalho2025(double k0h, int formula)
{
    // Input validation
    if (k0h < 0 || formula < 1 || formula > 20)
    {
        return -1.0;
    }

    // Special edge case
    if (k0h == 0.0)
    {
        return 0.0;
    }

    // Implement exact formulas from Python version
    switch (formula)
    {
    case 1:
    {
        double kh_Carv = k0h / tanh(k0h / (sqrt(tanh(sqrt(sinh(k0h)))) * pow(tanh(k0h), 0.25)));
        return (kh_Carv * kh_Carv + k0h * pow(cosh(kh_Carv), 2)) / (kh_Carv + sinh(kh_Carv) * cosh(kh_Carv));
    }

    case 2:
    {
        if (k0h <= 1.2)
        {
            return sqrt(1 / k0h - exp(pow(k0h, 1.962983) - 6.242035)) / (1 / k0h - 0.168659434);
        }
        else if (1.2 < k0h && k0h <= 2.35)
        {
            return (k0h + pow((k0h / 70.13327717), pow(k0h, 3))) / exp(log(pow(4.89859, k0h)) / (1.134674 - pow(10, k0h)));
        }
        else
        { // k0h > 2.35
            return k0h * exp(1.596671172 * k0h / pow(10, k0h));
        }
    }

    case 3:
        return k0h / tanh(
                         k0h / (sqrt(tanh(sqrt(sinh(k0h)))) * pow(tanh(k0h), 0.25)));

    case 4:
        return k0h / (tanh(k0h / tanh(k0h / tanh(k0h / sinh(tanh(sqrt(k0h)))))));

    case 5:
        return k0h / (tanh(pow(1.199315, pow(k0h, 1.047086)) * pow(k0h, 0.499947)));

    case 6:
        return k0h / tanh(pow(1.1999, pow(k0h, 1.045)) * sqrt(k0h));

    case 7:
        return k0h / tanh(k0h / tanh(k0h / sinh(tanh(sqrt(k0h)))));

    case 8:
        return k0h / tanh(
                         sinh(sqrt(k0h >= 3.04425 ? 3.04425 : k0h)) * cosh(k0h / 5.194671));

    case 9:
        return k0h / (sqrt(tanh(sqrt(sinh(k0h)))) * pow(tanh(k0h), 0.25));

    case 10:
        return k0h / tanh(pow((6.0 / 5.0), k0h) * sqrt(k0h));

    case 11:
        return k0h / tanh(sqrt(pow(1.438995, k0h) * k0h));

    case 12:
        return k0h / tanh(k0h / tanh(sinh(sqrt(k0h))));

    case 13:
        return k0h + sqrt(k0h) / (pow(4.35144, k0h) + 0.718409 / pow((1.0 / k0h), 0.437408));

    case 14:
        return k0h / (pow(tanh(sqrt(k0h)), (1.0 / cosh(k0h))));

    case 15:
        return k0h / (sqrt(tanh(k0h)) * tanh(k0h + 1.0 / sqrt(k0h)));

    case 16:
        return k0h / pow(tanh(k0h), ((k0h + 4.0) / 8.0));

    case 17:
        return k0h / pow((pow(tanh(k0h), (k0h / tanh(k0h)))), 0.5);

    case 18:
        return k0h / tanh(sinh(sqrt(k0h)));

    case 19:
        return sqrt(k0h) + pow(k0h, 2) / (k0h + 4.0);

    case 20:
        return k0h / pow((pow(sqrt(tanh(k0h)), (tanh(k0h) + 4.0))), 0.25);

    default:
        return -1.0;
    }
}

// =============================================================================
// YAMAGUCHI & NONANKA (2007) family of explicit solutions
// =============================================================================

/**
 * Yamaguchi & Nonaka (2007) family of explicit solutions (YN1–YN10).
 *
 * Yamaguchi and Nonaka (2007) introduced a set of ten explicit formulas (YN1–YN10) to approximate
 * the wave dispersion relation in linear wave theory. These formulas provide accurate and
 * computationally efficient alternatives to the implicit dispersion equation, which usually requires
 * iterative numerical methods for exact solutions.
 *
 * Overview:
 * The Yamaguchi & Nonaka (2007) approximations eliminate the need for iterative procedures,
 * enhancing computational efficiency in fields like coastal engineering and oceanography.
 * Each formula employs mathematical functions such as hyperbolic cotangent (coth), hyperbolic tangent (tanh),
 * and exponential functions to approximate the dispersion relation across various depth conditions.
 *
 * **Parameters:**
 * k0h (float): Nondimensional deep-water parameter, defined as the product of the wavenumber (k0) and water depth (h).
 * formula (int): An integer (1 to 10) specifying which Yamaguchi & Nonaka (2007) formula (YN1–YN10) to use.
 *
 * References:
 *   - Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
 *     *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
 *   - Yamaguchi, M. and H. Nonaka: Comparative study of explicit solutions to wave dispersion equation,
 *     *Annu. Jour. Eng.*, Ehime Univ., Vol. 6, 2007 in CD-ROM.
 *
 * **Returns:**
 * Approximated nondimensional wavenumber, kh.
 */
double YamaguchiNonaka(double k0h, int formula)
{
    if (k0h == 0)
    {
        return 0.0;
    }

    if (formula < 1 || formula > 10)
    {
        return -1.0;
    }

    auto coth = [](double x) -> double
    {
        if (x == 0)
        {
            return numeric_limits<double>::infinity();
        }
        else
        {
            return cosh(x) / sinh(x);
        }
    };

    if (formula == 1)
    {
        return k0h * pow(coth(pow(k0h, 1.485 / 2)), (1 / 1.485));
    }
    else if (formula == 2)
    {
        return k0h / tanh(k0h * pow(coth(pow(k0h, 1.378 / 2)), (1 / 1.378)));
    }
    else if (formula == 3)
    {
        return k0h / tanh(sqrt(k0h) * (1.0 + sqrt(k0h) / (2.0 * PI)));
    }
    else if (formula == 4)
    {
        return k0h * pow(1.0 + 1.0 / (k0h * k0h), 0.25);
    }
    else if (formula == 5)
    {
        return k0h * pow(coth(pow(k0h, 1.434 / 2)), (1 / 1.434));
    }
    else if (formula == 6)
    {
        return k0h / tanh(sqrt(sinh(k0h)));
    }
    else if (formula == 7)
    {
        return k0h / pow((1 - exp(-pow(k0h, 2.445 / 2))), (1 / 2.445));
    }
    else if (formula == 8)
    {
        return k0h / tanh(k0h * pow(coth(pow(k0h, 1.310 / 2)), (1 / 1.310)));
    }
    else if (formula == 9)
    {
        return k0h / tanh(pow(1.1965, k0h) * sqrt(k0h));
    }
    else if (formula == 10)
    {
        return k0h / tanh(k0h / (sqrt(tanh(sqrt(sinh(k0h)))) * pow(tanh(k0h), 0.25)));
    }
    else
    {
        return -1.0; // Should never get here with the checks at the top
    }
}

/**
 * Beji (2013) improved explicit approximation.
 *
 * **Formulation:**
 *   kh ≈ [k0h/√(tanh(k0h))] · [1 + k0h^1.09 · exp(–1.55 – 1.30·k0h – 0.216·k0h²)].
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double beji2013(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }

    // Compute the exponential correction term
    double exp_term = exp(-(1.55 + 1.30 * k0h + 0.216 * k0h * k0h));

    return (k0h * (1 + pow(k0h, 1.09) * exp_term)) / sqrt(tanh(k0h));
}

/**
 * Simarro & Orfila (2013) two-step Newton-corrected approximation.
 *
 * **Concept:**
 *   Uses Beji's approximation as the initial guess and then applies one Newton–Raphson correction:
 *       kh* = [ (kh_B)² + k0h·cosh²(kh_B) ] / [ kh_B + sinh(kh_B)·cosh(kh_B) ],
 *   where kh_B is the Beji estimate.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Corrected approximated nondimensional wavenumber, kh.
 */
double Simarro_2013(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double kh_Beji = beji2013(k0h);
    return (kh_Beji * kh_Beji + k0h * pow(cosh(kh_Beji), 2)) /
           (kh_Beji + sinh(kh_Beji) * cosh(kh_Beji));
}

/**
 * Vatankhah & Aghashariatmadari (2013) – Single-step explicit formula #1.
 *
 * **Formulation:**
 *   kh ≈ [k0h + k0h² · exp(–1.835 – 1.225·k0h^1.35)] / √(tanh(k0h)).
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double vatankhah2013_2(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    return (k0h + k0h * k0h * exp(-1.835 - 1.225 * pow(k0h, 1.35))) / sqrt(tanh(k0h));
}

/**
 * Vatankhah & Aghashariatmadari (2013) – Single-step explicit formula #2.
 *
 * **Concept:**
 *   Splits the approximation into two parts:
 *     partA = [k0h + k0h² · exp(–(3.2 + k0h^1.65))] / √(tanh(k0h)),
 *     partB = k0h · [1 – exp(–k0h^0.132)]^(5.0532 + 2.1584·k0h^1.505).
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double vatankhah2013_1(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double partA = (k0h + k0h * k0h * exp(-(3.2 + pow(k0h, 1.65)))) / sqrt(tanh(k0h));
    double partB = k0h * pow(1 - exp(-pow(k0h, 0.132)), 5.0532 + 2.1584 * pow(k0h, 1.505));
    return partA + partB;
}

/**
 * Hunt (1979) Padé-type rational approximation for the dispersion relation.
 *
 * **Purpose:**
 *   Provides an explicit expression for kh by approximating tanh(kh) via a rational (Padé) function.
 *
 * **Formulation:**
 *   kh ≈ √[ k0h² + k0h / (1 + c₁·k0h + c₂·k0h² + … + c₉·k0h⁹) ],
 *   where the coefficients c₁,…,c₉ are from Hunt (1979).
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double hunt1979_9(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }

    // Hunt's coefficients for the 9th-order approximation
    const double coeffs[] = {
        0.6666666667, 0.3555, 0.16084, 0.0632, 0.02174,
        0.00654, 0.00171, 0.00039, 0.00011};

    // Compute the sum of terms
    double sum = 0.0;
    for (int i = 0; i < 9; ++i)
    {
        sum += coeffs[i] * pow(k0h, i + 1);
    }

    return sqrt(k0h * k0h + k0h / (1 + sum));
}

/**
 * Hunt (1979) – 5th-order approximate solution (Hunt1) from Yamaguchi & Nonaka (2007).
 *
 * **Formulation:**
 *   (kₐ·h)² = α · [ α + 1 / (1 + 0.6522·α + 0.4622·α² + 0.0864·α⁴ + 0.0675·α⁵) ],
 *   where α = k0h.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter (α).
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber (kₐ·h).
 */
double hunt1979_5(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }

    double alpha = k0h;
    // Compute denominator polynomial
    double denom = 1.0 + 0.6522 * alpha + 0.4622 * (alpha * alpha) + 0.0864 * pow(alpha, 4) + 0.0675 * pow(alpha, 5);

    return sqrt(alpha * (alpha + 1.0 / denom));
}

/**
 * Fenton & McKee (1990) iterative-type approximation for kh.
 *
 * **Formulation:**
 *   Compute βₐ = k0h · [coth(k0h)]^(1/2), then
 *       kh = [ k0h + βₐ² · sech²(βₐ) ] / [ tanh(βₐ) + βₐ · sech²(βₐ) ].
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximate nondimensional wavenumber, kh.
 */
double fenton_mckee1990_1(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double alpha = k0h;
    double coth_alpha = (alpha == 0) ? numeric_limits<double>::infinity() : cosh(alpha) / sinh(alpha);
    double beta_a = alpha * sqrt(coth_alpha);
    double sech_beta_a = 1.0 / cosh(beta_a);
    double tanh_beta_a = tanh(beta_a);
    double beta_a_sq = beta_a * beta_a;
    double sech_sq = sech_beta_a * sech_beta_a;
    double numerator = alpha + beta_a_sq * sech_sq;
    double denominator = tanh_beta_a + beta_a * sech_sq;
    return numerator / denominator;
}

/**
 * Fenton & McKee (1990) all-depth empirical approximation for kh.
 *
 * **Formulation:**
 *   kh ≈ k0h / [ tanh(k0h^(3/4)) ]^(2/3).
 *
 * **Parameters:**
 *   k0h (float): Deep-water parameter.
 *
 * **Returns:**
 *   float: Approximate nondimensional wavenumber, kh.
 */
double fenton_mckee1990_2(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    return k0h / pow(tanh(pow(k0h, 0.75)), 2.0 / 3.0);
}

/**
 * Wu & Thornton (1986) explicit approximation for the dispersion relation.
 *
 * **Overview:**
 *   Provides a piecewise approximation:
 *     - For shallow water: kh ≈ √(k0h)[1 + (k0h/6)(1 + k0h/5)].
 *     - For deeper water: kh is adjusted using an exponential decay so that kh → k0h.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double wu_thornton1986(double k0h)
{
    double threshold = 0.2 * 2 * PI;
    if (k0h <= threshold)
    {
        return sqrt(k0h) * (1 + (k0h / 6.0) * (1 + k0h / 5.0));
    }
    else
    {
        double y = k0h * (1 + 1.26 * exp(-1.84 * k0h));
        return k0h * (1 + 2 * exp(-2 * y) * (1 + exp(-2 * y)));
    }
}

/**
 * Nielsen (1982) approximation for kh.
 *
 * **Formulation:**
 *   For k0h ≤ 2, use a series expansion; for k0h > 2, use an exponential adjustment.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double nielsen1982(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    if (k0h <= 2)
    {
        return sqrt(k0h) * sqrt(1 + (1.0 / 3.0) * k0h + (4.0 / 45.0) * k0h * k0h + (16.0 / 945.0) * pow(k0h, 3));
    }
    else
    {
        return k0h * (1 + 2 * exp(-2 * k0h));
    }
}

/**
 * You solution for shallow water.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double you2002(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    if (k0h <= 2)
    {
        return sqrt(k0h) * sqrt(1 + (1.0 / 3.0) * k0h +
                                (4.0 / 45.0) * k0h * k0h +
                                (16.0 / 945.0) * pow(k0h, 3));
    }
    else
    {
        return k0h * (1 + 2 * exp(-2 * k0h));
    }
}

/**
 * Yu (2014) explicit approximation using trigonometric identity.
 *
 * **Formulation:**
 *   Let a = k0h, then:
 *       kh = a/√(tanh(a)) + 0.0527 · sin( arccos(2·tanh(a)^2.5 – 1) ).
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double yu2014(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double a = k0h;
    double term = 2.0 * pow(tanh(a), 2.5) - 1.0;
    // Ensure term is within the valid range for acos
    term = max(-1.0, min(1.0, term));
    return (a / sqrt(tanh(a))) + 0.0527 * sin(acos(term));
}

/**
 * Gilbert (circa 1989, publ. 2000) empirical approximation (USACE version).
 *
 * **Background:**
 *   A simple curve-fitted formula based on experimental data. Uses a piecewise definition:
 *     - For k0h ≤ 1: kh ≈ √(k0h) · (1 + 0.2·k0h).
 *     - For k0h > 1: kh ≈ k0h · [1 + 0.2·exp(2 – 2·k0h)].
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter.
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double gilbert2000(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    if (k0h <= 1)
    {
        return sqrt(k0h) * (1 + 0.2 * k0h);
    }
    else
    {
        return k0h * (1 + 0.2 * exp(2 - 2 * k0h));
    }
}

/**
 * Guo (2002) explicit solution via logarithmic matching.
 *
 * **Formulation:**
 *   kh = α / [1 – exp(–α^(m/2))]^(1/m), with m ≈ 2.4901 and α = k0h.
 *
 * **Parameters:**
 *   k0h (float): Nondimensional deep-water parameter (α).
 *
 * **Returns:**
 *   float: Approximated nondimensional wavenumber, kh.
 */
double guo2002(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    double m = 2.4901;
    return k0h / pow((1.0 - exp(-pow(k0h, m / 2))), (1.0 / m));
}

/**
 * Guan & Ju (2005) explicit formula.
 *
 * Empirical formula widely used for finite-depth conditions.
 *
 * **Parameters:**
 * k0h Nondimensional deep-water parameter.
 * Approximated nondimensional wavenumber, kh.
 */
double guan2005(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    return sqrt(k0h) * exp(-1.115 * k0h) + k0h * tanh(1.325 * sqrt(k0h));
}

/**
 * Iwagaki (1987) solution.
 *
 * **Parameters:**
 * k0h Nondimensional deep-water parameter (α).
 * Approximated nondimensional wavenumber, kh (β_a).
 */
double iwagaki2007(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }
    return k0h / tanh(sqrt(k0h) * (1.0 + sqrt(k0h) / (2.0 * PI)));
}

/**
 * Eckart (1951) classical approximation for the dispersion relation.
 *
 * **Formulation:**
 *   kh ≈ k0h / √(tanh(k0h))
 *
 * **Parameters:**
 * k0h Nondimensional deep-water parameter.
 * **Returns:**
 * Approximated nondimensional wavenumber, kh.
 */
double eckart1951(double k0h)
{
    if (k0h == 0)
    {
        return 0.0;
    }

    return k0h / sqrt(tanh(k0h));
}

/**
 * Main function implementing demonstration and comparison of various approximation methods.
 */
int main()
{
    // Create test values - use 10000 points like in Python version
    const int num_points = 10000;
    vector<double> k0h_vals;
    k0h_vals.reserve(num_points);

    for (int i = 0; i < num_points; ++i)
    {
        double val = 0.0001 + i * (2 * PI - 0.0001) / (num_points - 1);
        k0h_vals.push_back(val);
    }

    // Define approximation methods to test - match the exact ordering from the Python code
    map<string, function<double(double)>> approximations;

    // Reference exact solution
    approximations["kh_numeric"] = [](double k0h)
    { return kh_numeric(k0h); };

    // Add all the Padé approximants
    for (int i = 1; i <= 13; i++)
    {
        string name = "Pade(2025)_" + to_string(i);
        int formula = i;
        approximations[name] = [formula](double k0h)
        { return pade2025(k0h, formula); };
    }

    // Add all Carvalho methods
    for (int i = 1; i <= 20; i++)
    {
        string name = "Carvalho(2025)_" + to_string(i);
        int formula = i;
        approximations[name] = [formula](double k0h)
        { return carvalho2025(k0h, formula); };
    }

    // Add all Yamaguchi & Nonaka methods
    for (int i = 1; i <= 10; i++)
    {
        string name = "Yamaguchi(2007)_" + to_string(i);
        int formula = i;
        approximations[name] = [formula](double k0h)
        { return YamaguchiNonaka(k0h, formula); };
    }

    // Add all other methods
    approximations["Beji(2013)"] = beji2013;
    approximations["Eckart(1951)"] = eckart1951;
    approximations["Fenton&McKee(1990)_1"] = fenton_mckee1990_1;
    approximations["Fenton&McKee(1990)_2"] = fenton_mckee1990_2;
    approximations["Gilbert(2000)"] = gilbert2000;
    approximations["Guo(2002)"] = guo2002;
    approximations["Guan&Ju(2005)"] = guan2005;
    approximations["Hunt(1979)_5"] = hunt1979_5;
    approximations["Hunt(1979)_9"] = hunt1979_9;
    approximations["Iwagaki(2007)"] = iwagaki2007;
    approximations["Nielsen(1982)"] = nielsen1982;
    approximations["Simarro&Orfila(2013)"] = Simarro_2013;
    approximations["Wu&Thornton(1986)"] = wu_thornton1986;
    approximations["You(2002)"] = you2002;
    approximations["Yu(2014)"] = yu2014;
    approximations["Vatankhah(2013)_1"] = vatankhah2013_1;
    approximations["Vatankhah(2013)_2"] = vatankhah2013_2;

    // Results storage
    vector<tuple<string, double, double, double, double>> results;

    for (const auto &pair : approximations)
    {
        const string &name = pair.first;
        const auto &func = pair.second;

        // Measure execution time
        double random_k0h = static_cast<double>(rand()) / RAND_MAX * 2 * PI;
        auto start = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
        {
            func(random_k0h);
        }
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = end - start;
        double timing = 1000 * duration.count(); // Scaled to match Python time1m

        // Compute errors
        if (name == "kh_numeric")
        {
            // Reference method has zero error by definition
            results.push_back(make_tuple(name, 0.0, 0.0, 0.0, timing));
            continue;
        }

        vector<double> errors;
        errors.reserve(k0h_vals.size());

        for (double val : k0h_vals)
        {
            double exact_kh = kh_numeric(val);
            double approx_kh = func(val);
            double rel_error;

            if (exact_kh == 0)
            {
                rel_error = (approx_kh != 0) ? 100.0 * fabs(approx_kh) : 0.0;
            }
            else
            {
                rel_error = 100.0 * fabs((exact_kh - approx_kh) / exact_kh); // Percent error directly
            }

            errors.push_back(rel_error);
        }

        // Calculate statistics
        double max_error = *max_element(errors.begin(), errors.end());

        // Find k0h value at maximum error
        auto max_it = max_element(errors.begin(), errors.end());
        auto max_idx = static_cast<size_t>(distance(errors.begin(), max_it));
        double k0h_max = k0h_vals[max_idx];

        // Calculate average error
        double sum = 0.0;
        for (const double &e : errors)
        {
            sum += e;
        }
        double mean_error = sum / static_cast<double>(errors.size());

        results.push_back(make_tuple(name, mean_error, max_error, k0h_max, timing));
    }

    // Sort results by average error (ascending)
    sort(results.begin(), results.end(),
         [](const auto &a, const auto &b)
         {
             // First by average error
             if (abs(get<1>(a) - get<1>(b)) > 1e-10)
             {
                 return get<1>(a) < get<1>(b);
             }
             // Then by max error
             return get<2>(a) < get<2>(b);
         });

    // Output header
    cout << "Approximation Errors (absolute %, relative to kh_numeric) for k0h in [0.0001, 2π]" << endl;
    cout << endl;
    cout << "Rank Method                  AvgErr       MaxErr        k0h_MaxErr  Time1M" << endl;

    // Print results - include kh_numeric as the reference (first in rankings)
    int rank = 1;
    for (const auto &result : results)
    {
        string name = get<0>(result);
        double avg_error = get<1>(result);
        double max_error = get<2>(result);
        double k0h_max = get<3>(result);
        double timing = get<4>(result);

        // Formatted output with exact spacing to match the example
        cout << setw(4) << rank << " "
             << left << setw(24) << name
             << right << fixed << setprecision(7) << avg_error << "%   "
             << fixed << setprecision(7) << max_error << "%    "
             << fixed << setprecision(4) << k0h_max << "      "
             << fixed << setprecision(2) << timing << endl;

        rank++;
    }

    // Save to file
    ofstream outfile("wave-disp-equation_output.txt");
    if (outfile.is_open())
    {
        outfile << "Approximation Errors (absolute %, relative to kh_numeric) for k0h in [0.0001, 2π]" << endl;
        outfile << endl;
        outfile << "Rank Method                  AvgErr       MaxErr        k0h_MaxErr  Time1M" << endl;

        rank = 1;
        for (const auto &result : results)
        {
            string name = get<0>(result);
            double avg_error = get<1>(result);
            double max_error = get<2>(result);
            double k0h_max = get<3>(result);
            double timing = get<4>(result);

            // Formatted output
            outfile << setw(4) << rank << " "
                    << left << setw(24) << name
                    << right << fixed << setprecision(7) << avg_error << "%   "
                    << fixed << setprecision(7) << max_error << "%    "
                    << fixed << setprecision(4) << k0h_max << "      "
                    << fixed << setprecision(2) << timing << endl;

            rank++;
        }

        outfile.close();
        cout << "\nResults saved to 'wave-disp-equation_output.txt'" << endl;
    }
    else
    {
        cerr << "Error: Unable to open output file." << endl;
    }

    return 0;
}
