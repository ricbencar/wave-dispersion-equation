"""
============================================================
WAVE DISPERSION EQUATION – Linear Gravity Waves (Airy Theory)
============================================================

This module provides a comprehensive suite of solutions and approximations for analyzing wave dispersion,
essential for wave prediction, oceanographic calculations, and coastal engineering design. It includes a
reference "exact" solution using the Newton-Raphson method, classical and contemporary explicit approximations,
and high-order Padé approximants for high precision.

**Background:**

The **linear wave dispersion equation** relates wave frequency (or period) to wavenumber and water depth
for gravity waves:

    ω² = g · k · tanh(k · h)

where:
    - ω (omega) is the angular frequency (ω = 2π/T, T = wave period),
    - g is the gravitational acceleration,
    - k is the wavenumber (k = 2π/L, L = wavelength),
    - h is the water depth.

The transcendental nature of this equation prevents closed-form solutions for *k*. Therefore, it is
nondimensionalized to:

    k₀h = kh · tanh(kh)

This equation is solved iteratively for the dimensionless wavenumber *kh* (k₀ = ω²/g). An accurate
*kh* evaluation is vital for computing wave phase speed, group velocity, and understanding various
nearshore processes. Explicit approximations bypass the need for iteration but must be chosen
carefully based on accuracy requirements.

**Module Contents:**

  - Reference "Exact" Solution: kh_numeric() implements the Newton-Raphson iteration method for a
    highly precise solution of wave dispersion, acting as a benchmark for other techniques.

  - Classical Approximations: Established methods from researchers like Hunt, Eckart, Nielsen, and Gilbert.

  - Contemporary Approximations: Recent techniques from researchers such as Guo, Beji, Vatankhah &
    Aghashariatmadari, Simarro & Orfila, Yu, Fenton & McKee, Guan & Ju, and Iwagaki.

  - High-Order Padé Approximations: Carvalho's 2025 high-order Padé approximants deliver exceptional
    precision, addressing increasing complexity in wave calculations as a robust alternative to simpler
    methods.

**References:**

 1. Wikipedia. "Airy wave theory". [Wikipedia](https://en.wikipedia.org/wiki/Airy_wave_theory).
 2. Wikipedia. "Dispersion (water waves)". [Wikipedia](https://en.wikipedia.org/wiki/Dispersion_(water_waves)).
 3. Yu, J. (2014). "A Note on Approximations of the Dispersion Relationship of Water Waves",
    *Journal of Engineering Mechanics (ASCE)*, 140(1), 233–237.
 4. Simarro, G. & Orfila, A. (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves:
    Another discussion", *Coastal Engineering*, 80, 15–16.
 5. Vatankhah, A.R. & Aghashariatmadari, Z. (2013). "Improved explicit approximation of linear dispersion relationship
    for gravity waves: a discussion", *Coastal Engineering*, 78, 21–22.
 6. Beji, S. (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves",
    *Coastal Engineering*, 73, 11–12.
 7. **You, Z.J.** (2008). "A close approximation of wave dispersion relation for direct calculations",
    *Applied Ocean Research*, 30(2), 141–143.
 8. Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
    *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
 9. Yamaguchi, M. and H. Nonaka: Comparative study of explicit solutions to wave dispersion equation,
    *Annu. Jour. Eng.*, Ehime Univ., Vol. 6, 2007 in CD-ROM.
10. Carvalho, R. (2006). Unpublished work based on gene expression programming for wave dispersion equations.
11. You, Z.J. "Discussion of 'Simple and explicit solution to the wave dispersion equation'",
    [Coastal Engineering 45 (2002) 71-74], Coastal Eng., Vol. 48, pp.133-135, 2003.
12. Guo, J. (2002). "Simple and explicit solution of the wave dispersion equation",
    *Coastal Engineering*, 45, 71–74.
13. Fenton, J.D. & McKee, W.D. (1990). "On calculating the lengths of water waves",
    *Coastal Engineering*, 14, 499–513.
14. Fenton, J.D. "The numerical solution of steady water wave problems", *Computers & Geosciences*,
    Vol. 4, No. 3, pp.357-368, 1988.
15. Fenton, J.D. (1972). "A ninth-order solution for the solitary wave",
    *Journal of Fluid Mechanics*, 53, 257–271.
16. Wu, C. S. and E. B. Thornton. "Wave numbers of linear progressive waves",
    *Journal of Waterway, Port, Coastal and Ocean Engineering*, ASCE, Vol. 112, No. 4, pp.536-540, 1986.
17. Nielsen, P. "Explicit solutions to practical wave problems", Proc. 19th ICCE, Vol. 1, pp.968-982, 1984.
18. Nielsen, P. "Explicit formulae for practical wave calculations", Coastal Eng., No. 6, pp.389-398, 1982.
19. Hunt, J.N. (1979). "A simple approximation for the dispersion relation of water waves",
    *Journal of Waterway, Port, Coastal and Ocean Engineering*, 105(4), 457–459.
20. Hunt, J. N. "Direct solution of wave dispersion equation", *J. Waterway, Port, Coastal and Ocean Div.*,
    Proc. ASCE, Vol. 105, No. WW4, pp.457-459, 1979.
21. Eckart, C. (1951). "The propagation of gravity waves from deep to shallow water",
    U.S. Department of Commerce, National Bureau of Standards Circular 521.

These references have been arranged by both historical and topical relevance. They span classical
methods (e.g., Eckart, 1951; Hunt, 1979), modern explicit approximations (e.g., Guo, 2002; Beji, 2013;
Vatankhah & Aghashariatmadari, 2013; Simarro & Orfila, 2013; Yu, 2014), as well as pivotal contributions
from Fenton and colleagues (Fenton & McKee, 1990; Fenton, 1972) and innovative computational approaches
(Carvalho, 2006 & 2025) to improve the dispersion relation accuracy.

"""

# Importing the math module to access mathematical functions like sqrt, sinh, tanh, etc.
import math

# Importing numpy for numerical operations and array manipulations.
import numpy as np

# Importing pyplot from matplotlib for data visualization.
import matplotlib.pyplot as plt

# Importing time module to track the execution time of code or create delays.
import time

# =============================================================================
# EXACT SOLUTION (NEWTON–RAPHSON) - Reference Implementation
# =============================================================================

def kh_numeric(k0h, tol=1e-15, max_iter=100):
    """
    Compute the 'exact' nondimensional wavenumber *kh* by numerically solving the dispersion
    relation using the Newton–Raphson method.

    **Equation Solved:**
    The nondimensional dispersion relation (derived from Airy wave theory) is:
        f(kh) = k0h - kh * tanh(kh) = 0
    where k0h = k₀·h (with k₀ = ω²/g) and kh = k·h.

    **Method:**
      - Newton–Raphson iteration:
            kh_new = kh - f(kh) / f'(kh)
      - f'(kh) = -tanh(kh) - kh * sech²(kh)
      - Initial guess: kh₀ ≈ k0h / tanh((6/5)^k0h * sqrt(k0h)) (Carvalho, 2006 style).  This initialization provides
        a reasonable starting point and promotes quicker convergence.

    **Convergence Criteria:**
      Iteratively adjust *kh* until the relative change |Δkh/kh| is below the tolerance `tol`.

    **References:**
      - Fenton & McKee (1990); Yamaguchi & Nonaka (2007); Press et al. (1992).

    **Parameters:**
      k0h (float): Nondimensional deep-water parameter (k₀·h).  Must be non-negative.
      tol (float): Relative convergence tolerance (default: 1e-15). Adjust for desired precision.
      max_iter (int): Maximum number of iterations (default: 100). Increase if convergence fails.

    **Returns:**
      float: Computed nondimensional wavenumber *kh*. Returns 0.0 if `k0h` is 0.

    """
    if k0h == 0:
        return 0.0
    kh = k0h / math.tanh((6/5)**k0h * math.sqrt(k0h))
    for _ in range(max_iter):
        f = k0h - kh * math.tanh(kh)
        df = -math.tanh(kh) - kh / math.cosh(kh)**2
        dkh = f / df
        kh_new = kh - dkh
        if abs(dkh/kh) < tol:
            return kh_new
        kh = kh_new
    return kh

# =============================================================================
# PADE APPROXIMANT - a ratio of two power series
# =============================================================================

def pade2025(k0h, formula):
    """
    Approximations using Padé approximants for the nondimensional wavenumber.

    **Description:**
    Padé approximants are rational functions that approximate a given function by matching its Taylor
    series expansion up to a specified order. They use ratios of polynomials instead of solely
    polynomial expansions as in Taylor series, offering better accuracy, particularly around
    singularities and for functions exhibiting poles. This suite of approximations corresponds to
    approximations derived in Carvalho (2006) using gene expression programming, to fit a ratio
    of power series for efficient wave dispersion estimates.

    **Advantages of Padé Approximants:**
    - Improved Accuracy: Approximations are often superior to Taylor series, mainly in regions where behavior
      cannot be fully explained through polynomial formulations.
    - Convergence Characteristics: Demonstrate accelerated convergence with cases that may be divergent for Taylor series.
    - Singularity Representations: Allows effective characterization close to points involving singular values or poles.
    - Analytic Continuations: It makes provision to broaden approximation to areas beyond standard convergence boundaries with conventional Taylor functions.

    **Error Characteristics:**
      Potential errors can arise from numerical instability in the polynomial evaluation for `k0h > 2π.
      Performance assessments must be conducted for choosing ideal formulation to the given input.

    **Parameters:**
      k0h (float): Nondimensional deep-water parameter (k₀·h). Should be >=0 and <= 2π.
      formula (int): An integer (1 to 13) indicating which formula to compute and use.

    **Returns:**
      float: An approximation to nondimensional wavenumber, kh.
      returns -1.0 when 'formula' parameters are inappropriate / non compliant of input regulations.

    **References:**
        - R. Carvalho (2025). Work published on GitHub, actually the code you're reading right now.
    """

    if formula == 1:
        return ( (1.00649052194019*k0h**0.5 + 0.423646282789217*k0h**1.5 + 0.175406661440005*k0h**2.5) ) / ( (0.306955955676234*k0h**1.0 + 0.0328975279727171*k0h**2.0 + 1) )
    elif formula == 2:
        return ( (0.998980252114366*k0h**0.5 + 0.0240176797055886*k0h**1.5 + 0.102524886754552*k0h**2.5 + 0.0317327085938995*k0h**3.5) ) / ( (-0.150350405960952*k0h**1.0 + 0.112157962910113*k0h**2.0 + 0.00294483072586115*k0h**3.0 + 1) )
    elif formula == 3:
        return ( (1.00006668638419*k0h**0.5 + 0.322645945302282*k0h**1.5 + 0.0860384450810725*k0h**2.5 + 0.051143347041175*k0h**3.5 + 0.0153420957423937*k0h**4.5) ) / ( (0.157166943736625*k0h**1.0 + 0.0245168267924732*k0h**2.0 + 0.0462567432956417*k0h**3.0 + 0.00175392506101448*k0h**4.0 + 1) )
    elif formula == 4:
        return ( (0.999996682596798*k0h**0.5 - 0.0889915717930786*k0h**1.5 + 0.147076211695128*k0h**2.5 + 0.0123471280480147*k0h**3.5 + 0.00866458140843225*k0h**4.5 + 0.00204463718201973*k0h**5.5) ) / ( (-0.255723982020183*k0h**1.0 + 0.159493904911975*k0h**2.0 - 0.0106101311382749*k0h**3.0 + 0.00784491418150148*k0h**4.0 + 0.000184273251439305*k0h**5.0 + 1) )
    elif formula == 5:
        return ( (0.999998218345888*k0h**0.5 - 0.424362176674708*k0h**1.5 + 0.171875463304611*k0h**2.5 - 0.0357487982640122*k0h**3.5 + 0.00410625374333464*k0h**4.5 - 0.000978753693904127*k0h**5.5 - 0.000636955605769902*k0h**6.5) ) / ( (-0.591069429462395*k0h**1.0 + 0.240083348894323*k0h**2.0 - 0.0617593442909405*k0h**3.0 + 0.0104920694265126*k0h**4.0 - 0.00231970889331938*k0h**5.0 - 5.65924775627923e-5*k0h**6.0 + 1) )
    elif formula == 6:
        return ( (1.0000000012405*k0h**0.5 - 0.350251200743747*k0h**1.5 + 0.229153326540668*k0h**2.5 - 0.0205204312544928*k0h**3.5 + 0.0133231478358294*k0h**4.5 + 0.0010401274983046*k0h**5.5 + 0.00048671850792775*k0h**6.5 + 8.40088474488992e-5*k0h**7.5) ) / ( (-0.516917882097882*k0h**1.0 + 0.284751410622371*k0h**2.0 - 0.0555622365621819*k0h**3.0 + 0.0161071584333013*k0h**4.0 - 0.000808341017586247*k0h**5.0 + 0.000381511960690599*k0h**6.0 + 6.40735447518177e-6*k0h**7.0 + 1) )
    elif formula == 7:
        return ( (0.999999995257458*k0h**0.5 - 0.543811114837314*k0h**1.5 + 0.297774393256421*k0h**2.5 - 0.0648661921727468*k0h**3.5 + 0.0174768559302056*k0h**4.5 - 0.00151793039097231*k0h**5.5 + 0.000295750461715408*k0h**6.5 - 4.1567697098083e-6*k0h**7.5 - 1.62498860684328e-5*k0h**8.5) ) / ( (-0.710477969437912*k0h**1.0 + 0.385633880896532*k0h**2.0 - 0.110812474410256*k0h**3.0 + 0.0270500418196004*k0h**4.0 - 0.00394602590941703*k0h**5.0 + 0.000553921270262525*k0h**6.0 - 6.65533846723705e-5*k0h**7.0 - 1.24656671282763e-6*k0h**8.0 + 1) )
    elif formula == 8:
        return ( (1.00000000020126*k0h**0.5 - 0.388439115555858*k0h**1.5 + 0.310223332529737*k0h**2.5 - 0.0496321949056331*k0h**3.5 + 0.0293825301580729*k0h**4.5 - 0.000149900084396432*k0h**5.5 + 0.00139739652490532*k0h**6.5 + 0.000171592322622253*k0h**7.5 + 4.94349555930422e-5*k0h**8.5 + 5.36329658499187e-6*k0h**9.5) ) / ( (-0.555105771884377*k0h**1.0 + 0.372185262898829*k0h**2.0 - 0.0980736559689507*k0h**3.0 + 0.036689986163131*k0h**4.0 - 0.00440840853967443*k0h**5.0 + 0.00139722437126806*k0h**6.0 - 3.64106868131082e-6*k0h**7.0 + 2.72315576473091e-5*k0h**8.0 + 3.80576928768955e-7*k0h**9.0 + 1) )
    elif formula == 9:
        return ( (1.00000000054683*k0h**0.5 - 0.302517970258141*k0h**1.5 + 0.216194173804703*k0h**2.5 - 0.00911867675112112*k0h**3.5 + 0.0131923444114312*k0h**4.5 + 0.00197230469011907*k0h**5.5 + 0.000583685774943952*k0h**6.5 + 0.000117061426950459*k0h**7.5 + 2.16577683514269e-5*k0h**8.5 - 4.06436904033371e-6*k0h**9.5 - 1.35655741732812e-7*k0h**10.5) ) / ( (-0.469184609052338*k0h**1.0 + 0.263835664763448*k0h**2.0 - 0.0421256860069538*k0h**3.0 + 0.0141905686336732*k0h**4.0 + 0.000171234247175628*k0h**5.0 + 0.000284037180073895*k0h**6.0 + 7.23728845968762e-5*k0h**7.0 - 3.33925175755566e-6*k0h**8.0 - 1.1196042731312e-6*k0h**9.0 - 6.10361897619335e-9*k0h**10.0 + 1) )
    elif formula == 10:
        return ( (1.00000000069543*k0h**0.5 - 0.334954381847524*k0h**1.5 + 0.229997336203832*k0h**2.5 - 0.0170932674073579*k0h**3.5 + 0.0142286437913579*k0h**4.5 + 0.00157766182461221*k0h**5.5 + 0.000555414435408452*k0h**6.5 + 0.000109241385068583*k0h**7.5 + 2.04181414663277e-5*k0h**8.5 - 4.47936721320148e-6*k0h**9.5 + 7.20245847805242e-8*k0h**10.5 + 2.27591359161482e-9*k0h**11.5) ) / ( (-0.501621014271701*k0h**1.0 + 0.283044819037831*k0h**2.0 - 0.0523102821262626*k0h**3.0 + 0.0164455383396673*k0h**4.0 - 0.000365622576684072*k0h**5.0 + 0.000304882740985283*k0h**6.0 + 7.15596032735346e-5*k0h**7.0 - 5.64260090297976e-6*k0h**8.0 - 7.22847827757772e-7*k0h**9.0 + 3.16954475149092e-8*k0h**10.0 - 4.18523172095159e-11*k0h**11.0 + 1) )
    elif formula == 11:
        return ( (1.00000000021134*k0h**0.5 - 0.439538511010958*k0h**1.5 + 0.262071966075091*k0h**2.5 - 0.0352260757847662*k0h**3.5 + 0.0136888861362354*k0h**4.5 + 0.00119983590894612*k0h**5.5 + 0.000306132196963262*k0h**6.5 + 9.3344593984067e-5*k0h**7.5 + 1.13485236045952e-5*k0h**8.5 - 2.19033671564094e-6*k0h**9.5 - 1.52303393432862e-7*k0h**10.5 + 3.29680588537e-8*k0h**11.5 - 3.59648926971857e-9*k0h**12.5) ) / ( (-0.606205166733527*k0h**1.0 + 0.332550449569891*k0h**2.0 - 0.0755002415059252*k0h**3.0 + 0.0186169006337531*k0h**4.0 - 0.000624409536143469*k0h**5.0 + 0.000110698204767146*k0h**6.0 + 7.31839474973857e-5*k0h**7.0 - 4.17221912299122e-6*k0h**8.0 - 1.00985652913037e-6*k0h**9.0 + 1.15600808535017e-7*k0h**10.0 - 7.31093937337407e-9*k0h**11.0 - 3.576922099426e-10*k0h**12.0 + 1) )
    elif formula == 12:
        return ( (1.00000000034658*k0h**0.5 - 0.35765423836608*k0h**1.5 + 0.220474157851537*k0h**2.5 - 0.0143101981637556*k0h**3.5 + 0.0106882664277003*k0h**4.5 + 0.00177809831685665*k0h**5.5 + 0.000425551341246499*k0h**6.5 + 8.39454668343144e-5*k0h**7.5 + 1.68918586510378e-5*k0h**8.5 - 3.47736023040874e-6*k0h**9.5 + 1.1564308170262e-7*k0h**10.5 - 5.38824115168048e-9*k0h**11.5 + 2.55400838802905e-9*k0h**12.5 - 5.77246087402385e-10*k0h**13.5) ) / ( (-0.524320887192481*k0h**1.0 + 0.277305164506462*k0h**2.0 - 0.0478782214771677*k0h**3.0 + 0.0124223110237922*k0h**4.0 + 0.00037129622072912*k0h**5.0 + 0.000158246646088704*k0h**6.0 + 5.24876551697954e-5*k0h**7.0 - 1.38557325635002e-6*k0h**8.0 - 7.96232034767084e-7*k0h**9.0 + 8.35464703802028e-9*k0h**10.0 + 1.1834526426828e-8*k0h**11.0 - 1.48873774124982e-9*k0h**12.0 - 5.59846037648501e-11*k0h**13.0 + 1) )
    elif formula == 13:
        return ( (1.00000000043044*k0h**0.5 - 0.341214787680155*k0h**1.5 + 0.216029315236116*k0h**2.5 - 0.0116443516054976*k0h**3.5 + 0.0108812744435703*k0h**4.5 + 0.00184089184824533*k0h**5.5 + 0.000466873470691525*k0h**6.5 + 8.83178779786945e-5*k0h**7.5 + 1.82969284264395e-5*k0h**8.5 - 3.70984687913519e-6*k0h**9.5 + 1.15252770559743e-7*k0h**10.5 + 4.05866020101675e-10*k0h**11.5 + 1.00682972256747e-9*k0h**12.5 - 1.03386709606535e-10*k0h**13.5 + 2.54165346162395e-11*k0h**14.5) ) / ( (-0.507881432557414*k0h**1.0 + 0.27012036160695*k0h**2.0 - 0.0445169260696995*k0h**3.0 + 0.0122182321499095*k0h**4.0 + 0.00039913233839911*k0h**5.0 + 0.000190210573213751*k0h**6.0 + 5.3696821111172e-5*k0h**7.0 - 1.69377956737975e-6*k0h**8.0 - 7.45889615784955e-7*k0h**9.0 + 1.5632250659559e-9*k0h**10.0 + 1.19811786331838e-8*k0h**11.0 - 1.3440816876962e-9*k0h**12.0 + 1.13151579925971e-10*k0h**13.0 + 1.46210486272321e-12*k0h**14.0 + 1) )
    else:
        return -1

# =============================================================================
# CARVALHO (2025) GEP-based approximations
# =============================================================================

def carvalho2025(k0h, formula):
    """
    Approximations using Carvalho's (2025) Gene Expression Programming (GEP) solutions
    for estimating the nondimensional wavenumber *kh*.

    **Description:**
    These approximations, derived via Gene Expression Programming (GEP) (Carvalho, 2006 & 2025),
	represent one of the early attempts to utilize GEP to fit a closed-form,
    explicit expression for the dispersion relation over the whole range of nondimensional water depths.

    **Method:**
    This function employs a suite of pre-computed GEP formulas. Each formula provides a different
    algebraic expression estimating *kh* based on *k0h*. The *formula* parameter selects the
    specific GEP expression to use. Many expressions incorporate combinations of hyperbolic functions,
    power functions, and exponential terms designed to mimic the behavior of the exact dispersion relation.

    **Parameters:**
      k0h (float): Nondimensional deep-water parameter (k₀·h). Must be non-negative.
      formula (int): An integer ranging from 1 to 20, inclusive, selecting which GEP formula to use.

    References:
      - Carvalho, R. (2006). Unpublished work; see also Yamaguchi & Nonaka (2007) for discussion.
      - Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
        *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
      - Ferreira, C. (2006). *Gene Expression Programming: Mathematical Modeling by an Artificial Intelligence*.
        2nd Edition. Springer-Verlag, Germany. This book provides a comprehensive introduction to Gene Expression
        Programming (GEP), a novel genetic algorithm that uses the principles of evolution to solve complex problems
        in mathematics and modeling. Available at: [Springer Link](https://link.springer.com/book/10.1007/3-540-32849-1).
      - Ferreira, C. (2004). "Gene Expression Programming and the Evolution of Computer Programs". 
        In Leandro N. de Castro and Fernando J. Von Zuben (Eds.), *Recent Developments in Biologically Inspired Computing*, 
        pp. 82-103. Idea Group Publishing.
      - Ferreira, C. (2002). "Gene Expression Programming in Problem Solving". In R. Roy, M. Káppen, 
        S. Ovaska, T. Furuhashi, and F. Hoffmann (Eds.), *Soft Computing and Industry: Recent Applications*, 
        pp. 635-654. Springer-Verlag.
      - Ferreira, C. (2002). *Gene Expression Programming: Mathematical Modeling by an Artificial Intelligence*. 
        Angra do Heroismo, Portugal. Online version.
      - Ferreira, C. (2001). "Gene Expression Programming: A New Adaptive Algorithm for Solving Problems", 
        *Complex Systems*, 13(2), 87-129.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if formula == 1:
        kh_Carv = k0h / math.tanh(k0h / (math.sqrt(math.tanh(math.sqrt(math.sinh(k0h)))) * math.tanh(k0h) ** 0.25))
        return (kh_Carv ** 2 + k0h * math.cosh(kh_Carv) ** 2) / (kh_Carv + math.sinh(kh_Carv) * math.cosh(kh_Carv))
    elif formula == 2:
        if k0h <= 1.2:
            return math.sqrt(1 / k0h - math.exp(k0h ** 1.962983 - 6.242035)) / (1 / k0h - 0.168659434)
        elif 1.2 < k0h <= 2.35:
            return (k0h + (k0h / 70.13327717) ** (k0h ** 3)) / math.exp(math.log(4.89859 ** k0h) / (1.134674 - 10 ** k0h))
        elif k0h > 2.35:
            return k0h * math.exp(1.596671172 * k0h / (10 ** k0h))
    elif formula == 3:
        return k0h / math.tanh(
            k0h / (math.sqrt(math.tanh(math.sqrt(math.sinh(k0h)))) * math.tanh(k0h) ** 0.25))
    elif formula == 4:
        return k0h / (math.tanh(k0h / math.tanh(k0h / math.tanh(k0h / math.sinh(math.tanh(math.sqrt(k0h)))))))
    elif formula == 5:
        return k0h / (math.tanh(1.199315 ** (k0h ** 1.047086) * k0h ** 0.499947))
    elif formula == 6:
        return k0h / math.tanh(math.pow(1.1999, k0h ** 1.045) * math.sqrt(k0h))
    elif formula == 7:
        return k0h / math.tanh(k0h / math.tanh(k0h / math.sinh(math.tanh(math.sqrt(k0h)))))
    elif formula == 8:
        return k0h / math.tanh(
            math.sinh(math.sqrt(3.04425 if k0h >= 3.04425 else k0h)) * math.cosh(k0h / 5.194671))
    elif formula == 9:
        return k0h / (math.sqrt(math.tanh(math.sqrt(math.sinh(k0h)))) * math.tanh(k0h) ** 0.25)
    elif formula == 10:
        return k0h / math.tanh(((6 / 5) ** k0h) * math.sqrt(k0h))
    elif formula == 11:
        return k0h / math.tanh(math.sqrt((1.438995 ** k0h) * k0h))
    elif formula == 12:
        return k0h / math.tanh(k0h / math.tanh(math.sinh(math.sqrt(k0h))))
    elif formula == 13:
        return k0h + math.sqrt(k0h) / (4.35144 ** k0h + 0.718409 / (1 / k0h) ** 0.437408)
    elif formula == 14:
        return k0h / (math.tanh(math.sqrt(k0h)) ** (1 / math.cosh(k0h)))
    elif formula == 15:
        return k0h / (math.sqrt(math.tanh(k0h)) * math.tanh(k0h + 1 / math.sqrt(k0h)))
    elif formula == 16:
        return k0h / (math.tanh(k0h) ** ((k0h + 4) / 8))
    elif formula == 17:
        return k0h / ((math.tanh(k0h) ** (k0h / math.tanh(k0h))) ** 0.5)
    elif formula == 18:
        return k0h / math.tanh(math.sinh(math.sqrt(k0h)))
    elif formula == 19:
        return math.sqrt(k0h) + k0h ** 2 / (k0h + 4)
    elif formula == 20:
        return k0h / ((math.sqrt(math.tanh(k0h)) ** (math.tanh(k0h) + 4)) ** 0.25)
    else:
        return -1

# =============================================================================
# YAMAGUCHI & NONANKA (2007) family of explicit solutions
# =============================================================================

def YamaguchiNonaka(k0h, formula):
    """
    Yamaguchi & Nonaka (2007) family of explicit solutions (YN1–YN10).

    Yamaguchi and Nonaka (2007) introduced a set of ten explicit formulas (YN1–YN10) to approximate
    the wave dispersion relation in linear wave theory. These formulas provide accurate and
    computationally efficient alternatives to the implicit dispersion equation, which usually requires
    iterative numerical methods for exact solutions.

    Overview:
    The Yamaguchi & Nonaka (2007) approximations eliminate the need for iterative procedures,
    enhancing computational efficiency in fields like coastal engineering and oceanography.
    Each formula employs mathematical functions such as hyperbolic cotangent (coth), hyperbolic tangent (tanh),
    and exponential functions to approximate the dispersion relation across various depth conditions.

    Parameters:
      k0h (float): Nondimensional deep-water parameter, defined as the product of the wavenumber (k0) and water depth (h).
      formula (int): An integer (1 to 10) specifying which Yamaguchi & Nonaka (2007) formula (YN1–YN10) to use.

    References:
      - Yamaguchi, M. & Nonaka, H. (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation",
        *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
      - Yamaguchi, M. and H. Nonaka: Comparative study of explicit solutions to wave dispersion equation,
        *Annu. Jour. Eng.*, Ehime Univ., Vol. 6, 2007 in CD-ROM.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """

    if k0h == 0:
        return 0.0
    coth = lambda x: math.cosh(x) / math.sinh(x) if x != 0 else float('inf')
    if formula == 1:
        return k0h * coth(k0h**(1.485/2)) ** (1/1.485)
    elif formula == 2:
        return k0h / math.tanh( k0h * (coth(k0h**(1.378/2)))**(1/1.378) )
    elif formula == 3:
        return k0h / math.tanh(math.sqrt(k0h) * (1.0 + math.sqrt(k0h)/(2.0 * math.pi)))
    elif formula == 4:
        return k0h * (1.0 + 1.0/(k0h**2))**0.25
    elif formula == 5:
        return k0h * ((coth(k0h**(1.434/2))) ** (1/1.434))
    elif formula == 6:
        return k0h / math.tanh(math.sqrt(math.sinh(k0h)))
    elif formula == 7:
        return k0h / ((1 - math.exp(-k0h**(2.445/2))) ** (1/2.445))
    elif formula == 8:
        return k0h / math.tanh(k0h * (coth(k0h**(1.310/2)))**(1/1.310))
    elif formula == 9:
        return k0h / math.tanh((1.1965**k0h)*math.sqrt(k0h))
    elif formula == 10:
        return k0h / math.tanh(k0h / (math.sqrt(math.tanh(math.sqrt(math.sinh(k0h)))) * math.tanh(k0h)**0.25))
    else:
        return -1

def Simarro_2013(k0h):
    """
    Simarro & Orfila (2013) two-step Newton-corrected approximation.

    **Concept:**
      Uses Beji’s approximation as the initial guess and then applies one Newton–Raphson correction:
          kh* = [ (kh_B)² + k0h·cosh²(kh_B) ] / [ kh_B + sinh(kh_B)·cosh(kh_B) ],
      where kh_B is the Beji estimate.

    References:
      - Simarro, G. & Orfila, A. (2013). "Improved explicit approximation of linear dispersion relationship
        for gravity waves: Another discussion", *Coastal Engineering*, 80, 15–16.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Corrected approximated nondimensional wavenumber, kh.
    """
    if k0h == 0:
        return 0.0
    kh_Beji = beji2013(k0h)
    return (kh_Beji**2 + k0h * math.cosh(kh_Beji)**2) / (kh_Beji + math.sinh(kh_Beji) * math.cosh(kh_Beji))

def vatankhah2013_1(k0h):
    """
    Vatankhah & Aghashariatmadari (2013) – Single-step explicit formula #2.

    **Concept:**
      Splits the approximation into two parts:
        partA = [k0h + k0h² · exp(–(3.2 + k0h^1.65))] / √(tanh(k0h)),
        partB = k0h · [1 – exp(–k0h^0.132)]^(5.0532 + 2.1584·k0h^1.505).

    Reference:
      - Vatankhah, A.R. & Aghashariatmadari, Z. (2013). "Improved explicit approximation of linear dispersion
        relationship for gravity waves: a discussion", *Coastal Engineering*, 78, 21–22.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h == 0:
        return 0.0
    partA = (k0h + k0h**2 * math.exp(-(3.2 + k0h**1.65))) / math.sqrt(math.tanh(k0h))
    partB = k0h * (1 - math.exp(-k0h**0.132))**(5.0532 + 2.1584*(k0h**1.505))
    return partA + partB

def vatankhah2013_2(k0h):
    """
    Vatankhah & Aghashariatmadari (2013) – Single-step explicit formula #1.

    **Formulation:**
      kh ≈ [k0h + k0h² · exp(–1.835 – 1.225·k0h^1.35)] / √(tanh(k0h)).

    Reference:
      - Vatankhah, A.R. & Aghashariatmadari, Z. (2013). "Improved explicit approximation of linear dispersion
        relationship for gravity waves: a discussion", *Coastal Engineering*, 78, 21–22.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h == 0:
        return 0.0
    return (k0h + k0h**2 * math.exp(-1.835 - 1.225 * k0h**1.35)) / math.sqrt(math.tanh(k0h))

def hunt1979_9(k0h):
    """
    Hunt (1979) Padé-type rational approximation for the dispersion relation.

    **Purpose:**
      Provides an explicit expression for kh by approximating tanh(kh) via a rational (Padé) function.

    **Formulation:**
      kh ≈ √[ k0h² + k0h / (1 + c₁·k0h + c₂·k0h² + … + c₉·k0h⁹) ],
      where the coefficients c₁,…,c₉ are from Hunt (1979).

    References:
      - Hunt, J.N. (1979). "A simple approximation for the dispersion relation of water waves",
        *Journal of Waterway, Port, Coastal and Ocean Engineering*, 105(4), 457–459.
      - Hunt, J. N.: Direct solution of wave dispersion equation, J. Waterway, Port, Coastal and Ocean Div.,
        Proc. ASCE, Vol. 105, No. WW4, pp.457-459, 1979.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    coeffs = [0.6666666667, 0.3555, 0.16084, 0.0632, 0.02174,
              0.00654, 0.00171, 0.00039, 0.00011]
    s = sum(coeffs[i] * (k0h ** (i+1)) for i in range(len(coeffs)))
    return math.sqrt(k0h**2 + k0h / (1 + s))

def hunt1979_5(k0h):
    """
    Hunt (1979) – 5th-order approximate solution (Hunt1) from Yamaguchi & Nonaka (2007).

    **Formulation:**
      (kₐ·h)² = α · [ α + 1 / (1 + 0.6522·α + 0.4622·α² + 0.0864·α⁴ + 0.0675·α⁵) ],
      where α = k0h.

    References:
      - Hunt, J.N. (1979). "A simple approximation for the dispersion relation of water waves",
        *Journal of Waterway, Port, Coastal and Ocean Engineering*, 105(4), 457–459.
      - Hunt, J. N.: Direct solution of wave dispersion equation, J. Waterway, Port, Coastal and Ocean Div.,
        Proc. ASCE, Vol. 105, No. WW4, pp.457-459, 1979.

    Parameters:
      k0h (float): Nondimensional deep-water parameter (α).

    Returns:
      float: Approximated nondimensional wavenumber (kₐ·h).
    """
    alpha = k0h
    denom = 1.0 + 0.6522*alpha + 0.4622*(alpha**2) + 0.0864*(alpha**4) + 0.0675*(alpha**5)
    return math.sqrt(alpha * (alpha + 1.0/denom))

def fenton_mckee1990_1(k0h):
    """
    Fenton & McKee (1990) iterative-type approximation for kh, as described in Yamaguchi & Nonaka (2007).

    **Formulation:**
      Compute βₐ = k0h · [coth(k0h)]^(1/2), then
          kh = [ k0h + βₐ² · sech²(βₐ) ] / [ tanh(βₐ) + βₐ · sech²(βₐ) ].

    References:
      - Fenton, J.D. & McKee, W.D. (1990). "On calculating the lengths of water waves",
        *Coastal Engineering*, 14, 499–513.
      - Fenton, J.D.: The numerical solution of steady water wave problems, Computers & Geosciences,
        Vol. 4, No. 3, pp.357-368, 1988.
      - Wu, C. S. and E. B. Thornton: Wave numbers of linear progressive waves,
        *Journal of Waterway, Port, Coastal and Ocean Engineering*, ASCE, Vol. 112, No. 4, pp.536-540, 1986.
      - Fenton, J.D. (1972). "A ninth-order solution for the solitary wave",
        *Journal of Fluid Mechanics*, 53, 257–271.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximate nondimensional wavenumber, kh.
    """
    alpha = k0h
    beta_a = alpha * (math.cosh(alpha) / math.sinh(alpha))**0.5
    sech_beta_a = 1.0 / math.cosh(beta_a)
    tanh_beta_a = math.tanh(beta_a)
    beta_a_sq = beta_a * beta_a
    sech_sq = sech_beta_a * sech_beta_a
    numerator = alpha + beta_a_sq * sech_sq
    denominator = tanh_beta_a + beta_a * sech_sq
    return numerator / denominator

def fenton_mckee1990_2(k0h):
    """
    Fenton & McKee (1990) all-depth empirical approximation for kh.

    **Formulation:**
      kh ≈ k0h / [ tanh(k0h^(3/4)) ]^(2/3).

    References:
      - Fenton, J.D. & McKee, W.D. (1990). "On calculating the lengths of water waves",
        *Coastal Engineering*, 14, 499–513.
      - Fenton, J.D.: The numerical solution of steady water wave problems, Computers & Geosciences,
        Vol. 4, No. 3, pp.357-368, 1988.
      - Wu, C. S. and E. B. Thornton: Wave numbers of linear progressive waves,
        *Journal of Waterway, Port, Coastal and Ocean Engineering*, ASCE, Vol. 112, No. 4, pp.536-540, 1986.
      - Fenton, J.D. (1972). "A ninth-order solution for the solitary wave",
        *Journal of Fluid Mechanics*, 53, 257–271.

    Parameters:
      k0h (float): Deep-water parameter.

    Returns:
      float: Approximate nondimensional wavenumber, kh.
    """
    return k0h / (math.tanh(k0h**(3/4))**(2/3))

def wu_thornton1986(k0h):
    """
    Wu & Thornton (1986) explicit approximation for the dispersion relation.

    **Overview:**
      Provides a piecewise approximation:
        - For shallow water: kh ≈ √(k0h)[1 + (k0h/6)(1 + k0h/5)].
        - For deeper water: kh is adjusted using an exponential decay so that kh → k0h.

    References:
      - Wu, C. S. and E. B. Thornton: Wave numbers of linear progressive waves,
        J. Waterway, Port, Coastal, and Ocean Eng., ASCE, Vol. 112, No. 4, pp.536-540, 1986.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h <= 0.2 * 2 * math.pi:
        return math.sqrt(k0h) * (1 + (k0h/6) * (1 + k0h/5))
    else:
        y = k0h * (1 + 1.26 * math.exp(-1.84 * k0h))
        return k0h * (1 + 2 * math.exp(-2*y) * (1 + math.exp(-2*y)))

def beji2013(k0h):
    """
    Beji (2013) improved explicit approximation.

    **Formulation:**
      kh ≈ [k0h/√(tanh(k0h))] · [1 + k0h^1.09 · exp(–1.55 – 1.30·k0h – 0.216·k0h²)].

    Reference:
      - Beji, S. (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves",
        *Coastal Engineering*, 73, 11–12.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h == 0:
        return 0.0
    return (k0h * (1 + k0h**1.09 * math.exp(-(1.55 + 1.30*k0h + 0.216*k0h**2)))) / math.sqrt(math.tanh(k0h))

def nielsen1982(k0h):
    """
    Nielsen (1982) approximation for kh.

    **Formulation:**
      For k0h ≤ 2, use a series expansion; for k0h > 2, use an exponential adjustment.

    References:
    - Nielsen, P.: Explicit formulae for practical wave calculations, Coastal Eng., No. 6, pp.389-398, 1982.
    - Nielsen, P.: Explicit solutions to practical wave problems, Proc. 19th ICCE, Vol. 1, pp.968-982, 1984.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h <= 2:
        return math.sqrt(k0h) * math.sqrt(1 + (1/3)*k0h + (4/45)*k0h**2 + (16/945)*k0h**3)
    else:
        return k0h * (1 + 2 * math.exp(-2 * k0h))

def you2002(k0h):
    """
    You solution for shallow water (from Yamaguchi & Nonaka (2007), Eq. (66)).

    Reference:
      - **You, Z.J.** (2008). "A close approximation of wave dispersion relation for direct calculations",
        *Applied Ocean Research*, 30(2), 141–143.
      - You, Z.J.: Discussion of “Simple and explicit solution to the wave dispersion equation”,
        [Coastal Engineering 45 (2002) 71-74], Coastal Eng., Vol. 48, pp.133-135, 2003.

    Parameters:
      k0h : Nondimensional deep-water parameter (α).

    Returns:
      Approximated nondimensional wavenumber, kh.
    """
    if k0h <= 2:
        return math.sqrt(k0h) * math.sqrt(1 + (1/3) * k0h +
                                          (4/45) * k0h**2 +
                                          (16/945) * k0h**3)
    else:
        return k0h * (1 + 2 * math.exp(-2 * k0h))

def yu2014(k0h):
    """
    Yu (2014) explicit approximation using trigonometric identity.

    **Formulation:**
      Let a = k0h, then:
          kh = a/√(tanh(a)) + 0.0527 · sin( arccos(2·tanh(a)^2.5 – 1) ).

    Reference:
      - Yu, J. (2014). "A Note on Approximations of the Dispersion Relationship of Water Waves",
        *Journal of Engineering Mechanics (ASCE)*, 140(1), 233–237.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    a = k0h
    if a == 0:
        return 0.0
    term = 2.0 * (math.tanh(a) ** 2.5) - 1.0
    return (a / math.sqrt(math.tanh(a))) + 0.0527 * math.sin(math.acos(term))

def gilbert2000(k0h):
    """
    Gilbert (circa 1989, publ. 2000) empirical approximation (USACE version).

    **Background:**
      A simple curve-fitted formula based on experimental data. Uses a piecewise definition:
        - For k0h ≤ 1: kh ≈ √(k0h) · (1 + 0.2·k0h).
        - For k0h > 1: kh ≈ k0h · [1 + 0.2·exp(2 – 2·k0h)].

    Reference:
      - USACE Technical Note CETN-I-17; Soulsby (2006); Wiberg & Sherwood (2008).

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h <= 1:
        return math.sqrt(k0h) * (1 + 0.2 * k0h)
    else:
        return k0h * (1 + 0.2 * math.exp(2 - 2 * k0h))

def guo2002(k0h):
    """
    Guo (2002) explicit solution via logarithmic matching.

    **Formulation:**
      kh = α / [1 – exp(–α^(m/2))]^(1/m), with m ≈ 2.4901 and α = k0h.

    Reference:
      - Guo, J. (2002). "Simple and explicit solution of the wave dispersion equation",
        *Coastal Engineering*, 45, 71–74..

    Parameters:
      k0h (float): Nondimensional deep-water parameter (α).

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    if k0h == 0:
        return 0.0
    m = 2.4901
    return k0h / ((1.0 - math.exp(-k0h**(m/2))) ** (1.0/m))

def guan2005(k0h):
    """
    Guan & Ju (2005) explicit formula.

    Empirical formula widely used for finite-depth conditions.

    Reference:
      - Guan, C. & Ju, H. (2005), "Empirical formula for wavelength of ocean wave in finite depth water",
        Chinese Journal of Oceanology and Limnology, Vol.23(1), pp.17–22.

    Parameters:
      k0h : Nondimensional deep-water parameter.

    Returns:
      Approximated nondimensional wavenumber, kh.
    """
    return math.sqrt(k0h) * math.exp(-1.115 * k0h) + k0h * math.tanh(1.325 * math.sqrt(k0h))

def iwagaki2007(k0h):
    """
    Iwagaki (1987) solution [Eq. (13) in Yamaguchi & Nonaka (2007)].

    The original formula is:
       β_a = α * coth{ sqrt(α) * [ 1 + sqrt(α)/(2π) ] }
    or equivalently
       β_a = α / tanh{ sqrt(α) * [ 1 + sqrt(α)/(2π) ] }

    Parameters:
      k0h : Nondimensional deep-water parameter (α).

    Returns:
      Approximated nondimensional wavenumber, kh (β_a).
    """
    return k0h / math.tanh( math.sqrt(k0h) * (1.0 + math.sqrt(k0h)/(2.0*math.pi)) )

def eckart1951(k0h):
    """
    Eckart (1951/1990) early explicit approximation.

    Reference:
      - Eckart, C. (1951). "The propagation of gravity waves from deep to shallow water",
        U.S. Department of Commerce, National Bureau of Standards Circular 521.

    Parameters:
      k0h (float): Nondimensional deep-water parameter.

    Returns:
      float: Approximated nondimensional wavenumber, kh.
    """
    return k0h / (math.tanh(k0h)**0.5)

# =============================================================================
# ORDERED APPROXIMATIONS (ALL FORMULAS ARE RANKED)
# =============================================================================
ordered_approx = {
    "kh_numeric": kh_numeric,

    "Pade(2025)_1": lambda x: pade2025(x, 1),
    "Pade(2025)_2": lambda x: pade2025(x, 2),
    "Pade(2025)_3": lambda x: pade2025(x, 3),
    "Pade(2025)_4": lambda x: pade2025(x, 4),
    "Pade(2025)_5": lambda x: pade2025(x, 5),
    "Pade(2025)_6": lambda x: pade2025(x, 6),
    "Pade(2025)_7": lambda x: pade2025(x, 7),
    "Pade(2025)_8": lambda x: pade2025(x, 8),
    "Pade(2025)_9": lambda x: pade2025(x, 9),
    "Pade(2025)_10": lambda x: pade2025(x, 10),
    "Pade(2025)_11": lambda x: pade2025(x, 11),
    "Pade(2025)_12": lambda x: pade2025(x, 12),
    "Pade(2025)_13": lambda x: pade2025(x, 13),

    "Carvalho(2025)_1": lambda x: carvalho2025(x, 1),
    "Carvalho(2025)_2": lambda x: carvalho2025(x, 2),
    "Carvalho(2025)_3": lambda x: carvalho2025(x, 3),
    "Carvalho(2025)_4": lambda x: carvalho2025(x, 4),
    "Carvalho(2025)_5": lambda x: carvalho2025(x, 5),
    "Carvalho(2025)_6": lambda x: carvalho2025(x, 6),
    "Carvalho(2025)_7": lambda x: carvalho2025(x, 7),
    "Carvalho(2025)_8": lambda x: carvalho2025(x, 8),
    "Carvalho(2025)_9": lambda x: carvalho2025(x, 9),
    "Carvalho(2025)_10": lambda x: carvalho2025(x, 10),
    "Carvalho(2025)_11": lambda x: carvalho2025(x, 11),
    "Carvalho(2025)_12": lambda x: carvalho2025(x, 12),
    "Carvalho(2025)_13": lambda x: carvalho2025(x, 13),
    "Carvalho(2025)_14": lambda x: carvalho2025(x, 14),
    "Carvalho(2025)_15": lambda x: carvalho2025(x, 15),
    "Carvalho(2025)_16": lambda x: carvalho2025(x, 16),
    "Carvalho(2025)_17": lambda x: carvalho2025(x, 17),
    "Carvalho(2025)_18": lambda x: carvalho2025(x, 18),
    "Carvalho(2025)_19": lambda x: carvalho2025(x, 19),
    "Carvalho(2025)_20": lambda x: carvalho2025(x, 20),

    "Yamaguchi(2007)_1": lambda x: YamaguchiNonaka(x, 1),
    "Yamaguchi(2007)_2": lambda x: YamaguchiNonaka(x, 2),
    "Yamaguchi(2007)_3": lambda x: YamaguchiNonaka(x, 3),
    "Yamaguchi(2007)_4": lambda x: YamaguchiNonaka(x, 4),
    "Yamaguchi(2007)_5": lambda x: YamaguchiNonaka(x, 5),
    "Yamaguchi(2007)_6": lambda x: YamaguchiNonaka(x, 6),
    "Yamaguchi(2007)_7": lambda x: YamaguchiNonaka(x, 7),
    "Yamaguchi(2007)_8": lambda x: YamaguchiNonaka(x, 8),
    "Yamaguchi(2007)_9": lambda x: YamaguchiNonaka(x, 9),
    "Yamaguchi(2007)_10": lambda x: YamaguchiNonaka(x, 10),

    "Beji(2013)": beji2013,
    "Eckart(1951)": eckart1951,
    "Fenton&McKee(1990)_1": fenton_mckee1990_1,
    "Fenton&McKee(1990)_2": fenton_mckee1990_2,
    "Gilbert(2000)": gilbert2000,
    "Guo(2002)": guo2002,
    "Guan&Ju(2005)": guan2005,
    "Hunt(1979)_5": hunt1979_5,
    "Hunt(1979)_9": hunt1979_9,
    "Iwagaki(2007)": iwagaki2007,
    "Nielsen(1982)": nielsen1982,
    "Simarro&Orfila(2013)": Simarro_2013,
    "Wu&Thornton(1986)": wu_thornton1986,
    "You(2002)": you2002,
    "Yu(2014)": yu2014,
    "Vatankhah(2013)_1": vatankhah2013_1,
    "Vatankhah(2013)_2": vatankhah2013_2,
}

# =============================================================================
# ERROR ANALYSIS
# =============================================================================
#
# Approximation Errors (absolute %, relative to kh_numeric) for k0h in [0.0001, 2π]
#
# Rank    Method                  AvgErr       MaxErr        k0h_MaxErr    Time1M
#    1    kh_numeric              0.0000000%   0.0000000%    0.0001        5.13
#    2    Pade(2025)_11           0.0000000%   0.0000000%    0.0001        4.82
#    3    Pade(2025)_8            0.0000000%   0.0000000%    0.0001        3.63
#    4    Pade(2025)_12           0.0000000%   0.0000000%    0.0001        5.15
#    5    Pade(2025)_13           0.0000000%   0.0000000%    0.0001        5.72
#    6    Pade(2025)_9            0.0000000%   0.0000001%    0.0001        4.01
#    7    Pade(2025)_10           0.0000000%   0.0000001%    0.0001        4.48
#    8    Pade(2025)_7            0.0000000%   0.0000075%    5.3306        3.27
#    9    Pade(2025)_6            0.0000000%   0.0000002%    0.7686        2.93
#   10    Carvalho(2025)_1        0.0000003%   0.0000042%    0.3256        2.08
#   11    Simarro&Orfila(2013)    0.0000005%   0.0000082%    0.0359        2.00
#   12    Pade(2025)_4            0.0000136%   0.0003311%    0.0001        2.36
#   13    Pade(2025)_5            0.0000186%   0.0064444%    2.9490        2.57
#   14    Vatankhah(2013)_1       0.0004901%   0.0017590%    0.9515        1.62
#   15    Pade(2025)_3            0.0006188%   0.0066566%    0.0001        1.77
#   16    Carvalho(2025)_2        0.0028329%   0.0167221%    1.2003        0.68
#   17    Hunt(1979)_9            0.0035762%   0.0081503%    3.7892        3.49
#   18    Fenton&McKee(1990)_1    0.0055981%   0.0511137%    0.4393        1.09
#   19    Pade(2025)_2            0.0059286%   0.1018976%    0.0001        1.40
#   20    Wu&Thornton(1986)       0.0069848%   0.0342738%    0.7089        0.87
#   21    Vatankhah(2013)_2       0.0085449%   0.0188902%    0.0705        0.96
#   22    Carvalho(2025)_3        0.0094391%   0.0444630%    1.1274        1.35
#   23    Yamaguchi(2007)_10      0.0094391%   0.0444630%    1.1274        1.64
#   24    Carvalho(2025)_4        0.0097978%   0.0501488%    0.3463        1.19
#   25    Beji(2013)              0.0164113%   0.0442262%    0.2998        1.02
#   26    Carvalho(2025)_5        0.0293525%   0.0760247%    1.5603        0.91
#   27    Carvalho(2025)_6        0.0296085%   0.0654110%    0.5769        0.94
#   28    Carvalho(2025)_7        0.0308926%   0.1351853%    1.3863        1.11
#   29    Hunt(1979)_5            0.0385834%   0.0783780%    1.8098        0.89
#   30    Nielsen(1982)           0.0451224%   0.5804859%    1.9996        0.37
#   31    You(2002)               0.0451224%   0.5804859%    1.9996        0.38
#   32    Carvalho(2025)_8        0.0642401%   0.1698012%    2.6550        1.03
#   33    Carvalho(2025)_9        0.0668541%   0.2037817%    2.6569        1.30
#   34    Carvalho(2025)_10       0.0729712%   0.2712270%    0.3941        0.85
#   35    Carvalho(2025)_11       0.0742522%   0.2609919%    0.3872        0.97
#   36    Yu(2014)                0.0863942%   0.3291155%    0.0843        1.10
#   37    Yamaguchi(2007)_9       0.0883941%   0.3122414%    1.6590        0.95
#   38    Carvalho(2025)_12       0.0892311%   0.4009798%    0.9621        1.06
#   39    Pade(2025)_1            0.1226442%   0.6485218%    0.0001        1.07
#   40    Carvalho(2025)_13       0.1267865%   0.3714239%    0.0007        1.00
#   41    Yamaguchi(2007)_2       0.1888776%   0.7321854%    0.1792        1.26
#   42    Carvalho(2025)_14       0.2085362%   1.2784690%    0.2200        1.12
#   43    Gilbert(2000)           0.2292764%   0.7930431%    0.4858        0.41
#   44    Yamaguchi(2007)_8       0.2554179%   1.1092228%    1.0287        1.37
#   45    Guo(2002)               0.2863639%   0.7571934%    1.7822        0.73
#   46    Carvalho(2025)_15       0.2925759%   1.3110942%    0.1421        1.17
#   47    Carvalho(2025)_16       0.3101789%   1.4216558%    0.1365        0.98
#   48    Guan&Ju(2005)           0.3393654%   0.4970393%    3.5938        0.67
#   49    Carvalho(2025)_17       0.3571345%   1.8310303%    0.3281        1.21
#   50    Yamaguchi(2007)_7       0.3597263%   1.0313632%    1.6848        0.93
#   51    Carvalho(2025)_18       0.3724391%   1.1289804%    1.4912        1.03
#   52    Carvalho(2025)_19       0.4540868%   1.3892626%    0.3067        0.99
#   53    Carvalho(2025)_20       0.5191037%   2.0235808%    1.1161        1.51
#   54    Fenton&McKee(1990)_2    0.6716352%   1.6314251%    0.3400        0.56
#   55    Yamaguchi(2007)_1       0.7077884%   1.5428393%    1.9814        1.05
#   56    Yamaguchi(2007)_5       0.8704391%   2.0361583%    1.8695        1.14
#   57    Eckart(1951)            1.1926080%   4.9797834%    0.6938        0.37
#   58    Yamaguchi(2007)_6       1.2101780%   5.2587798%    0.7403        0.84
#   59    Yamaguchi(2007)_3       1.4889834%   3.1473467%    1.8060        0.90
#   60    Iwagaki(2007)           1.4889834%   3.1473467%    1.8060        0.62
#   61    Yamaguchi(2007)_4       1.5872652%   3.1772841%    0.4268        0.73

def compute_errors(func, k0h_vals):
    """
    For a given approximation function 'func' and an array of k0h values,
    compute the average and maximum absolute relative errors (%) compared
    to the exact solution 'kh_numeric'.

    Absolute relative error (%) = 100 * |approx - exact| / |exact|

    Returns:
      (average_error, maximum_error, k0h_max)
      where k0h_max is the k0h value at which the maximum error occurs.
    """
    errors = []
    for val in k0h_vals:
        exact_kh = kh_numeric(val)
        approx_kh = func(val)
        if exact_kh == 0:
            err = 0.0
        else:
            err = 100.0 * abs(approx_kh - exact_kh) / abs(exact_kh)
        errors.append(err)
    avg_error = np.mean(errors)
    max_error = max(errors)
    idx_max = errors.index(max_error)
    k0h_max = k0h_vals[idx_max]
    return avg_error, max_error, k0h_max

# =============================================================================
# MAIN DRIVER
# =============================================================================
def main():
    """
    Main driver function to:
      - Evaluate each dispersion approximation method over k0h in [0.0001, 2π] using 10,000 points.
      - Compare each method's output to the exact solution computed by kh_numeric.
      - Compute for each method:
          * Average relative error (%) and maximum relative error (%).
          * The k0h value where the maximum error occurs ("k0h_max").
          * The time (in seconds, rounded to 1 decimal) it takes to compute the approximation
            1 million times at a random k0h value ("Time1M").
      - Sort methods by average relative error (primary) and maximum relative error (secondary).
      - Print a ranking table with detailed error statistics and write it to wave-disp-equation_output.txt.
      - Plot a chart comparing the average errors of all approximation methods.
    """
    k0h_vals = np.linspace(0.0001, 2 * math.pi, 10000)
    results = []
    for name, func in ordered_approx.items():
        try:
            avg_e, max_e, k0h_max = compute_errors(func, k0h_vals)
        except Exception as exc:
            avg_e, max_e, k0h_max = None, None, None
            print(f"Error computing {name}: {exc}")
        try:
            random_k0h = np.random.uniform(0.0001, 2*math.pi)
            start = time.perf_counter()
            for _ in range(1000):
                func(random_k0h)
            end = time.perf_counter()
            time1m = 1000*(end - start)
        except Exception as exc:
            time1m = None
            print(f"Error timing {name}: {exc}")
        results.append((name, avg_e, max_e, k0h_max, time1m))

    results_sorted = sorted(
        results,
        key=lambda x: (float('inf') if x[1] is None else x[1],
                       float('inf') if x[2] is None else x[2])
    )

    header = "Approximation Errors (absolute %, relative to kh_numeric) for k0h in [0.0001, 2π]\n\n"
    header += f"{'Rank':4s} {'Method':17s} {'AvgErr':>12s} {'MaxErr':>12s} {'k0h_MaxErr':>17s} {'Time1M':>7s}"
    lines = [header]
    for idx, (meth, av, mx, k0h_max, t1m) in enumerate(results_sorted, start=1):
        if av is None:
            line = (
                f"{idx:4d} "
                f"{meth:20s} "
                f"{'ERROR':>12s}"
                f"{'':>13s} "
                f"{'':>10s} "
                f"{'':>10s}"
            )
        else:
            line = (
                f"{idx:4d} "
                f"{meth:20s} "
                f"{av:12.7f}%"
                f"{mx:12.7f}%"
                f"{k0h_max:10.4f}"
                f"{t1m:10.2f}"
            )
        lines.append(line)
    output_str = "\n".join(lines)
    print(output_str)
    with open("wave-disp-equation_output.txt", "w", encoding="utf-8") as f:
        f.write(output_str)

    # =============================================================================
    # PLOTTING: Plot a chart of average absolute error (%) for each method
    # =============================================================================

    # Filter out methods with valid error values
    valid_results = [r for r in results_sorted if r[1] is not None]
    methods = [r[0] for r in valid_results]
    avg_errors = [r[1] for r in valid_results]

    # Create a figure
    plt.figure(figsize=(16, 8))

    # Create bar plot
    bars = plt.bar(methods, avg_errors, color='skyblue', edgecolor='black')

    # Set the y-axis major ticks to have a step of 0.1
    plt.yticks(np.arange(0, max(avg_errors)+0.1, 0.1))

    # Add labels and title
    plt.xlabel('Approximation Method', fontsize=12)
    plt.ylabel('Average Absolute Relative Error (%)', fontsize=12)
    plt.title('Average Absolute Relative Errors (%) of Wave Dispersion Approximations', fontsize=12)

    # Enhance the x-tick labels
    plt.xticks(rotation=90, fontsize=11)

    # Set the limits for x-axis to occupy the whole space
    plt.xlim(-0.7, len(methods) - 0.0)

    # Optional: Add gridlines for better readability
    plt.grid(axis='both', which='major', linestyle='--', linewidth=1.5, alpha=0.7)  # Main x and y gridlines
	
    # Adding vertical data value labels on top of the bars
    for i, bar in enumerate(bars):
        yval = bar.get_height()
        # Choose the format based on the index
        if i < 56:
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.01, f"{yval:.15f}%", ha='center', va='bottom', rotation=90,color='red')
        else:
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.01, f"{yval:.3f}%", ha='center', va='bottom', rotation=90,color='red')

    # Text to be added at the top
    explanatory_text = (
        "The linear wave dispersion equation relates wave frequency (or period) to wavenumber and water depth\n\n"
        "\n\n\nwhere:\n\n"
        "    - ω (omega) is the angular frequency (ω = 2π/T, T = wave period),\n"
        "    - k is the wavenumber (k = 2π/L, L = wavelength),\n"
        "    - index 0 = offshore conditions (k₀ = 2π/L₀, L₀ = gT²/(2π)),\n"
        "    - h is the water depth,\n"
        "    - g is the gravitational acceleration."
    )

    # Adding the explanatory text at the top of the chart
    plt.text(0.45, 0.95, explanatory_text, ha='center', va='top', fontsize=12, fontweight='normal', transform=plt.gca().transAxes)

    # Now adding the bold equation separately
    plt.text(0.5, 0.90, "\nω² = g · k · tanh(k · h)   or   k₀ · h = k · h · tanh(k · h)", ha='center', va='top', fontsize=14, fontweight='bold', transform=plt.gca().transAxes)

    # Adjust the margins to reduce blank space
    plt.subplots_adjust(left=0.05, right=0.95)

    # Ensure layout is tight
    plt.tight_layout()

    # Save the figure
    plt.savefig("wave-disp-equation_errors.png", dpi=300, bbox_inches='tight')

    # Show the plot
    plt.show()

if __name__ == "__main__":
    main()
