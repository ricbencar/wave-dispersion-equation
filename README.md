# WAVE DISPERSION EQUATION: Theoretical and Computational Documentation 

## Introduction to the Computational Framework 

**The repository wave-disp-equation.py represents a computational suite engineered to solve this specific physical relationship through a spectrum of methodologies, ranging from exact numerical iteration to state-of-the-art analytical approximations.** Because the fundamental relationship governing surface gravity waves is transcendental, the extraction of the wavenumber from a known wave period and water depth introduces a significant computational bottleneck in large-scale phase-resolving and spectral wave models.

This documentation serves as a comprehensive theoretical and mathematical companion to the `wave-disp-equation.py` module. It systematically deconstructs every algorithm, explicit formulation, and Padé approximant codified within the repository, explaining the hydrodynamic origins, mathematical derivations, and computational optimization strategies. The analysis traverses classical historical approximations, nested minimax optimizations, restricted domain expansions, and the highly advanced half-integer fractional Padé approximants. 

---

## Theoretical Foundations: Linear Gravity Waves (Airy Theory) 
To fully contextualize the algorithms within the `wave-disp-equation.py` repository, it is imperative to establish the fluid dynamic principles from which the governing equations are derived. 
**Airy wave theory, formally recognized as linear wave theory, provides the foundational mathematical description of surface gravity wave propagation over a homogeneous fluid layer.** 

### Fluid Dynamic Assumptions and the Laplace Equation 
The derivation of the linear dispersion relation begins with the Navier-Stokes equations, subjected to a specific set of simplifying assumptions. The fluid is assumed to be inviscid (frictionless), incompressible (constant density), and undergoing irrotational flow. The assumption of irrotationality ($\nabla \times \mathbf{u} = 0$) allows the fluid velocity vector $\mathbf{u} = (u, w)$ to be expressed as the gradient of a scalar velocity potential $\Phi(x, z, t)$ such that $\mathbf{u} = \nabla \Phi$. 
Due to the incompressibility constraint ($\nabla \cdot \mathbf{u} = 0$), the velocity potential must satisfy the Laplace equation throughout the entire fluid domain:

$$\frac{\partial^2 \Phi}{\partial x^2} + \frac{\partial^2 \Phi}{\partial z^2} = 0$$ 

While potential flow theory typically fails in boundary-dominated fluid mechanics, it succeeds exceptionally well for surface gravity waves because wave-induced vorticity is strictly confined to microscopic oscillatory Stokes boundary layers at the seabed and the free surface, leaving the vast majority of the water column irrotational. 

### Boundary Value Problem and Linearization 
The Laplace equation is a linear, elliptic partial differential equation that requires strict boundary conditions to close the system. The fluid domain is bounded horizontally but vertically constrained by an impermeable flat seabed at $z = -h$ and a free-moving fluid surface at $z = \eta(x,t)$. 
* **Kinematic Bed Boundary Condition:** Fluid particles cannot penetrate the solid seabed. Therefore, the vertical velocity at the bottom must be zero:
 $$\left. \frac{\partial \Phi}{\partial z} \right|_{z = -h} = 0$$ 
* **Kinematic Free-Surface Boundary Condition:** A fluid particle that is on the free surface must remain on the free surface. The vertical velocity of the surface dictates the vertical velocity of the fluid:
 $$\frac{\partial \eta}{\partial t} + \frac{\partial \Phi}{\partial x}\frac{\partial \eta}{\partial x} = \frac{\partial \Phi}{\partial z} \quad \text{at} \quad z = \eta(x,t)$$ 
* **Dynamic Free-Surface Boundary Condition:** The pressure at the free surface must equal the uniform atmospheric pressure. AApplying the unsteady Bernoulli equation yields:
$$\frac{\partial \Phi}{\partial t} + \frac{1}{2}\left[ \left(\frac{\partial \Phi}{\partial x}\right)^2 + \left(\frac{\partial \Phi}{\partial z}\right)^2 \right] + g\eta = 0 \quad \text{at} \quad z = \eta(x,t)$$

Because the free-surface conditions are evaluated at an unknown moving boundary $z = \eta(x,t)$ and contain nonlinear terms, the system is analytically intractable. 
**Airy theory resolves this by assuming the wave amplitude $a$ is infinitesimally small compared to the wavelength $\lambda$ and water depth $h$ (i.e., small wave steepness $ak \ll 1$).** This permits the linearization of the boundary conditions and their evaluation at the mean water level ($z = 0$):

$$\frac{\partial \eta}{\partial t} = \frac{\partial \Phi}{\partial z} \quad \text{at} \quad z = 0$$ 

$$\frac{\partial \Phi}{\partial t} + g\eta = 0 \quad \text{at} \quad z = 0$$ 

### Derivation of the Dispersion Relation 
Assuming a monochromatic, progressive harmonic wave, the free surface elevation is defined as: 

$$\Large \eta(x,t) = a \cos(kx - \omega t)$$ 

Where $a$ is the wave amplitude, $k$ is the angular wavenumber ($2\pi/\lambda$), and $\omega$ is the angular frequency ($2\pi/T$). 

Applying separation of variables to the Laplace equation and enforcing the kinematic bed condition yields a velocity potential of the form:

$$\Large \Phi(x, z, t) = A \cosh(k(z+h)) \sin(kx - \omega t)$$ 

Substituting this potential and the surface elevation into the linearized dynamic and kinematic free-surface boundary conditions establishes an eigenvalue problem. For a non-trivial solution to exist, the relationship between the temporal parameter ($\omega$) and spatial parameter ($k$) must satisfy the fundamental linear wave dispersion equation:

$$\Large \omega^2 = g \cdot k \cdot \tanh(k \cdot h)$$ 

**This equation explicitly states that surface gravity waves are dispersive media; waves of different lengths travel at different speeds**. 
**The phase velocity ($c_p = \omega/k$) and group velocity ($c_g = \partial\omega/\partial k$) are wholly dictated by the depth-to-wavelength ratio**. Evaluating $k$ from a known $\omega$ and $h$ is precisely the core utility of the `wave-disp-equation.py` module.

### The Transcendental Nature of the Problem and Nondimensionalization 

In virtually all practical applications spanning oceanography and coastal engineering, the wave period $T$ (and thus $\omega$) is measured via wave buoys or synthesized by spectral models, and the bathymetric depth $h$ is acquired via hydrographic surveys. The engineer must calculate the wavenumber $k$ (and wavelength $\lambda$) to subsequently derive phase speed, group celerity, refraction indices, and shoaling coefficients. 

**However, the fundamental equation $\omega^2 = gk \tanh(kh)$ cannot be algebraically inverted to isolate $k$**. 
**The simultaneous presence of $k$ as a linear multiplier and as the argument within the hyperbolic tangent function renders the equation mathematically transcendental**. Such equations explicitly defy closed-form analytical solutions utilizing standard elementary functions, establishing a massive computational bottleneck in large-scale models. 

To computationally attack the problem, the `wave-disp-equation.py` repository first standardizes the equation into a universally applicable dimensionless form. By defining the deep-water wavenumber $k_0 = \omega^2/g$ (which corresponds to a deep-water wavelength $L_0 = gT^2/2\pi$), the dispersion relation transforms into:

$$\Large k_0 h = kh \cdot \tanh(kh)$$ 

To simplify the mathematical notation utilized throughout the module's algorithms, the independent dimensionless deep-water parameter is defined as $\alpha = k_0h$, and the dependent dimensionless wavenumber is defined as $\beta = kh$. 
The canonical equation thus becomes:

$$\Large \alpha = \beta \tanh \beta$$ 

**The entire algorithmic architecture of `wave-disp-equation.py` is singularly dedicated to mapping the known input $\alpha$ to the unknown output $\beta$ efficiently, accurately, and stably across all possible hydrodynamic depth regimes**.

### Asymptotic Limits: The Physical and Mathematical Boundaries

The design of the computational algorithms relies heavily on the asymptotic behavior of the hyperbolic tangent function ($\tanh \beta \rightarrow \beta$ as $\beta \rightarrow 0$, and $\tanh \beta \rightarrow 1$ as $\beta \rightarrow \infty$). **Because the transitional intermediate depths lack a closed-form analytical solution, explicit approximations and fractional Padé polynomials (like Carvalho's 2025 formulas) are structurally anchored to these two extreme physical boundaries to guarantee global mathematical stability**. 

* **Deep Water Limit ($\alpha \rightarrow \infty$):** Physically, this hydrodynamic regime occurs when the water depth is significantly greater than half the wavelength ($h/L_0 > 0.5$). In this domain, the sub-surface wave orbital motions decay exponentially with depth and dissipate entirely before interacting with the seabed. Mathematically, as the argument $\beta \rightarrow \infty$, the hyperbolic tangent geometrically flattens and approaches unity ($\tanh \beta \rightarrow 1$). 
    The canonical dispersion equation thus reduces to:

    $$\Large \beta \approx \alpha$$

    **This strict equality implies that $kh \approx k_0h$, meaning the local wavenumber $k$ converges entirely to the deep-water wavenumber $k_0$, and the local wavelength $L$ exactly equals $L_0$**. 
    By substituting $k = \omega^2/g$, the phase speed becomes:

    $$\Large c_p = \frac{\omega}{k} \approx \frac{\omega}{k_0} = \frac{\omega}{\omega^2 / g} = \frac{g}{\omega} = \frac{g}{2\pi}T$$

    **In the deep-water limit, the phase speed is directly proportional to the wave period $T$ and completely independent of the bathymetric depth $h$**. This frequency-dependent velocity represents pure dispersion and explains why long-period ocean swells outrun short-period wind chop across open ocean basins.

* **Shallow Water Limit ($\alpha \rightarrow 0$):** Physically, this regime occurs when waves propagate into very thin nearshore water layers ($h/L_0 < 0.05$). Here, the wave orbital motions become highly compressed and elliptical, dragging heavily against the solid seabed. Mathematically, by applying a Maclaurin series expansion for small arguments, $\tanh(\beta) \approx \beta - \beta^3/3 + \dots \approx \beta$. 
    Substituting this linear approximation into the canonical dispersion relation yields:

    $$\Large \alpha \approx \beta \cdot (\beta) = \beta^2$$

    Taking the square root mathematically isolates the dependent variable:

    $$\Large \beta \approx \sqrt{\alpha}$$

    **By unpacking the dimensional parameters ($k \cdot h \approx \sqrt{\frac{\omega^2}{g} \cdot h}$), we extract the shallow-water wavenumber. In this regime, phase speed becomes:**

    $$\Large c_p = \frac{\omega}{k} = \frac{\omega}{\frac{\omega}{\sqrt{gh}}} = \sqrt{gh}$$

    **Crucially, this phase speed is entirely dictated by depth ($h$) and is completely independent of wavelength or frequency, resulting in a strictly non-dispersive wave profile**. In shallow water, all waves (regardless of their period) travel at the exact same velocity, which is the fundamental mechanism driving wave shoaling and forcing the uniform parallel alignment of wave crests to the coastline via refraction.

---

## Exact Numerical Resolution: The Newton-Raphson Implementation 
**The module establishes its baseline of mathematical truth through the `kh_numeric(k0h, tol=1e-15, max_iter=100)` function, which utilizes the Newton-Raphson iterative root-finding algorithm to converge upon a numerically exact solution.** 

### Algorithmic Mechanics 
The Newton-Raphson method isolates the root of an objective function $f(\beta) = 0$. 
Using the non-dimensionalized variables, the residual function is established as:

$$f(\beta_n) = \alpha - \beta_n \tanh(\beta_n)$$ 

To apply the Newton-Raphson update, the analytical derivative (Jacobian) of the residual function with respect to the unknown $\beta$ must be derived. Applying the product rule and acknowledging that $\frac{d}{d\beta}\tanh(\beta) = \text{sech}^2(\beta)$ yields:

$$f'(\beta_n) = -\tanh(\beta_n) - \beta_n \text{sech}^2(\beta_n)$$ 

Within `wave-disp-equation.py`, this derivative is computationally implemented utilizing the fundamental identity $\text{sech}(\beta) = 1/\cosh(\beta)$:

```python
df = -math.tanh(kh) - kh / math.cosh(kh)**2
``` 

The iterative scheme advances by projecting the tangent line of the current guess to the zero intercept:

$$\beta_{n+1} = \beta_n - \frac{f(\beta_n)}{f'(\beta_n)} = \beta_n - \frac{\alpha - \beta_n \tanh(\beta_n)}{-\tanh(\beta_n) - \beta_n \text{sech}^2(\beta_n)}$$ 

The iteration loop rigorously evaluates the relative step size $|\Delta \beta / \beta_n|$. 
When this fractional change drops below the strict tolerance parameter `tol=1e-15`, the loop terminates, having achieved double-precision floating-point exactness. A safety break is enforced via `max_iter=100` to prevent infinite looping in the event of pathological input data. 

### Stability Topography and Carvalho's Seed: Evolutionary Optimization via GEP

**The standard residual function $f(\beta) = \alpha - \beta \tanh(\beta)$ inherently possesses an inflection point at approximately $\beta \approx (6/5)$ (where $\alpha \approx 1.003$ or $h/L_0 \approx 0.159$) .** At an inflection point, the second derivative $f''(\beta) = 0$. 
If an intermediate iteration step lands precisely near this topological feature, the tangent slope becomes momentarily flat, which violently throws the subsequent estimate $\beta_{n+1}$ far outside the physical domain, inducing numerical instability or total algorithmic divergence. 
Historically, researchers like Goda advocated mutating the residual function to $f(\beta) = \beta - \alpha \coth(\beta)$ or $f(\beta) = \alpha/\beta - \tanh(\beta)$ to globally eliminate the inflection point. **However, the `wave-disp-equation.py` repository circumvents this geometrical hazard entirely by supplying an exceptionally sophisticated initial guess $\beta_0$, documented within the source code as the "Carvalho, 2006 style" seed.** 

**This specific formula was developed by R. Carvalho in 2006 using Gene Expression Programming (GEP), a high-performance evolutionary algorithm.** GEP differs from traditional genetic algorithms by using linear chromosomes (the genotype) that are expressed as branched expression trees (the phenotype), allowing the system to discover complex nonlinear algebraic relationships that human derivation might overlook. **Carvalho’s GEP-derived seed represents a breakthrough in performance engineering, achieving superior precision (0.27% maximum error) compared to traditional explicit bridges.**

$$\Large \beta_0 = \frac{\alpha}{\tanh\left( \left(\frac{6}{5}\right)^\alpha \sqrt{\alpha} \right)}$$ 

This specific seed algorithm is a masterclass in asymptotic matching. 
**In the shallow limit ($\alpha \rightarrow 0$), the fractional exponential $(6/5)^0 \rightarrow 1$, leaving $\tanh(\sqrt{\alpha})$.** For small arguments, $\tanh(x) \approx x$, meaning the denominator resolves to $\sqrt{\alpha}$. 
**The entire expression flawlessly reduces to $\beta_0 = \alpha / \sqrt{\alpha} = \sqrt{\alpha}$, precisely matching shallow water wave theory.** Conversely, in the deep-water limit ($\alpha \rightarrow \infty$), the argument of the hyperbolic tangent approaches infinity, forcing the denominator to exactly $1$, causing the expression to reduce to $\beta_0 = \alpha$, perfectly tracking the deep-water boundary. **The empirical constant $\frac{6}{5}$ (1.2) derived by GEP smoothly morphs the transition zone, ensuring the initial guess is so geometrically close to the true root that the Newton-Raphson method skips the inflection point entirely, typically converging to $10^{-15}$ precision within 3 iterations.**

---

## Group I1: Classical and Contemporary Explicit Approximations (Full Domain) 
**While the `kh_numeric` function provides flawless mathematical precision, the repetitive evaluation of expensive transcendental functions ($\tanh$, $\cosh$) inside a `while` loop is highly deleterious to the performance of large-scale, phase-resolving coastal models and global spectral grids, which evaluate the dispersion relation millions of times per temporal step.** To alleviate this computational burden, scientists have engineered explicit, non-iterative approximations. 
The repository categorizes its explicit formulas into specific groups. Group I1 details analytical solutions valid across the entire physical domain ($0 < h/L_0 < \infty$). 
These solutions vary in algebraic complexity and maximum relative error $\tilde{\epsilon}$. The formulas encoded within the repository framework are mathematically deconstructed below. 

**1. The Eckart Solution (1951)** 
 Carl Eckart formulated the foundational explicit approximation by recognizing that multiplying the deep water linear limit by a hyperbolic tangent containing the shallow water limit generates a globally continuous curve:

 $$k_a h = \alpha (\coth \alpha)^{1/2}$$ 

 Equivalently expressed in terms of wavelength: 

 $$L_a = L_0 (\tanh \alpha)^{1/2}$$ 

 This pioneering equation guarantees boundary compliance. However, its rigid geometrical structure lacks tunable degrees of freedom in intermediate depths. 
 Consequently, the approximation bows away from the exact solution, producing a maximum positive error of $+5.24\%$ centered precisely at $h/L_0 = 0.111$. While computationally lightweight, a 5% error fundamentally distorts wave celerity in the shoaling zone, making it unsuitable for modern quantitative analysis. 

**2. The Iwagaki Solution** 
 To constrain the intermediate depth error, Iwagaki injected a corrective linear growth term into the hyperbolic argument:

 $$k_a h = \alpha \coth\left\{ \alpha^{1/2} \left( 1 + \frac{\alpha^{1/2}}{2\pi} \right) \right\}$$ 

 The added internal degrees of freedom force the curve closer to the exact solution, producing an oscillatory error bound. 
 The error spans from a negative maximum of $-3.05\%$ at $h/L_0 = 0.287$ to a positive maximum of $+3.14\%$ at $h/L_0 = 0.023$. 

**3. Yamaguchi's 4th Solution (`Yamaguchi(2007)_4`)** Yamaguchi and Nonaka identified that eliminating transcendental functions entirely could vastly accelerate processing times on legacy hardware. This explicit solution operates solely on algebraic fractions and roots:

$$k_a h = \alpha(1 + \alpha^{-2})^{1/4}$$ 

Despite bypassing hyperbolic geometry entirely, the formula manages to restrict the maximum relative error to 3.177% (occurring at $\alpha = 0.4268$).

**4. Fenton and MacKee (FM) Solution (1990)** 
 Fenton and MacKee revolutionized explicit modeling by generalizing Eckart's static $1/2$ exponent into a tunable parameter $1/m$:

 $$k_a h = \alpha (\coth \alpha^{m/2})^{1/m}$$ 

 Through minimax error optimization, they deduced that setting $m = 1.5$ dramatically tightened the curve across the transitional depths. 
 The FM solution confines the relative error strictly between $-1.39\%$ (at $h/L_0 = 0.321$) and $+1.66\%$ (at $h/L_0 = 0.054$). 

**5. Yamaguchi and Nonaka Modifications (YN1 and YN2)** 
 Yamaguchi and Nonaka (2007) observed that the FM solution produced an asymmetrical error domain. 
 * **YN1 Solution:** By perturbing the exponent to $m = 1.485$, Yamaguchi and Nonaka perfectly balanced the positive and negative maximum errors, resulting in an error band strictly bounded by $\pm 1.55\%$. 
 * **YN2 Solution:** To achieve sub-1% accuracy, they proposed nesting the YN1 explicit estimate directly back into a rearranged form of the analytical dispersion relation. By calculating a temporary variable $\beta_F = \alpha(\coth \alpha^{m/2})^{1/m}$ with a customized parameter $m = 1.378$, the final wavenumber is extracted via:
 $$k_a h = \alpha \coth \beta_F$$ 
 This nested architecture simulates a single quasi-Newton corrective step, halving the maximum error to a symmetrical $\pm 0.73\%$. 

**6. The Guo Solution (2002)** 
 Jie Guo engineered a fundamentally different geometric approach, abandoning the hyperbolic cotangent in favor of a standard exponential decay function:

 $$k_a h = \frac{\alpha}{\{ 1 - \exp(-\alpha^{m/2}) \}^{1/m}}$$ 

 Utilizing a highly optimized parameter $m = 2.4901$, Guo's formula produces an exact $\pm 0.75\%$ error bound. **The mathematical superiority of this formula lies in its computational velocity; the standard exponential function exp() is evaluated significantly faster by ALUs (Arithmetic Logic Units) than the hyperbolic functions tanh() or coth(), resulting in massive speedups in vectorized array processing.** 

**7. Carvalho's 18th, 10th, and 9th Solutions** The module addresses a spectrum of formulas derived by Carvalho utilizing a variety of mathematical structures, ranging from simple dual-nested hyperbolic functions to exponential asymptotic geometries. 
* **Carvalho(2025)_18:** $k_a h = \alpha \coth(\sinh \alpha^{1/2})$. By cleanly nesting a single hyperbolic sine inside a hyperbolic cotangent, this configuration yields a strictly negative error profile, bounding the maximum relative error to 1.129% (occurring at $\alpha = 1.4912$). 
* **Carvalho(2025)_10:** $k_a h = \alpha \coth((6/5)^\alpha \cdot \alpha^{1/2})$. This formula provides excellent precision, with the maximum relative error constrained to 0.271% (occurring at $\alpha = 0.3941$). As previously discussed, this specific $1.2^\alpha$ (or $(6/5)^\alpha$) geometry leverages an exponential fractional base rather than deep nesting. It is precisely what powers the highly stable initial guess for the `kh_numeric` Newton solver. 
* **Carvalho(2025)_9:**
$$k_a h = \frac{\alpha}{(\tanh \alpha)^{1/4} [\tanh \{ (\sinh \alpha)^{1/2} \}]^{1/2}}$$
Unlike the simpler configurations above, this specific formula aggressively triply nests hyperbolic functions alongside radical roots. This deep topological nesting suppresses the maximum relative error down to 0.204% (occurring at $\alpha = 2.6569$). However, the computational latency induced by evaluating three separate transcendental instructions renders this particular formula inefficient compared to fractional Padé alternatives.

**8. Carvalho's Evolutionary GEP Solutions (2025) and Yamaguchi Modifications** The module includes a vast array of explicit formulas derived via Gene Expression Programming (GEP) by Carvalho, alongside modifications detailed by Yamaguchi & Nonaka. By strictly analyzing the `wave-disp-equation.py` Python array outputs, the most notable mathematically distinct explicit formulas are deconstructed below:

**The Deeply Nested Fraction (`Carvalho(2025)_4`):**
  $$\Large k_a h = \frac{\alpha}{\tanh\left( \frac{\alpha}{\tanh\left( \frac{\alpha}{\tanh\left( \frac{\alpha}{\sinh(\tanh(\sqrt{\alpha}))} \right)} \right)} \right)}$$
  **By deeply nesting hyperbolic functions in a recursive continuous fraction architecture, the maximum relative error is drastically suppressed down to an exceptional $0.050\%$ (occurring at $\alpha = 0.3463$)**. However, the extreme computational latency induced by evaluating six separate sequential transcendental instructions per grid point renders this specific formula mathematically inefficient when compared to modern fractional Padé alternatives.

**The Optimized Exponential (`Carvalho(2025)_5`):**
  $$\Large k_a h = \frac{\alpha}{\tanh\left(1.199315^{\alpha^{1.047086}} \cdot \alpha^{0.499947}\right)}$$
  **This GEP derived formula highly tuned by empirical coefficients provides exceptional precision, sharply bounding the maximum relative error to just $0.076\%$ (occurring at $\alpha = 1.5603$)**. It represents an extremely optimized evolution of the $(6/5)^\alpha$ seed geometry.

**The Triple Nested Root (`Carvalho(2025)_9`):** $$\Large k_a h = \frac{\alpha}{(\tanh \alpha)^{1/4} \sqrt{\tanh\left(\sqrt{\sinh \alpha}\right)}}$$
  By triply nesting hyperbolic functions and radical roots, the curve is forced to tightly hug the exact transitional depths, restricting the maximum relative error to $0.204\%$ (occurring at $\alpha = 2.6569$). 

**The Asymptotic Seed (`Carvalho(2025)_10`):** $$\Large k_a h = \frac{\alpha}{\tanh\left( \left(\frac{6}{5}\right)^\alpha \sqrt{\alpha} \right)}$$
  This is the masterful GEP asymptotic bridge that powers the highly stable initial guess for the `kh_numeric` Newton solver. **It seamlessly connects shallow and deep asymptotes, limiting the maximum relative error across all depths to $0.271\%$ (occurring at $\alpha = 0.3941$)**.

* **The Hyperbolic Sine Root (`Carvalho(2025)_18`):** $$\Large k_a h = \frac{\alpha}{\tanh(\sinh \sqrt{\alpha})}$$
  This configuration eliminates positive error overshoot, yielding a strictly negative error profile bounded by a maximum relative error of $1.129\%$ (occurring at $\alpha = 1.4912$).

* **The Algebraic Solution (`Yamaguchi(2007)_4`):** $$\Large k_a h = \alpha(1 + \alpha^{-2})^{1/4}$$
  This lightweight equation bypasses hyperbolic geometry entirely. **It operates strictly on algebraic fractions and roots, achieving a maximum relative error of $3.177\%$ (occurring at $\alpha = 0.4268$)**. It is ideal for legacy systems where hardware ALUs execute simple roots faster than transcendentals.

---

### Formulas Maximum Errors Ranking Table

| Model Name | Equation Format | MaxErr |
| :--- | :--- | :--- |
| Carvalho(2025)_4 | Deep continuous fraction | $0.050\%$ |
| Carvalho(2025)_5 | $\alpha / \tanh(1.199^{\alpha^{1.047}} \alpha^{0.499})$ | $0.076\%$ |
| Carvalho(2025)_9 | $\alpha / \{(\tanh \alpha)^{1/4} \sqrt{\tanh(\sqrt{\sinh \alpha})}\}$ | $0.204\%$ |
| Carvalho(2025)_10 (Seed)| $\alpha / \tanh((6/5)^\alpha \alpha^{1/2})$ | $0.271\%$ |
| Yamaguchi(2007)_2 (YN2) | Nested FM substitution $(m=1.378)$ | $0.732\%$ |
| Guo(2002) | $\alpha / \{ 1 - \exp(-\alpha^{m/2}) \}^{1/m} \quad (m=2.4901)$ | $0.757\%$ |
| Carvalho(2025)_18 | $\alpha / \tanh(\sinh \alpha^{1/2})$ | $1.129\%$ |
| Yamaguchi(2007)_1 (YN1) | $\alpha(\coth \alpha^{m/2})^{1/m} \quad (m=1.485)$ | $1.543\%$ |
| Fenton&McKee(1990)_2 | $\alpha / (\tanh \alpha^{3/4})^{2/3}$ | $1.631\%$ |
| Iwagaki(2007) | $\alpha / \tanh\{\alpha^{1/2}(1 + \alpha^{1/2}/2\pi)\}$ | $3.147\%$ |
| Yamaguchi(2007)_4 | $\alpha(1 + \alpha^{-2})^{1/4}$ | $3.177\%$ |
| Eckart(1951) | $\alpha / (\tanh \alpha)^{1/2}$ | $4.980\%$ |
---

## Group I2: High-Order Padé Approximants and Rational Functions 
When resolving complex nearshore wave-current interactions or phase-resolving Boussinesq hydrodynamics, standard explicit bounds of $\sim 0.5\%$ introduce unacceptable phase errors over long temporal simulations. **To achieve arbitrary sub-0.1% precision without the algorithmic branching of Newton iterations, mathematicians deploy Padé approximants.** 
A Padé approximant represents a given function as the ratio of two polynomials. 
It generally outperforms Taylor series truncations because rational functions can seamlessly simulate poles, asymptotes, and nonlinear saturation geometries, making them the perfect analytical tool for modeling the hyperbolic tangent transition zone. 

### The Hunt Solutions (1979) 
J.N. Hunt published the definitive early application of Padé approximants for wave dispersion, expanding the squared dimensionless wavenumber $(k_ah)^2$ as a rational polynomial of $\alpha$. 

* **Hunt's 5th Order Approximant (Hunt1):** 
 $$(k_a h)^2 = \alpha^2 + \frac{\alpha}{1 + 0.6522\alpha + 0.4622\alpha^2 + 0.0864\alpha^4 + 0.0675\alpha^5}$$ 
 Notice the deliberate mathematical omission of the cubic parameter ($\alpha^3$); this forces the polynomial to warp across the transitional intermediate depths smoothly, achieving an impressive error boundary between $-0.070\%$ and $+0.078\%$. 

* **Hunt's 9th Order Approximant (Hunt2):** 
 To push precision to the extreme, Hunt formulated a 9th order denominator polynomial: 
 $$(k_a h)^2 = \alpha^2 + \frac{\alpha}{1 + \sum_{n=1}^9 d_n \alpha^n}$$ 
 Where the precisely tuned coefficients are: $d_1 = 0.66667$, $d_2 = 0.35550$, $d_3 = 0.16084$, $d_4 = 0.06320$, $d_5 = 0.02174$, $d_6 = 0.00654$, $d_7 = 0.00170$, $d_8 = 0.00039$, and $d_9 = 0.00010$. This expansion compresses the relative error to an exceptionally tight margin of $-0.0082\%$ to $+0.0054\%$. Despite the high degree, evaluating a 9th order polynomial using Horner's method requires only basic multiplication and addition, entirely eliminating CPU-intensive transcendental logic. 

### Carvalho's 2025 Fractional Padé Approximations 
The repository `wave-disp-equation.py` deploys state-of-the-art fractional Padé approximants within the `pade2025(k0h, formula)` function, mapping formulas 1 through 13. **Padé approximants are rational functions that approximate a given function by matching its Taylor series expansion up to a specified order**. **They use ratios of polynomials instead of solely polynomial expansions, offering vastly superior accuracy, particularly around singularities, asymptotes, and for functions exhibiting nonlinear saturation geometries like the hyperbolic tangent transition zone**.

**Derived utilizing gene expression programming and analytical minimax optimization, Carvalho's novel mathematical insight was constructing the Padé polynomials utilizing half-integer powers ($\alpha^{0.5}, \alpha^{1.5}$) instead of standard integer powers**. 
This fractional architecture natively embeds the $\sqrt{\alpha}$ shallow-water asymptote directly into the fundamental basis vectors of the polynomial. 
This completely nullifies the necessity of squaring the left-hand side and computing an expensive square root at the end of the calculation, as required by the Hunt formulation. 
The function defines the following specific high-accuracy fractional polynomials:

* **Formula 1 (`Pade(2025)_1`):** $$k_a h = \frac{n_1 \alpha^{0.5} + n_2 \alpha^{1.5} + n_3 \alpha^{2.5}}{1 + d_1 \alpha^{1.0} + d_2 \alpha^{2.0}}$$ 
  
  Where the numerator coefficients are:
  * $n_1 = 1.00649052194019$
  * $n_2 = 0.423646282789217$
  * $n_3 = 0.175406661440005$
  
  And the denominator coefficients are:
  * $d_1 = 0.306955955676234$
  * $d_2 = 0.0328975279727171$

  As $\alpha \rightarrow 0$, the lowest-order terms dominate, collapsing the fraction to $1.006 \alpha^{0.5} / 1$, correctly modeling the shallow limit. As $\alpha \rightarrow \infty$, the highest-order terms dominate, balancing out to a linear $\alpha^{1.0}$ trajectory to flawlessly mirror the deep-water limit. **This fundamental configuration provides an excellent baseline, yielding a maximum relative error of just $0.6485218\%$ (occurring at $\alpha = 0.0001$)**.

* **Formula 2 (`Pade(2025)_2`):** $$k_a h = \frac{n_1 \alpha^{0.5} + n_2 \alpha^{1.5} + n_3 \alpha^{2.5} + n_4 \alpha^{3.5}}{1 + d_1 \alpha^{1.0} + d_2 \alpha^{2.0} + d_3 \alpha^{3.0}}$$
  
  Where the numerator coefficients are:
  * $n_1 = \phantom{-}0.998980252114366$
  * $n_2 = \phantom{-}0.0240176797055886$
  * $n_3 = \phantom{-}0.102524886754552$
  * $n_4 = \phantom{-}0.0317327085938995$
  
  And the denominator coefficients are:
  * $d_1 = -0.150350405960952$
  * $d_2 = \phantom{-}0.112157962910113$
  * $d_3 = \phantom{-}0.00294483072586115$
  
  The deliberate introduction of a negative linear coefficient ($-0.15035...\alpha^{1.0}$) in the denominator acts as a topological inflection constraint, forcing the polynomial curve to hug the transcendental $\tanh$ geometry with supreme accuracy throughout intermediate depth transitions. **This structural refinement dramatically constrains the maximum relative error down to a highly accurate $0.1018976\%$ (occurring at $\alpha = 0.0001$)**.

* **Formula 3 (`Pade(2025)_3`):** $$k_a h = \frac{n_1 \alpha^{0.5} + n_2 \alpha^{1.5} + n_3 \alpha^{2.5} + n_4 \alpha^{3.5} + n_5 \alpha^{4.5}}{1 + d_1 \alpha^{1.0} + d_2 \alpha^{2.0} + d_3 \alpha^{3.0} + d_4 \alpha^{4.0}}$$
  
  Where the numerator coefficients are:
  * $n_1 = 1.00006668638419$
  * $n_2 = 0.322645945302282$
  * $n_3 = 0.0860384450810725$
  * $n_4 = 0.051143347041175$
  * $n_5 = 0.0153420957423937$
  
  And the denominator coefficients are:
  * $d_1 = 0.157166943736625$
  * $d_2 = 0.0245168267924732$
  * $d_3 = 0.0462567432956417$
  * $d_4 = 0.00175392506101448$
  
  This expansion achieves phenomenal numerical stability. **By expanding the numerator to the 4.5-degree and the denominator to the 4.0-degree, this formulation pushes the precision even further, limiting the maximum relative error to an exceptional $0.0066566\%$ (occurring at $\alpha = 0.0001$)**. Because the mathematical operations rely exclusively on basic arithmetic and a single primitive square root, the `pade2025` function can be mapped directly onto highly parallelized tensor structures or GPU arrays without encountering the thread divergence issues caused by iterative conditional checking. Note that the programmatic array restricts `k0h` inputs to between $0$ and $2\pi$ to prevent unbounded polynomial resonance.

---

## Restricted Domain Approximations: Groups II and III 
For specialized oceanographic models that process hydrodynamics exclusively in the far open ocean (deep water) or exclusively in the nearshore surf zone (shallow water), computing a global explicit formula carries unnecessary mathematical overhead. 
Groups II and III encompass solutions defined strictly by localized Taylor and asymptotic series centered around $\alpha \rightarrow 0$ or $\alpha \rightarrow \infty$. 

### Shallow Water Formulations ( $h/L_0 \rightarrow 0$ ) 
These formulations leverage the Maclaurin series expansion of the fundamental transcendental equation, tailored for nearshore environments. 
* **Nielsen's 1st Solution (Niell):** By taking the first-order correction, $k_a h = \alpha^{1/2} \{1 + (5/8\pi)\alpha\}$. This limits the error to $\pm 0.74\%$, but mathematically diverges if applied past the surf zone bound of $h/L_0 \le 0.192$. 
* **Nielsen's 2nd Solution (Niel2):** Adding second-order terms extends the domain seaward: $k_a h = \alpha^{1/2} \{1 + (1/6)\alpha + (11/360)\alpha^2\}$. Error remains $\pm 0.44\%$ out to $h/L_0 \le 0.401$. 
* **Venezian's 1st Solution (Venel):** Applies a simple rational denominator shift: $k_a h = \alpha^{1/2} / (1 - \alpha/6)$. Maintains a tight $\pm 0.048\%$ error constraint up to $h/L_0 \le 0.165$. 
* **Wu and Thornton 1st Solution (WT1):** $k_a h = \alpha^{1/2} \{1 + \frac{\alpha}{6}(1 + \frac{\alpha}{5})\}$. Maintains precision out to $h/L_0 \le 0.219$. 
* **Olson (1973) and You (2008):** Olson aggressively pushed the polynomial to the 7th order: $k_a h = \alpha^{1/2} \{1 - \frac{1}{3}\alpha + \frac{1}{45}\alpha^2 + \frac{1}{189}\alpha^3 + 0.000776\alpha^4 - 0.000044\alpha^5 - 0.000071\alpha^6 - 0.000022\alpha^7 \}^{-1/2}$. This drives the error down to a microscopic $\pm 0.00003\%$ within the shallow zone $h/L_0 \le 0.186$. 
* **Venezian's 2nd Solution (Vene2):** Employs a specific, highly optimized nearshore Padé rational fraction:
 $$k_a h = \alpha^{1/2} \frac{1 - 0.428868\alpha + 0.093928\alpha^2 - 0.002694\alpha^3}{1 - 0.595534\alpha + 0.162628\alpha^2 - 0.014975\alpha^3}$$ 
 This yields extreme sub-surface accuracy ($-0.0002\%$ to $+0.000006\%$), but suffers catastrophic divergence if pushed past $h/L_0 > 0.159$. 

### Deep Water Formulations ( $h/L_0 \rightarrow \infty$ ) 
As depth increases towards infinity, $\tanh \beta \rightarrow 1$ following an asymptotic exponential decay trajectory proportional to $1 - 2e^{-2\beta} + \dots$. 
* **Nielsen's 3rd Solution (Niel3):** $k_a h = \alpha \{1 + 2\exp(-2\alpha)\}$. This effectively models the deep-water limit, maintaining an error of $-0.55\%$ down to $h/L_0 \ge 0.300$, but diverges into immense errors in shallow water. 
* **Wu and Thornton 2nd Solution (WT2):** Integrates nested exponentials to force the curve to hug the transitional depths slightly closer to shore: $k_a h = \alpha \{1 + 2t(1+t)\}$, utilizing the dynamic variable $t = \exp\{-2\alpha(1 + 1.26e^{-1.84\alpha})\}$. This restricts error to $\pm 0.025\%$ for conditions where $h/L_0 \ge 0.195$. 

**Algorithmic Warning:** Modeling systems frequently attempt to merge a shallow formulation (like WT1) with a deep formulation (like WT2) via a piecewise conditional logic switch at a specific depth boundary (e.g., $h/L_0 = 0.2$). This piecewise architecture inevitably introduces a mathematical discontinuity (a discrete "jump" in phase speed mapping) at the transition point, causing severe non-physical wave scattering artifacts in refraction-diffraction grids. 
Therefore, seamless full-domain approximations (Group I) or high-order Padé functions remain mathematically superior for generalized computing environments. 

---

## Contemporary Optimization: The Seed-and-Step Paradigm 
Recent theoretical advancements documented in the repository eschew the derivation of progressively complex polynomial expressions in favor of an elegant hybrid approach known as the "seed-and-step" method. 

### The Beji Formulation (2013) 
Serdar Beji established a highly accurate modification of Eckart's foundational equation. 
By proving that Eckart's basic square-root geometry fundamentally lacked the degrees of freedom necessary to mimic intermediate depth transitions, Beji injected a sophisticated empirical exponential correction function:

$$k_a h = \alpha \exp(-\beta_0 + \beta_1\alpha + \beta_2\alpha^2)$$ 

Through rigorous least-squares curve fitting across the operational domain, Beji determined the optimum constants to constrain the absolute maximum relative error to $0.044\%$ across all physical depths. 

### Vatankhah and Aghashariatmadari Refinement (2013) 
Almost immediately following Beji's publication, researchers Vatankhah and Aghashariatmadari applied advanced geometric curve-fitting algorithms directly to Beji’s core exponential expression. **By optimizing the exponent constants, they drove the maximum absolute error down to $0.019\%$ without imposing any additional computational latency upon the host hardware**. To accommodate varying precision requirements, they proposed two distinct single-step explicit formulas:

* **Single-Step Explicit Formula #1:**
 This formulation directly refines Beji's exponent terms to smoothly approximate the intermediate depth transition.

 $$\Large k_a h \approx \frac{\alpha + \alpha^2 \exp(-1.835 - 1.225\alpha^{1.35})}{\sqrt{\tanh(\alpha)}}$$

 **This streamlined geometric equation strictly bounds the maximum relative error to $0.0189\%$ (occurring at $\alpha = 0.0705$)**.

* **Single-Step Explicit Formula #2:**
 For scenarios demanding even tighter topological constraints, they split the approximation into a two-part expansion, introducing a sophisticated, high-order exponential corrective growth term:

 $$\Large k_a h \approx \frac{\alpha + \alpha^2 \exp(-(3.2 + \alpha^{1.65}))}{\sqrt{\tanh(\alpha)}} + \alpha \left[1 - \exp(-\alpha^{0.132})\right]^{(5.0532 + 2.1584\alpha^{1.505})}$$

 **This advanced split-form configuration aggressively suppresses the maximum relative error down to an exceptional $0.00176\%$ (occurring at $\alpha = 0.9515$), yielding near-exact precision while remaining fully explicit and avoiding iterative branching**.

### Simarro and Orfila's Single-Step Method (2013) 
The pinnacle of contemporary optimization was realized by Simarro and Orfila. They postulated that pursuing mathematically convoluted explicit analytical formulas (such as Carvalho's triple-nested hyperbolic roots or Hunt's massive 9th-order fractions) to achieve near-exact accuracy was fundamentally inefficient. 
**Instead, they utilized Beji's robust $0.044\%$ accurate explicit formulation not as a final solution, but exclusively as the $\beta_0$ initial seed for exactly one fixed Newton-Raphson iteration.** 
Because the Beji approximation places the algorithmic state extraordinarily close to the true root, taking exactly one Newton step yields: 

$$\beta_{\text{final}} = \beta_{\text{Beji}} - \frac{\alpha - \beta_{\text{Beji}} \tanh(\beta_{\text{Beji}})}{-\tanh(\beta_{\text{Beji}}) - \beta_{\text{Beji}} \text{sech}^2(\beta_{\text{Beji}})}$$ 

**This singular mathematical maneuver plunges the relative error down to an astonishing $0.0000082\%$.** Further refinement utilizing the Vatankhah modified seed pushed this precision to $0.00000028\%$. This hybrid "seed-and-step" topology represents a computational Holy Grail: it provides double-precision analytical exactness while entirely circumventing the need for while loops, conditional convergence thresholds, and algorithmic branching logic. This makes the logic exceptionally efficient when executing highly vectorized wave forecasting models operating on massively parallel GPU (CUDA/OpenCL) architectures where warp divergence imposes severe computational penalties. 

---

## Physical Extensions to the Dispersion Framework 
The standard Airy dispersion relation simulated in `wave-disp-equation.py` ($\omega^2 = gk \tanh(kh)$) defines the fundamental kinematic matrix. 
However, highly complex ocean waves phenomena force crucial theoretical extensions to this core relationship, which modelers must integrate. 

### 1. Wave-Current Interaction and Doppler Shifting 
When wave trains propagate across an underlying mean fluid current characterized by a depth-averaged velocity vector $\mathbf{V}$, the dispersion relation must undergo a Galilean transformation, resulting in a hydrodynamic Doppler shift. The absolute frequency observed from a stationary geographic reference frame ($\omega$) separates from the intrinsic relative frequency ($\sigma$) observed in a frame moving alongside the fluid. 
The generalized dispersion relation becomes:

$$(\omega - \mathbf{k} \cdot \mathbf{V})^2 = g k \tanh(kh) \quad \implies \quad \sigma^2 = g k \tanh(kh)$$ 

Where $\mathbf{k}$ is the spatial wavenumber vector. If the fluid current flows directly opposite to the wave propagation vector ($\mathbf{k} \cdot \mathbf{V} < 0$), the waves compress, steepen, and accumulate massive kinetic energy. 
If the adverse current velocity is sufficiently strong, total wave blocking occurs, physically preventing forward energy propagation. Models resolving wave-current interaction rely continuously on the rapid explicit algorithms documented above to resolve the intrinsic spatial parameters. 

### 2. Capillary Waves and Surface Tension Effects 
At microscopic spatial domains (specifically wavelengths $\lambda < 0.07$ meters), the hydrostatic restoring force of gravity is matched and eventually superseded by the cohesive fluid surface tension ($\gamma$ or $\sigma_t$). To accommodate this force, the gravity-capillary dispersion relation introduces a cubic wavenumber parameter:

$$\omega^2 = \left(gk + \frac{\gamma}{\rho} k^3\right) \tanh(kh)$$ 

Where $\rho$ defines the fluid density. Capillary waves exhibit inverted frequency dispersion physics: unlike gravity waves where longer waves traverse faster, capillary waves propagate faster as their wavelength becomes shorter. The presence of the $k^3$ parameter invalidates standard polynomial approximations, necessitating specific Padé expansions or dedicated root-finding mechanisms for microwave remote sensing and wind-wave generation models. 

### 3. Amplitude Dispersion and Stokes Non-Linearity 
The Airy wave framework is strictly a first-order linearized analytical solution, fundamentally assuming that wave steepness ($ak$) is infinitesimally small. Under real oceanic conditions, as waves steepen (e.g., hurricane-forced chaotic sea states), the nonlinear convective acceleration terms within the Navier-Stokes equations physically alter the wave's phase celerity. Applying Stokes' perturbation expansion theory to the third order in deep water modifies the governing relation to:

$\omega^2 = gk [1 + (ka)^2]$

This equation establishes the phenomenon of "amplitude dispersion." It mathematically proves that large-amplitude crests propagate tangibly faster than small-amplitude waves sharing the identical frequency. **This $(ka)^2$ phase dependency fundamentally governs nonlinear wave-wave interaction physics, directly precipitating the Benjamin-Feir modulational instability, which drives the spontaneous generation of extreme rogue (freak) waves in the open ocean.** 

### 4. Interfacial and Internal Wave Dispersion 
When evaluating stratified, multi-layered fluid bodies (such as a buoyant freshwater estuarine plume sliding over dense saline oceanic water), internal waves propagate silently along the subsurface pycnocline. 
For two semi-infinite homogeneous layers with lower density $\rho$ and upper density $\rho'$, the deep-water dispersion relation transitions into an internal density gradient expression:

$$\omega^2 = gk \left( \frac{\rho - \rho'}{\rho + \rho'} \right)$$ 

Because the physical density variation across the interface is extremely small (fractional percentages), the effective restoring force (reduced gravity) is weak. Consequently, internal waves propagate at a fraction of the speed of surface waves but can attain colossal vertical amplitudes spanning hundreds of meters. 

---

## Synthesis of Approximations and Model Accuracy

The extraction of the wavenumber from the dispersion equation is not merely an exercise in theoretical math; the value of $k$ completely dictates the calculation of second-order momentum flux, the radiation stress tensor, wave energy density ($E = \frac{1}{2}\rho g a^2$), and the shoaling and refraction coefficients that define all nearshore morphological evolution models.

To understand the immense physical weight of this mathematical inversion, one must recognize that **every kinematic and dynamic property of a surface gravity wave is fundamentally downstream of the wavenumber $k$**. Once $k$ is extracted, it immediately defines the phase velocity ($c_p$), which governs wave refraction over varying bathymetry:

$$\Large c_p = \frac{\omega}{k} = \frac{g}{\omega} \tanh(kh)$$

Furthermore, $k$ dictates the group celerity ($c_g$), which represents the actual velocity at which wave energy propagates across the ocean surface:

$$\Large c_g = \frac{\partial \omega}{\partial k} = \frac{1}{2} c_p \left[ 1 + \frac{2kh}{\sinh(2kh)} \right]$$

**Any fractional error introduced during the numerical estimation of $k$ will exponentially propagate through these secondary derivatives, fundamentally corrupting wave height predictions, energy flux balances, and sediment transport gradients in coastal models**. 

The extensive algorithmic landscape detailed within `wave-disp-equation.py`—from Eckart's elementary 1951 asymptotic bridge to the highly refined 2025 half-integer fractional Padé approximants by Carvalho—illustrates an evolutionary continuum of mathematical performance engineering. This evolution mirrors the historical progression of computational hardware: from the slide rules of the 1950s requiring simple algebraic bridges, to the single-core microprocessors of the 1990s optimizing minimax polynomials, and ultimately to the massively parallel Tensor and GPU architectures of the 2020s demanding strictly divergence-free, non-iterative logic.

By benchmarking these algorithms against a strict Newton-Raphson baseline, several distinct operational paradigms emerge for the modern oceanographer:

* **For pedagogical clarity or rapid generalized analysis, the purely algebraic Guo or Carv14 explicit formulations provide excellent computational velocity with minimal coding complexity**.
 * **The Guo (2002) formula** relies entirely on standard exponential decay, evaluating millions of arrays in mere fractions of a second (approx. 0.73 seconds per million operations) while globally constraining its error to a remarkably symmetrical $\pm 0.75\%$.
* **YN4 Solution** mathematically bypasses transcendentals (like `tanh` or `exp`) entirely, operating strictly on algebraic fractions and roots to achieve a maximum relative error of 3.177%. **These are ideal for legacy systems, lightweight embedded wave buoys, or rapid prototyping where hardware ALUs execute simple roots faster than complex hyperbolic geometry**.

* **For heavily constrained processing environments where iterative divergence poses a catastrophic failure risk, the Yamaguchi & Nonaka modifications or the robust Hunt Padé arrays guarantee strict topological bounds without requiring conditional looping logic**.
 * In massively parallel GPU computing (CUDA/OpenCL), traditional iterative `while` loops cause severe "warp divergence"—where some threads finish early while others stall, bottlenecking the entire array. 
 * **Hunt's 9th Order Approximant (Hunt2)** circumvents this by utilizing a massive 9th-degree polynomial ratio. Executed via Horner's method, it compresses maximum error to an exceptionally tight margin of just $0.00815\%$, utilizing only primitive multiplication and addition.
 * Similarly, the **Carvalho (2025) Fractional Padé Approximants** completely rewrite the underlying polynomial basis vectors using half-integer powers (e.g., $\alpha^{0.5}, \alpha^{1.5}$). **Formulas like Pade(2025)_8 and Pade(2025)_11 push the maximum relative error to absolute $0.0000000\%$ (matching double-precision numeric exactness) without a single conditional check, offering a true 1:1 replacement for iterative solvers on tensor cores**.

* **For modern, state-of-the-art predictive ocean waves modeling requiring supreme accuracy, the hybridization of an empirically exact seed (such as the Beji equation or the Carvalho initial guess) inextricably linked to a deterministic, single-step Newton-Raphson update yields sub-machine-epsilon accuracy ($< 10^{-7}\%$), entirely eradicating wavenumber-induced phase lag over multi-decadal simulation timelines**.
 * **This paradigm perfectly leverages Carvalho's 2006 Gene Expression Programming (GEP) derived seed, which elegantly maps the asymptotic limits to provide a geometrically perfect starting position**:

 $$\Large \beta_0 = \frac{\alpha}{\tanh\left( \left(\frac{6}{5}\right)^\alpha \sqrt{\alpha} \right)}$$

 * This is the computational "Holy Grail" realized by Simarro and Orfila (2013). By evaluating highly accurate explicit formulations not as an endpoint, but exclusively as a pre-calculated starting position ($\beta_{\text{seed}}$), they place the algorithm so close to the true root that exactly *one* rigid Newton step drops the error down to an astonishing $0.0000082\%$. 

 $$\Large \beta_{\text{final}} = \beta_{\text{seed}} - \frac{\alpha - \beta_{\text{seed}} \tanh(\beta_{\text{seed}})}{-\tanh(\beta_{\text{seed}}) - \beta_{\text{seed}} \text{sech}^2(\beta_{\text{seed}})}$$

 * **Because ocean wave models (like SWAN, WaveWatch III, or Boussinesq phase-resolving grids) integrate phase celerity over millions of time-steps, even a $0.1\%$ error in $k$ accumulates into massive, non-physical phase shifts, completely destroying wave interference patterns and refraction corridors**. This hybrid "seed-and-step" topology guarantees analytical exactness while retaining the vectorized speed necessary for global forecasting.

Through this suite of algorithms, hydrodynamical researchers are comprehensively equipped to balance processing execution latency against the uncompromising demands of numerical precision, ensuring robust wave phase and group celerity integration across every conceivable spectral forecasting environment. **Whether modeling short-crested capillary chop in a coastal flume or resolving deep-ocean rogue wave non-linearities, the `wave-disp-equation.py` repository provides the exact mathematical key to unlocking the dispersion relation**.

---

## Most Relevant Formulas Ranked

| Model Name | Equation Format | MaxErr |
| :--- | :--- | :--- |
| Pade(2025)_8 | High-order fractional Padé polynomial | $0.0000000\%$ |
| Pade(2025)_11 | High-order fractional Padé polynomial | $0.0000000\%$ |
| Simarro&Orfila(2013) | $\beta_{\text{Beji}} - f(\beta_{\text{Beji}})/f'(\beta_{\text{Beji}})$ | $0.0000082\%$ |
| Vatankhah(2013)_1 | Split-form exponential corrective growth | $0.00176\%$ |
| Pade(2025)_3 | $N_{4.5}(\alpha) / D_{4.0}(\alpha)$ | $0.00666\%$ |
| Hunt(1979)_9 | $[\alpha^2 + \alpha / (1 + \sum_{n=1}^9 d_n \alpha^n)]^{1/2}$ | $0.00815\%$ |
| Vatankhah(2013)_2 | $[\alpha + \alpha^2 \exp(-1.835 - 1.225\alpha^{1.35})] / \sqrt{\tanh(\alpha)}$ | $0.01889\%$ |
| Wu&Thornton(1986) | Piecewise shallow/deep bounds | $0.03427\%$ |
| Beji(2013) | $\alpha \exp(-\beta_0 + \beta_1\alpha + \beta_2\alpha^2)$ | $0.04423\%$ |
| Carvalho(2025)_4 | Deep continuous fraction | $0.050\%$ |
| Carvalho(2025)_5 | $\alpha / \tanh(1.199^{\alpha^{1.047}} \alpha^{0.499})$ | $0.076\%$ |
| Hunt(1979)_5 | $[\alpha^2 + \alpha / (1 + 0.652\alpha + \dots + 0.067\alpha^5)]^{1/2}$ | $0.078\%$ |
| Pade(2025)_2 | $N_{3.5}(\alpha) / D_{3.0}(\alpha)$ | $0.102\%$ |
| Carvalho(2025)_9 | $\alpha / \{(\tanh \alpha)^{1/4} \sqrt{\tanh(\sqrt{\sinh \alpha})}\}$ | $0.204\%$ |
| Carvalho(2025)_10 (Seed)| $\alpha / \tanh((6/5)^\alpha \alpha^{1/2})$ | $0.271\%$ |
| Nielsen(1982) | Piecewise Taylor series / exponential | $0.580\%$ |
| You(2002) | Piecewise Taylor series / exponential | $0.580\%$ |
| Pade(2025)_1 | $N_{2.5}(\alpha) / D_{2.0}(\alpha)$ | $0.649\%$ |
| Yamaguchi(2007)_2 (YN2) | Nested FM substitution $(m=1.378)$ | $0.732\%$ |
| Guo(2002) | $\alpha / \{ 1 - \exp(-\alpha^{m/2}) \}^{1/m} \quad (m=2.4901)$ | $0.757\%$ |
| Carvalho(2025)_18 | $\alpha / \tanh(\sinh \alpha^{1/2})$ | $1.129\%$ |
| Yamaguchi(2007)_1 (YN1) | $\alpha(\coth \alpha^{m/2})^{1/m} \quad (m=1.485)$ | $1.543\%$ |
| Fenton&McKee(1990)_2 | $\alpha / (\tanh \alpha^{3/4})^{2/3}$ | $1.631\%$ |
| Iwagaki(2007) | $\alpha / \tanh\{\alpha^{1/2}(1 + \alpha^{1/2}/2\pi)\}$ | $3.147\%$ |
| Yamaguchi(2007)_4 | $\alpha(1 + \alpha^{-2})^{1/4}$ | $3.177\%$ |
| Eckart(1951) | $\alpha / (\tanh \alpha)^{1/2}$ | $4.980\%$ |

---

## Glossary of Variables, Parameters, and Terms 

This glossary lists physical variables, mathematical parameters, algorithmic concepts, and theoretical terms used in this document.

## I. Fundamental Physical & Spatial Variables

* **$x$**: The primary horizontal spatial coordinate. It defines the geographic direction of progressive wave propagation across the fluid domain.
* **$z$**: The vertical spatial coordinate. It is constrained by the flat impermeable seabed at the bottom and the moving fluid free-surface at the top, with the origin ($z = 0$) strictly locked to the mean water level.
* **$t$**: The temporal parameter (time). It is essential for evaluating the unsteady, transient kinematics of surface gravity waves and computing angular frequencies.
* **$h$**: The bathymetric water depth. It represents the absolute vertical distance from the mean water level ($z = 0$) to the solid seabed ($z = -h$). It is typically acquired via physical hydrographic surveys and dictates the precise depth-regime of the wave.
* **$\eta(x,t)$**: The free-surface elevation. It represents the instantaneous vertical displacement of the fluid boundary. For a monochromatic harmonic wave, it is defined mathematically as $\eta(x,t) = a \cos(kx - \omega t)$.
* **$a$**: The wave amplitude. It defines the maximum vertical elevation of the wave crest above the mean water level. Under linear Airy theory, it is assumed to be infinitesimally small compared to the wavelength and depth.
* **$T$**: The wave period. It measures the time required for one complete wave cycle to pass a stationary reference point. It is typically empirical, gathered via offshore wave buoys or global spectral wind-wave models.
* **$\lambda$ (or $L$)**: The physical wavelength. It is the spatial distance between two consecutive wave crests. 
* **$L_0$**: The deep-water wavelength. Defined exclusively by the wave period and gravity as $L_0 = gT^2/2\pi$, representing the unshoaled length of a wave in the open ocean.
* **$g$**: The constant of gravitational acceleration. It acts as the primary downward hydrostatic restoring force that propagates surface gravity waves.
* **$\rho$**: The primary fluid density (typically seawater). It is assumed to be constant under the Navier-Stokes incompressibility constraint.
* **$\rho'$**: The density of a secondary, usually lighter, fluid layer (e.g., a buoyant freshwater estuarine plume). The density gradient between $\rho$ and $\rho'$ drives the propagation of internal interfacial waves.
* **$\gamma$ (or $\sigma_t$)**: Fluid surface tension. A cohesive microscopic force that supersedes gravity at wavelengths shorter than 0.07 meters, altering the dispersion relation for capillary waves.

## II. Derived Kinematic & Hydrodynamic Variables

* **$\Phi(x, z, t)$**: The scalar velocity potential. Enabled by the assumption of irrotational flow, it is a single mathematical scalar field whose gradient yields the entire fluid velocity vector.
* **$\mathbf{u} = (u, w)$**: The two-dimensional fluid velocity vector. It is composed of the horizontal velocity component $u$ and the vertical velocity component $w$, calculated as $\mathbf{u} = \nabla \Phi$.
* **$k$**: The angular wavenumber. Mathematically defined as $k = 2\pi/\lambda$, it is the spatial frequency parameter that must be solved for. It completely dictates all downstream dynamics, including wave energy density, momentum flux, and refraction coefficients.
* **$k_0$**: The deep-water angular wavenumber. Defined explicitly as $k_0 = \omega^2/g$, it represents the wavenumber strictly in the asymptotic deep-water limit.
* **$\omega$**: The absolute angular frequency. Defined as $\omega = 2\pi/T$, it is the temporal driving parameter of the wave system observed from a stationary geographic frame.
* **$c_p$**: The phase velocity or celerity. Calculated as $c_p = \omega/k$, it is the speed at which the physical shape of the wave (the crest) propagates through the domain.
* **$c_g$**: The group velocity. Calculated as the derivative $c_g = \partial\omega/\partial k$, it is the physical speed at which wave energy and macroscopic wave packages travel across the ocean surface.
* **$\mathbf{V}$**: The depth-averaged velocity vector of an ambient mean fluid current. It is required to compute Galilean transformations in wave-current interaction fields.
* **$\sigma$**: The intrinsic relative frequency. It represents the Doppler-shifted wave frequency observed by an entity floating along with the underlying mean fluid current $\mathbf{V}$.
* **$E$**: Wave energy density. Calculated as $E = \frac{1}{2}\rho g a^2$, it dictates the raw kinetic and potential energy carried by the wave.

## III. Dimensionless Mathematical Parameters & Algorithm Constants

* **$\alpha$**: The independent dimensionless deep-water parameter. It is defined strictly as $\alpha = k_0h$ and serves as the primary standardized input argument for every explicit algorithm in the repository.
* **$\beta$**: The dependent dimensionless local wavenumber. It is defined as $\beta = kh$ and represents the target output variable that algorithms map from $\alpha$.
* **$\beta_n, \beta_{n+1}$**: Iterative states. They represent the current ($\beta_n$) and subsequently updated ($\beta_{n+1}$) estimations of the root during a Newton-Raphson loop.
* **$\beta_0$ (or $\beta_{\text{seed}}$)**: The initial algorithmic guess. It is the starting position provided to the Newton-Raphson solver to bypass geometric instability. 
* **$\beta_{\text{Beji}}$**: A specific, highly accurate explicit estimate calculated using the Beji (2013) formula, utilized purely as the $\beta_0$ seed for the Simarro and Orfila single-step update.
* **$m$**: A tunable generalized exponent. Extensively utilized in minimax error optimizations by Fenton & MacKee ($m=1.5$), Yamaguchi & Nonaka ($m=1.485$), and Guo ($m=2.4901$) to tighten algebraic curves across transitional depths.
* **$\beta_F$**: A temporary nested variable utilized exclusively in the Yamaguchi and Nonaka (YN2) modification to simulate a quasi-Newton corrective step.
* **$d_n$ ($d_1$ through $d_9$)**: Precisely tuned analytical coefficients utilized within Hunt's 9th Order Padé denominator polynomial expansion.
* **$N_3(\alpha), D_3(\alpha)$**: The distinct fractional polynomial numerator and denominator arrays formulated in Carvalho's 2025 Padé functions.

## IV. Theoretical Concepts & Physical Phenomena

* **Linear Wave Theory (Airy Theory)**: The fundamental mathematical framework describing wave propagation over a fluid layer. It relies on linearizing boundary conditions by assuming infinitesimal wave steepness ($ak \ll 1$).
* **Laplace Equation ($\nabla^2 \Phi = 0$)**: The core linear elliptic partial differential equation derived from mass conservation (incompressibility) that the velocity potential must satisfy globally.
* **Irrotational Flow ($\nabla \times \mathbf{u} = 0$)**: A simplifying fluid dynamic assumption asserting the absence of macroscopic vorticity within the water column, allowing velocity to be mapped from a scalar potential.
* **Kinematic Bed/Free-Surface Boundary Conditions**: Strict bounding rules. The bed condition dictates zero vertical fluid penetration at the seabed, while the free-surface condition dictates that a particle on the surface must track with the surface's movement.
* **Dynamic Free-Surface Boundary Condition**: A pressure constraint applying the unsteady Bernoulli equation to ensure fluid surface pressure precisely matches atmospheric pressure.
* **Transcendental Equation**: An equation (like $\alpha = \beta \tanh \beta$) that cannot be resolved for a specific variable using standard inverse algebra, because the variable acts both as a linear multiplier and the argument of a non-algebraic function.
* **Deep Water Limit ($\alpha \rightarrow \infty$)**: The physical asymptotic boundary where depth greatly exceeds wavelength, causing $\tanh \beta \rightarrow 1$. Waves become purely dispersive, and speed is completely independent of depth.
* **Shallow Water Limit ($\alpha \rightarrow 0$)**: The physical asymptotic boundary where waves drag heavily on the seabed, causing $\tanh \beta \rightarrow \beta$. Waves become completely non-dispersive, and speed is strictly dependent on depth ($c_p = \sqrt{gh}$).
* **Amplitude Dispersion**: A third-order non-linear phenomenon where phase celerity becomes dependent on wave steepness via a $[1 + (ka)^2]$ multiplier, causing large crests to travel faster than small ones.
* **Benjamin-Feir Modulational Instability**: A direct chaotic consequence of amplitude dispersion driving nonlinear wave-wave interactions, precipitating extreme open-ocean rogue waves.
* **Hydrodynamic Doppler Shifting**: The Galilean phase shift observed when wave trains propagate across a moving background current, compressing or expanding the intrinsic frequency.
* **Wave Blocking**: A catastrophic energy accumulation event occurring when an adverse fluid current flows at a velocity strong enough to physically prevent the forward propagation of wave energy.
* **Pycnocline**: A subsurface density interface separating distinct fluid layers, across which internal waves slowly propagate with colossal vertical amplitudes.

## V. Computational & Algorithmic Methodologies

* **Newton-Raphson Implementation**: A standard iterative root-finding numerical method utilized to converge upon double-precision "mathematical truth" by advancing via tangent line projections.
* **Residual Function ($f(\beta)$)**: The objective equation rearranged to equal zero. In this module, it is defined as $f(\beta) = \alpha - \beta \tanh \beta$.
* **Jacobian ($f'(\beta)$)**: The analytical first derivative of the residual function with respect to $\beta$, used as the denominator in the Newton-Raphson step.
* **Stability Topography & Inflection Point**: A geometric hazard in the residual function occurring specifically at $\beta \approx 1.200$ (where $f''(\beta) = 0$). Tangent slopes here become momentarily flat, throwing subsequent estimates into divergence.
* **Asymptotic Matching**: An advanced mathematical technique used to engineer formulas that flawlessly decay into the exact theoretical boundaries of both the shallow water limit ($\sqrt{\alpha}$) and deep water limit ($\alpha$).
* **Gene Expression Programming (GEP)**: A high-performance evolutionary artificial intelligence algorithm utilizing linear chromosomes and branched expression trees to discover highly complex, non-intuitive algebraic relationships.
* **Padé Approximant**: A highly robust rational function built from the ratio of two separate polynomials. It is mathematically superior to Taylor expansions for simulating poles, asymptotes, and hyperbolic tangent saturation geometries.
* **Minimax Error Optimization**: A computational curve-fitting strategy designed to equally balance and minimize the absolute maximum positive and negative relative errors across a formula's domain.
* **Seed-and-Step Paradigm**: A contemporary optimization topology. It evaluates a highly accurate explicit formulation to place the algorithmic state geometrically close to the root, followed by exactly one deterministic Newton update to achieve sub-machine-epsilon precision without looping.
* **Warp Divergence**: A severe computational thread-stalling penalty occurring on highly parallelized GPU tensor architectures (CUDA/OpenCL) caused by evaluating iterative conditional logic (like `while` loops).
* **Piecewise Conditional Logic**: The practice of splitting the domain into multiple formulas (e.g., one for shallow, one for deep). It introduces dangerous mathematical discontinuities ("jumps") at the boundary that corrupt phase speed grids.
* **Maclaurin Series Expansion**: A specific mathematical Taylor series expansion centered at zero, leveraged primarily to derive restricted-domain approximations for the nearshore shallow water zone.

---

## References

1. **Wikipedia**. "Airy wave theory". [Wikipedia](https://en.wikipedia.org/wiki/Airy_wave_theory).
2. **Wikipedia**. "Dispersion (water waves)". [Wikipedia](https://en.wikipedia.org/wiki/Dispersion_(water_waves)).
3. **Yu, J.** (2014). "A Note on Approximations of the Dispersion Relationship of Water Waves", *Journal of Engineering Mechanics (ASCE)*, 140(1), 233–237.
4. **Simarro, G. & Orfila, A.** (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves: Another discussion", *Coastal Engineering*, 80, 15–16.
5. **Vatankhah, A.R. & Aghashariatmadari, Z.** (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves: a discussion", *Coastal Engineering*, 78, 21–22.
6. **Beji, S.** (2013). "Improved explicit approximation of linear dispersion relationship for gravity waves", *Coastal Engineering*, 73, 11–12.
7. **You, Z.J.** (2008). "A close approximation of wave dispersion relation for direct calculations",
 *Applied Ocean Research*, 30(2), 141–143.
8. **Yamaguchi, M. & Nonaka, H.** (2007). "Comparative Study of Explicit Solutions to Wave Dispersion Equation", *Journal of JSCE (Ocean Engineering)*, 63(1), 53–66.
9. **Yamaguchi, M. and H. Nonaka.** Comparative study of explicit solutions to wave dispersion equation, *Annu. Jour. Eng.*, Ehime Univ., Vol. 6, 2007 in CD-ROM.
10. **Carvalho, R.** (2006). Unpublished work based on gene expression programming for wave dispersion equations.
11. **You, Z.J.** "Discussion of 'Simple and explicit solution to the wave dispersion equation'", *Coastal Engineering 45 (2002) 71-74*, Coastal Eng., Vol. 48, pp. 133-135, 2003.
12. **Guo, J.** (2002). "Simple and explicit solution of the wave dispersion equation", *Coastal Engineering*, 45, 71–74.
13. **Fenton, J.D. & McKee, W.D.** (1990). "On calculating the lengths of water waves", *Coastal Engineering*, 14, 499–513.
14. **Fenton, J.D.** "The numerical solution of steady water wave problems", *Computers & Geosciences*, Vol. 4, No. 3, pp. 357–368, 1988.
15. **Fenton, J.D.** (1972). "A ninth-order solution for the solitary wave", *Journal of Fluid Mechanics*, 53, 257–271.
16. **Wu, C. S. and E. B. Thornton.** "Wave numbers of linear progressive waves", *Journal of Waterway, Port, Coastal and Ocean Engineering*, ASCE, Vol. 112, No. 4, pp. 536–540, 1986.
17. **Nielsen, P.** "Explicit solutions to practical wave problems", *Proc. 19th ICCE*, Vol. 1, pp. 968–982, 1984.
18. **Nielsen, P.** "Explicit formulae for practical wave calculations", *Coastal Eng.*, No. 6, pp. 389–398, 1982.
19. **Hunt, J.N.** (1979). "A simple approximation for the dispersion relation of water waves", *Journal of Waterway, Port, Coastal and Ocean Engineering*, 105(4), 457–459.
20. **Hunt, J. N.** "Direct solution of wave dispersion equation", *J. Waterway, Port, Coastal and Ocean Div.*, Proc. ASCE, Vol. 105, No. WW4, pp. 457–459, 1979.
21. **Eckart, C.** (1951). "The propagation of gravity waves from deep to shallow water", U.S. Department of Commerce, National Bureau of Standards Circular 521.

These references span classical methods (e.g., Eckart, 1951; Hunt, 1979), modern explicit approximations (e.g., Guo, 2002; Beji, 2013; Vatankhah & Aghashariatmadari, 2013; Simarro & Orfila, 2013; Yu, 2014), pivotal contributions from Fenton and colleagues (Fenton & McKee, 1990; Fenton, 1972), and innovative computational approaches (Carvalho, 2006 & 2025) to improve the dispersion relation accuracy.