# Wave Dispersion Equation – Linear Gravity Waves (Airy Theory)

This repository provides a comprehensive suite of solutions and approximations for analyzing wave dispersion, which is essential for wave prediction, oceanographic calculations, and coastal engineering design. It includes a reference "exact" solution using the Newton-Raphson method, classical and contemporary explicit approximations, and high-order Padé approximants for high precision.
![wave-disp-equation1](https://github.com/user-attachments/assets/7688bf2e-b6f1-4c7b-96c9-37b3078edb3e)
## Background

The **linear wave dispersion equation** relates wave frequency (or period) to wavenumber and water depth for gravity waves:

- **ω (omega)** is the angular frequency:  
  ω = 2π / T  
  (where T is the wave period)

- **g** is the gravitational acceleration.

- **k** is the wavenumber:  
  k = 2π / L  
  (where L is the wavelength)

- **h** is the water depth.

his equation is solved iteratively for the dimensionless wavenumber **kh** (where k_0 = ω²/g). An accurate **kh** evaluation is vital for computing wave phase speed, group velocity, and understanding various nearshore processes. Explicit approximations can bypass the need for iteration but must be chosen carefully based on accuracy requirements.
![wave-disp-equation2](https://github.com/user-attachments/assets/5dcadd0b-223b-461a-ab0f-b99fc3f55982)
## Module Contents

- **Reference "Exact" Solution:** `kh_numeric()` implements the Newton-Raphson iteration method for a highly precise solution of wave dispersion, acting as a benchmark for other techniques.

- **Classical Approximations:** Established methods from researchers like Hunt, Eckart, Nielsen, and Gilbert.

- **Contemporary Approximations:** Recent techniques from researchers such as Guo, Beji, Vatankhah & Aghashariatmadari, Simarro & Orfila, Yu, Fenton & McKee, Guan & Ju, and Iwagaki.

- **High-Order Padé Approximations:** Carvalho's 2025 high-order Padé approximants deliver exceptional precision, addressing increasing complexity in wave calculations as a robust alternative to simpler methods.

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
