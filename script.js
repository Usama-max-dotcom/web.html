// script.js
// Full self-contained quiz: 10 levels × 20 questions = 200 questions
// Config
const QUESTIONS_PER_LEVEL = 20;
const TIME_PER_LEVEL_SECONDS = 600; // 10 minutes

// Full question bank: 10 levels each with 20 math/physics questions
// Each question: { question: string, options: [4 strings], answer: index (0-3) }
const allLevels = [
  // -------------------- LEVEL 1 --------------------
  [
    { question: "Q1: If v = u + at, and u = 5 m/s, a = 2 m/s², t = 3 s, find v.", options: ["11 m/s", "10 m/s", "9 m/s", "6 m/s"], answer: 0 },
    { question: "Q2: Evaluate the definite integral ∫₀³ 2x dx.", options: ["9", "6", "12", "3"], answer: 0 },
    { question: "Q3: A projectile fired at 30° with speed 20 m/s. Using g = 10 m/s², what is the max height?", options: ["5 m", "10 m", "15 m", "20 m"], answer: 0 },
    { question: "Q4: If F = ma, m = 4 kg and a = 3 m/s², what is F?", options: ["7 N", "12 N", "9 N", "3 N"], answer: 1 },
    { question: "Q5: Compute d/dx[ln(x² + 1)] at x = 1.", options: ["1", "2", "1/2", "0"], answer: 0 },
    { question: "Q6: Smallest root of x² - 4x + 3 = 0 is:", options: ["1", "3", "-1", "0"], answer: 0 },
    { question: "Q7: Evaluate ∫₀¹ x³ dx.", options: ["1/3", "1/4", "1/2", "1/5"], answer: 1 },
    { question: "Q8: Block slides from height h on frictionless slope. Speed at bottom is:", options: ["√(gh)", "√(2gh)", "√(g/h)", "2gh"], answer: 1 },
    { question: "Q9: Period of y = sin(2x) is:", options: ["2π", "π", "π/2", "4π"], answer: 1 },
    { question: "Q10: Find d/dx [x² sin x] at x = 0.", options: ["0", "1", "2", "sin 0"], answer: 0 },
    { question: "Q11: Solve e^{2x} = e^4. What is x?", options: ["1", "2", "4", "0"], answer: 1 },
    { question: "Q12: Two identical resistors R in parallel → equivalent resistance:", options: ["R", "2R", "R/2", "R/4"], answer: 2 },
    { question: "Q13: Compute lim_{x→0} (1+x)^{1/x}.", options: ["1", "e", "∞", "0"], answer: 1 },
    { question: "Q14: Period T for angular frequency ω is:", options: ["T = 2π/ω", "T = ω/2π", "T = π/ω", "T = 1/ω"], answer: 0 },
    { question: "Q15: Photon energy E = hc/λ. If λ → 2λ, new energy is:", options: ["2E", "E/2", "E²", "√E"], answer: 1 },
    { question: "Q16: Positive root of x² + x - 6 = 0 is:", options: ["-3", "2", "3", "-2"], answer: 1 },
    { question: "Q17: d/dx[(sin x)/x] at x = π equals:", options: ["-1/π", "1/π", "0", "-π"], answer: 0 },
    { question: "Q18: If vector a ⟂ b, then a·b equals:", options: ["0", "1", "|a||b|", "cannot say"], answer: 0 },
    { question: "Q19: Angular frequency for mass-spring (m,k) is:", options: ["√(k/m)", "k/m", "√(m/k)", "m/k"], answer: 0 },
    { question: "Q20: Evaluate ∫₀^{π/2} cos x dx.", options: ["0", "1", "2", "√2"], answer: 1 }
  ],

  // -------------------- LEVEL 2 --------------------
  [
    { question: "Q1: Solve for x: ln x = 2.", options: ["x = 2", "x = e²", "x = 1/2", "x = 0"], answer: 1 },
    { question: "Q2: Sum of infinite geometric series with a=3, r=1/3 is:", options: ["9/2", "4.5", "3", "3/2"], answer: 0 },
    { question: "Q3: Work done by constant force F over displacement d is:", options: ["F·d", "F/d", "F+d", "Fd²"], answer: 0 },
    { question: "Q4: Derivative of cos x is:", options: ["sin x", "-sin x", "cos x", "-cos x"], answer: 1 },
    { question: "Q5: Solve quadratic x² - 2x -8 = 0, positive root:", options: ["-2", "4", "2", " -4"], answer: 1 },
    { question: "Q6: Magnetic force on charge q moving v in B is q(v×B). Units:", options: ["N", "J", "V", "T"], answer: 0 },
    { question: "Q7: If f(x)=x³, f'(2) is:", options: ["8", "4", "12", "6"], answer: 0 },
    { question: "Q8: Electric field of point charge falls as:", options: ["1/r", "1/r²", "1/r³", "r²"], answer: 1 },
    { question: "Q9: Acceleration in uniform circular motion with speed v and radius r:", options: ["v²/r", "v/r", "r/v", "v²r"], answer: 0 },
    { question: "Q10: Integral ∫ x e^x dx equals:", options: ["xe^x - e^x + C", "e^x + C", "xe^x + C", "e^{x}x²/2 + C"], answer: 0 },
    { question: "Q11: Distance between roots of x² - 6x + 8 = 0:", options: ["2", "sqrt(2)", "4", "0"], answer: 0 },
    { question: "Q12: Snell's law: n1 sinθ1 = n2 sinθ2. If n2>n1, transmitted angle is:", options: ["larger", "smaller", "equal", "zero"], answer: 1 },
    { question: "Q13: Entropy change ΔS for reversible heat Q at T is:", options: ["Q/T", "QT", "Q+T", "T/Q"], answer: 0 },
    { question: "Q14: If matrix A = [[2,0],[0,3]], det(A) =", options: ["6", "5", "1", "0"], answer: 0 },
    { question: "Q15: For y=ln(x), integral ∫ ln x dx = ", options: ["x ln x - x + C", "ln x + C", "x ln x + C", "x + C"], answer: 0 },
    { question: "Q16: If sum of probabilities = 1, P(A)=0.2, P(B)=0.5, other outcomes sum to:", options: ["0.3", "0.2", "0.5", "0.7"], answer: 0 },
    { question: "Q17: Kinetic energy of mass m with speed v is:", options: ["1/2 mv²", "mv", "mv²", "1/3 mv²"], answer: 0 },
    { question: "Q18: If f(x)=sin x, the Maclaurin series first term after 0 is:", options: ["x", "x²", "x³", "1"], answer: 0 },
    { question: "Q19: If two waves are 180° out of phase they:", options: ["interfere destructively", "constructive add", "no interference", "double amplitude"], answer: 0 },
    { question: "Q20: RMS of continuous set of values equal to sqrt(mean(square)). For sin wave amplitude A, RMS is:", options: ["A", "A/√2", "A/2", "√A"], answer: 1 }
  ],

  // -------------------- LEVEL 3 --------------------
  [
    { question: "Q1: Solve ∫₀^{1} (3x²) dx.", options: ["1", "3/2", "1/2", "3"], answer: 0 },
    { question: "Q2: If angular momentum L = Iω, for a solid cylinder I = (1/2)MR², L for cylinder is:", options: ["(1/2)MR² ω", "MR²ω", "MRω", "(1/4)MR² ω"], answer: 0 },
    { question: "Q3: Derivative of tan x is:", options: ["sec² x", "csc² x", "sec x tan x", "cos² x"], answer: 0 },
    { question: "Q4: If acceleration due to gravity is 9.8 m/s², free-fall distance in 2 s is:", options: ["19.6 m", "9.8 m", "39.2 m", "4.9 m"], answer: 0 },
    { question: "Q5: Area of triangle with base b and height h is:", options: ["(1/2) b h", "b h", "bh/3", "2bh"], answer: 0 },
    { question: "Q6: Divergence of constant vector field is:", options: ["0", "1", "depends", "infinite"], answer: 0 },
    { question: "Q7: For series 1 + 1/4 + 1/9 + ... the nth term behaves like:", options: ["1/n²", "1/n", "1/√n", "n²"], answer: 0 },
    { question: "Q8: If pH = 3, [H+] is:", options: ["10^-3 M", "10^3 M", "3 M", "1/3 M"], answer: 0 },
    { question: "Q9: For capacitor C, energy stored is:", options: ["(1/2)CV²", "CV", "C/V", "V/2C"], answer: 0 },
    { question: "Q10: For harmonic oscillator, restoring force is proportional to:", options: ["displacement", "velocity", "acceleration", "time"], answer: 0 },
    { question: "Q11: Solve quadratic 3x² - 12 = 0 positive x is:", options: ["2", "√(4)", "√(12)", " -2"], answer: 0 },
    { question: "Q12: Dot product of (1,2,3)·(4,-5,6) equals:", options: ["12", "4", "20", "0"], answer: 2 },
    { question: "Q13: If wave speed v = fλ, doubling f with fixed v halves:", options: ["λ", "v", "f", "amplitude"], answer: 0 },
    { question: "Q14: Euler's formula e^{iπ} + 1 = 0, so e^{iπ} = ?", options: ["-1", "1", "i", "0"], answer: 0 },
    { question: "Q15: If a resistor 10Ω and 5Ω in series, total is:", options: ["15Ω", "5Ω", "2Ω", "50Ω"], answer: 0 },
    { question: "Q16: For a parabola y = x², slope at x=3 is:", options: ["6", "3", "9", "2"], answer: 0 },
    { question: "Q17: Integrate ∫ cos x dx = ", options: ["sin x + C", "-sin x + C", "cos x + C", "-cos x + C"], answer: 0 },
    { question: "Q18: If change in internal energy ΔU = Q - W, it's first law of thermodynamics; W is:", options: ["work done by system", "heat", "work done on system", "energy lost"], answer: 0 },
    { question: "Q19: If two vectors equal magnitude but opposite direction, resultant is:", options: ["0", "2 times", "magnitude", "undefined"], answer: 0 },
    { question: "Q20: For function f(x)=e^{3x}, f'(x) =", options: ["3e^{3x}", "e^{x}", "e^{3x}/3", "9e^{3x}"], answer: 0 }
  ],

  // -------------------- LEVEL 4 --------------------
  [
    { question: "Q1: Solve ∫ (1/x) dx = ", options: ["ln|x| + C", "x + C", "1/x + C", "e^x + C"], answer: 0 },
    { question: "Q2: For projectile, range R = (v² sin 2θ)/g. Max range when θ = ?", options: ["45°", "30°", "60°", "90°"], answer: 0 },
    { question: "Q3: Determinant of [[1,2],[3,4]] is:", options: ["-2", "2", "10", "0"], answer: 0 },
    { question: "Q4: For uniform rod rotating about center, moment of inertia I ∝:", options: ["MR²", "M/R", "M", "R"], answer: 0 },
    { question: "Q5: If acceleration is derivative of velocity, integrate acceleration to get:", options: ["velocity", "position", "force", "mass"], answer: 0 },
    { question: "Q6: Solve equation cos x = 0 at smallest positive x:", options: ["π/2", "π", "π/3", "2π"], answer: 0 },
    { question: "Q7: Sum of first n natural numbers = ", options: ["n(n+1)/2", "n²", "n(n-1)/2", "n/2"], answer: 0 },
    { question: "Q8: Thermal expansion coefficient α relates ΔL = αLΔT, if α small increase is:", options: ["linear with ΔT", "quadratic", "inversely", "logarithmic"], answer: 0 },
    { question: "Q9: If I = current, V = IR, then R = ", options: ["V/I", "I/V", "VI", "V+I"], answer: 0 },
    { question: "Q10: Area under curve y = x from 0 to 2 is:", options: ["2", "1", "4", "0"], answer: 0 },
    { question: "Q11: For box moving constant velocity, net force =", options: ["0", "mass×acceleration", "weight", "friction"], answer: 0 },
    { question: "Q12: The Laplace transform of 1 is:", options: ["1/s", "s", "0", "e^s"], answer: 0 },
    { question: "Q13: Covalent bond formed by:", options: ["sharing electrons", "transfer electrons", "metallic bonding", "ionic"], answer: 0 },
    { question: "Q14: If wavefunction normalized, total probability = ", options: ["1", "0", "infinite", "depends"], answer: 0 },
    { question: "Q15: Sine series sin(nπx) are orthogonal on interval", options: ["0 to 1", "-1 to 1", "0 to π", "π to 2π"], answer: 0 },
    { question: "Q16: If f(x)=x^n, ∫ f'(x) dx = ", options: ["x^n + C", "nx^{n-1}", "nx^n", "n+x"], answer: 0 },
    { question: "Q17: If an ideal gas doubles T at constant volume, pressure:", options: ["doubles", "halves", "unchanged", "quadruples"], answer: 0 },
    { question: "Q18: For harmonic motion x(t)=A cos(ωt), max velocity is:", options: ["Aω", "A/ω", "ω/A", "Aω²"], answer: 0 },
    { question: "Q19: If rod length L, moment of inertia about center for thin rod is:", options: ["(1/12) mL²", "(1/2) mL²", "mL²", "(1/3) mL²"], answer: 0 },
    { question: "Q20: For function f(x)=1/(1+x²), integral from -∞ to ∞ is π. True or false?", options: ["True", "False", "Only finite", "Needs constant"], answer: 0 }
  ],

  // -------------------- LEVEL 5 --------------------
  [
    { question: "Q1: If f(x)=arctan x, derivative f'(x) =", options: ["1/(1+x²)", "1/x", "x/(1+x²)", " -1/(1+x²)"], answer: 0 },
    { question: "Q2: Kinetic energy of rotating solid cylinder: (1/2)Iω². If I=(1/2)MR², KE =", options: ["(1/4)MR²ω²", "(1/2)MR²ω²", "(1/8)MR²ω²", "MR²ω²"], answer: 0 },
    { question: "Q3: If cosh x = (e^x + e^{-x})/2, cosh 0 equals:", options: ["1", "0", "2", "1/2"], answer: 0 },
    { question: "Q4: The binomial coefficient C(n,k) counts:", options: ["k-combinations", "permutations", "arrangements with order", "functions"], answer: 0 },
    { question: "Q5: Work done in isothermal ideal gas compression from V1 to V2 is:", options: ["nRT ln(V1/V2)", "nRΔT", "PΔV", "0"], answer: 0 },
    { question: "Q6: For a diode forward biased, it conducts when V exceeds:", options: ["~0.7V (Si)", "~5V", "~0V", "~12V"], answer: 0 },
    { question: "Q7: In 2nd order ODE y'' + ω² y = 0, general solution is:", options: ["A cos ωt + B sin ωt", "Ae^{ωt}", "A/t", "constant"], answer: 0 },
    { question: "Q8: If pV = nRT, doubling T at constant p doubles:", options: ["V", "n", "R", "p"], answer: 0 },
    { question: "Q9: For simple pendulum small-angle period T ≈", options: ["2π√(L/g)", "√(L/g)", "2πL/g", "π√(L/g)"], answer: 0 },
    { question: "Q10: Covariance of independent variables is:", options: ["0", "1", "undefined", "infinite"], answer: 0 },
    { question: "Q11: Fourier transform converts time domain to:", options: ["frequency domain", "time domain", "spatial domain", "probability domain"], answer: 0 },
    { question: "Q12: For a uniform sphere moment of inertia is:", options: ["(2/5)MR²", "(1/2)MR²", "(3/5)MR²", "(1/5)MR²"], answer: 0 },
    { question: "Q13: If x solves sin x = 1, smallest positive x is:", options: ["π/2", "π", "2π", "π/4"], answer: 0 },
    { question: "Q14: The eigenvalues of rotation matrix in 2D are:", options: ["complex of unit modulus", "real >1", "real <1", "zero"], answer: 0 },
    { question: "Q15: If reflection coefficient r = (n1-n2)/(n1+n2), for n1=n2, r=", options: ["0", "1", "-1", "2"], answer: 0 },
    { question: "Q16: The RMS speed of ideal gas ∝", options: ["√T", "T", "T²", "1/T"], answer: 0 },
    { question: "Q17: Complex roots of polynomial with real coefficients come in:", options: ["conjugate pairs", "triples", "random", "not related"], answer: 0 },
    { question: "Q18: Sum of probabilities in discrete distribution equals:", options: ["1", "0", "depends", "≥1"], answer: 0 },
    { question: "Q19: If two capacitors C1 and C2 in series, equivalent Ceq =", options: ["1/(1/C1 + 1/C2)", "C1 + C2", "C1C2", "C1/C2"], answer: 0 },
    { question: "Q20: Heat capacity at constant volume Cv for ideal monoatomic gas = ", options: ["(3/2)R", "R", "(5/2)R", "0"], answer: 0 }
  ],

  // -------------------- LEVEL 6 --------------------
  [
    { question: "Q1: If f(x)=cos x, Maclaurin series first terms: 1 - x²/2! + x^4/4! + ... true or false?", options: ["True", "False", "Only odd terms", "Only constants"], answer: 0 },
    { question: "Q2: For non-relativistic particle, de Broglie wavelength λ = h/p, p is:", options: ["momentum", "energy", "mass", "velocity"], answer: 0 },
    { question: "Q3: For RLC circuit at resonance, impedance is:", options: ["minimum", "maximum", "infinite", "zero"], answer: 0 },
    { question: "Q4: Solve derivative of x^x at x=1: d/dx x^x |_{x=1} = ?", options: ["1 + ln1", "1", "0", "undefined"], answer: 1 },
    { question: "Q5: In uniform acceleration, displacement s = ut + 1/2 at²; if u=0, s after t is:", options: ["1/2 at²", "at", "at²", "ut"], answer: 0 },
    { question: "Q6: If eigenvectors of matrix are orthogonal for symmetric matrices: true or false?", options: ["True", "False", "Only diagonal", "Depends on size"], answer: 0 },
    { question: "Q7: If two waves of equal amplitude and phase meet, resultant amplitude is:", options: ["double", "zero", "same", "half"], answer: 0 },
    { question: "Q8: Stefan-Boltzmann law: radiated power ∝", options: ["T^4", "T", "√T", "ln T"], answer: 0 },
    { question: "Q9: If function has derivative zero everywhere, function is:", options: ["constant", "linear", "quadratic", "zero only"], answer: 0 },
    { question: "Q10: For Newton-Raphson root finding, convergence is usually:", options: ["quadratic", "linear", "logarithmic", "no convergence"], answer: 0 },
    { question: "Q11: The gradient points in direction of:", options: ["greatest increase", "greatest decrease", "zero change", "minimum"], answer: 0 },
    { question: "Q12: If a charge q at center of conducting shell, field inside conductor material is:", options: ["zero", "not zero", "depends", "infinite"], answer: 0 },
    { question: "Q13: Poisson distribution mean equals variance: true or false?", options: ["True", "False", "Only for large λ", "Only for λ=1"], answer: 0 },
    { question: "Q14: If y = ln(sin x), derivative y' equals:", options: ["cot x", "tan x", "csc x", "sec x"], answer: 0 },
    { question: "Q15: The Bernoulli equation is valid for:", options: ["incompressible, non-viscous steady flow", "all fluids", "gas flows only", "viscous flows"], answer: 0 },
    { question: "Q16: If an electron accelerates, it emits:", options: ["radiation (if accelerated)", "no radiation", "only gamma", "only heat"], answer: 0 },
    { question: "Q17: Probability density integrates to 1 for normalized wavefunction: true or false?", options: ["True", "False", "Only approximately", "Depends on potential"], answer: 0 },
    { question: "Q18: For a spherical shell, gravitational field outside is as if mass concentrated at:", options: ["center", "surface", "edge", "everywhere"], answer: 0 },
    { question: "Q19: The logistic map x_{n+1}=rx_n(1-x_n) shows chaos for r around:", options: ["3.57", "1", "2", "0.5"], answer: 0 },
    { question: "Q20: Fourier series requires function to be:", options: ["piecewise continuous", "derivative-free", "infinite", "only polynomial"], answer: 0 }
  ],

  // -------------------- LEVEL 7 --------------------
  [
    { question: "Q1: Solve ∫ x ln x dx = ?", options: ["(x²/2) ln x - x²/4 + C", "x ln x - x + C", "ln x + C", "x² ln x + C"], answer: 0 },
    { question: "Q2: If lens has focal length f, and object at distance u, image at v given by:", options: ["1/f = 1/u + 1/v", "f = u + v", "uv = f", "u = v"], answer: 0 },
    { question: "Q3: For damped oscillator, energy decays exponentially in time when damping present: true or false?", options: ["True", "False", "Only for critical damping", "Only for undamped"], answer: 0 },
    { question: "Q4: The eigenvalues of [[2,0],[0,3]] are:", options: ["2 and 3", "5 and -1", "0 and 6", "1 and 1"], answer: 0 },
    { question: "Q5: For Maxwell's equations, ∇·B = 0 means:", options: ["no magnetic monopoles", "no electric charge", "B constant", "magnetic sources exist"], answer: 0 },
    { question: "Q6: If function f has Taylor series around 0, remainder tends to 0 if:", options: ["function analytic", "function linear", "not differentiable", "discontinuous"], answer: 0 },
    { question: "Q7: In projectile motion, vertical and horizontal motions are:", options: ["independent", "dependent", "same", "inverses"], answer: 0 },
    { question: "Q8: The capacitor's reactance Xc = 1/(ωC). If ω increases, Xc:", options: ["decreases", "increases", "unchanged", "zero"], answer: 0 },
    { question: "Q9: For small oscillations potential ≈ quadratic; stable equilibrium means second derivative is:", options: ["positive", "negative", "zero", "infinite"], answer: 0 },
    { question: "Q10: A linear operator T is diagonalizable if:", options: ["has a basis of eigenvectors", "trace is zero", "determinant is 1", "is singular"], answer: 0 },
    { question: "Q11: If a mass m slides on frictionless track, mechanical energy is conserved: true or false?", options: ["True", "False", "Only if potential constant", "Only with friction"], answer: 0 },
    { question: "Q12: Resistivity ρ relates to resistance R by:", options: ["R = ρ L/A", "R = ρ A/L", "R = L/A ρ²", "R = ρ+L/A"], answer: 0 },
    { question: "Q13: For harmonic oscillator, phase difference between displacement and velocity is:", options: ["π/2", "π", "0", "2π"], answer: 0 },
    { question: "Q14: If f(x)=|x|, derivative at 0 exists? ", options: ["No", "Yes", "Only from left", "Only from right"], answer: 0 },
    { question: "Q15: If momentum p conserved in closed system, net external force:", options: ["is zero", "non-zero", "depends", "infinite"], answer: 0 },
    { question: "Q16: The heat equation ∂u/∂t = α ∂²u/∂x² is:", options: ["parabolic PDE", "elliptic", "hyperbolic", "algebraic"], answer: 0 },
    { question: "Q17: Inverse Laplace transform recovers time-domain from:", options: ["s-domain", "time-domain", "frequency-domain", "spatial-domain"], answer: 0 },
    { question: "Q18: For ideal transformer, Vp/Vs = Np/Ns. True or false?", options: ["True", "False", "Only DC", "Only AC"], answer: 0 },
    { question: "Q19: Shannon entropy measures:", options: ["average information per symbol", "energy", "mass", "force"], answer: 0 },
    { question: "Q20: If function g(x)=x/(1+x²), as x→∞ g(x)→", options: ["0", "1", "∞", "-1"], answer: 0 }
  ],

  // -------------------- LEVEL 8 --------------------
  [
    { question: "Q1: Solve ∫ sec² x dx = ?", options: ["tan x + C", "sec x + C", "-tan x + C", "cos x + C"], answer: 0 },
    { question: "Q2: For relativistic energy E = mc², if m doubles, E:", options: ["doubles", "quadruples", "halves", "unchanged"], answer: 0 },
    { question: "Q3: For simple harmonic motion, average potential energy over full cycle equals average kinetic energy: true or false?", options: ["True", "False", "Only at peaks", "Only at equilibrium"], answer: 0 },
    { question: "Q4: In a PN junction, minority carriers are:", options: ["electrons in p and holes in n", "holes in p", "electrons in n only", "ions"], answer: 0 },
    { question: "Q5: For matrix multiplication AB, (AB)^T equals:", options: ["B^T A^T", "A^T B^T", "A B", "B A"], answer: 0 },
    { question: "Q6: Blackbody peak wavelength λ_max ∝ 1/T (Wien's law): true or false?", options: ["True", "False", "∝ T", "no relation"], answer: 0 },
    { question: "Q7: Solve derivative of ln(x^2) = ?", options: ["2/x", "1/x", "x^2", "2x"], answer: 0 },
    { question: "Q8: If system is conservative, curl of force field is:", options: ["0", "non-zero", "infinite", "depends"], answer: 0 },
    { question: "Q9: In diffraction, single-slit minima occur when a sin θ = mλ. True or false?", options: ["True", "False", "only double-slit", "only large angles"], answer: 0 },
    { question: "Q10: The binomial expansion of (1+x)^n holds for real n using series if |x|<1: true or false?", options: ["True", "False", "Only integer n", "Only positive n"], answer: 0 },
    { question: "Q11: For two bodies gravitational force ∝", options: ["1/r²", "r", "r²", "1/r"], answer: 0 },
    { question: "Q12: If probability of event A is 0.6 and B is 0.5, maximum possible P(A∩B) is:", options: ["0.5", "0.6", "1.1", "0.1"], answer: 0 },
    { question: "Q13: For logistic regression output is:", options: ["probability (0-1)", "raw score", "complex", "negative only"], answer: 0 },
    { question: "Q14: If function is odd, integral symmetric about 0 equals:", options: ["0", "2 times", "undefined", "depends"], answer: 0 },
    { question: "Q15: Work-energy theorem states ΔK = W_net: true or false?", options: ["True", "False", "Only conservative", "Only steady"], answer: 0 },
    { question: "Q16: For function f(x)=x^4, f''(0) equals:", options: ["0", "12", "0.5", "4"], answer: 0 },
    { question: "Q17: The cross product a×b is perpendicular to:", options: ["both a and b", "a only", "b only", "none"], answer: 0 },
    { question: "Q18: In quantum mechanics, Heisenberg uncertainty relates Δx Δp ≥:", options: ["ħ/2", "0", "h", "1"], answer: 0 },
    { question: "Q19: For a rod stretched, stress = force/area; strain = ΔL/L. Young's modulus = ", options: ["stress/strain", "strain/stress", "force × length", "area × stress"], answer: 0 },
    { question: "Q20: If differential eqn y' = ky, solution is:", options: ["Ce^{kt}", "kt + C", "k/t", "0"], answer: 0 }
  ],

  // -------------------- LEVEL 9 --------------------
  [
    { question: "Q1: Evaluate ∫₀^{π} sin x dx =", options: ["2", "0", "1", "π"], answer: 0 },
    { question: "Q2: If two springs k1 and k2 in series, equivalent k =", options: ["1/(1/k1 + 1/k2)", "k1 + k2", "k1 k2", "k1/k2"], answer: 0 },
    { question: "Q3: For orthogonal vectors, dot product equals:", options: ["0", "1", "magnitude", "depends"], answer: 0 },
    { question: "Q4: If transform of cos ωt is π[δ(ω-ω0) + δ(ω+ω0)], true or false?", options: ["True", "False", "Only for sin", "Only discrete"], answer: 0 },
    { question: "Q5: The centripetal force needed is m v²/r. If v doubles, required force:", options: ["quadruples", "doubles", "halves", "stays"], answer: 0 },
    { question: "Q6: If A is invertible, det(A^{-1}) equals:", options: ["1/det(A)", "det(A)", "-det(A)", "0"], answer: 0 },
    { question: "Q7: In equipartition theorem each quadratic degree contributes (1/2)kT to average energy: true or false?", options: ["True", "False", "Only classical", "Only quantum"], answer: 0 },
    { question: "Q8: For entropy S, units are:", options: ["J/K", "J", "K", "dimensionless"], answer: 0 },
    { question: "Q9: If function f has discontinuity of first kind, it's piecewise continuous: true or false?", options: ["True", "False", "Undefined", "Depends on derivative"], answer: 0 },
    { question: "Q10: For small oscillations, frequency depends on curvature of potential: true or false?", options: ["True", "False", "Only amplitude", "Only mass"], answer: 0 },
    { question: "Q11: Kirchhoff's voltage law sums voltages around loop to:", options: ["0", "1", "depends", "total current"], answer: 0 },
    { question: "Q12: If polynomial degree n has at most how many real roots?", options: ["n", "n²", "infinite", "1"], answer: 0 },
    { question: "Q13: If binary entropy at p=0.5 it's maximal: true or false?", options: ["True", "False", "depends on base", "no"], answer: 0 },
    { question: "Q14: For capacitance, units are:", options: ["Farads (F)", "Ohms", "Henrys", "Volts"], answer: 0 },
    { question: "Q15: Euler-Lagrange equation yields equations of motion from:", options: ["action principle", "energy minimization", "force balance", "Newton only"], answer: 0 },
    { question: "Q16: Boltzmann constant k relates energy to temperature: units are J/K: true or false?", options: ["True", "False", "only eV", "only gas"], answer: 0 },
    { question: "Q17: For logistic map period-doubling leads to chaos: true or false?", options: ["True", "False", "only linear", "only cubic"], answer: 0 },
    { question: "Q18: If light enters denser medium, wavelength:", options: ["decreases", "increases", "unchanged", "becomes zero"], answer: 0 },
    { question: "Q19: Rank of identity matrix I_n is:", options: ["n", "1", "0", "n-1"], answer: 0 },
    { question: "Q20: For vector field F conservative iff curl F:", options: ["= 0", "≠ 0", "div F = 0", "depends"], answer: 0 }
  ],

  // -------------------- LEVEL 10 --------------------
  [
    { question: "Q1: Solve ∫ 0 to 1 of (1/(1+x^2)) dx = ?", options: ["π/4", "π/2", "1", "ln2"], answer: 0 },
    { question: "Q2: If Schrödinger equation is linear, superposition principle applies: true or false?", options: ["True", "False", "Only classical", "Only stationary"], answer: 0 },
    { question: "Q3: For special relativity time dilation: t = γ t0, γ = 1/√(1-(v²/c²)). If v→0 then γ→", options: ["1", "0", "∞", "2"], answer: 0 },
    { question: "Q4: Blackbody Stefan-Boltzmann constant σ appears in j = σT^4. Units of σ include:", options: ["W·m^{-2}·K^{-4}", "J/kg", "m/s", "W"], answer: 0 },
    { question: "Q5: If two mirrors face each other, multiple reflections produce cavities called:", options: ["Fabry–Pérot", "Michelson", "Wien", "Young"], answer: 0 },
    { question: "Q6: For complex number z = re^{iθ}, log z has branch cut: true or false?", options: ["True", "False", "Only real z", "No"], answer: 0 },
    { question: "Q7: If a particle confined in infinite potential well width L, ground-state energy ∝", options: ["1/L²", "L", "L²", "1/L"], answer: 0 },
    { question: "Q8: In Fourier transform, convolution in time domain corresponds to:", options: ["multiplication in frequency domain", "convolution in frequency", "addition", "subtraction"], answer: 0 },
    { question: "Q9: Pauli exclusion principle applies to fermions: true or false?", options: ["True", "False", "only bosons", "only photons"], answer: 0 },
    { question: "Q10: If Maxwell-Boltzmann distribution used for ideal gas at high temp and low density: true or false?", options: ["True", "False", "only quantum", "only solid"], answer: 0 },
    { question: "Q11: If two beams are coherent, interference pattern requires:", options: ["constant phase difference", "random phase", "different frequencies", "different polarizations"], answer: 0 },
    { question: "Q12: For a thin lens, sign convention matters: true or false?", options: ["True", "False", "only for mirrors", "only for thick lens"], answer: 0 },
    { question: "Q13: Group velocity can exceed c in anomalous dispersion but information velocity cannot: true or false?", options: ["True", "False", "both exceed", "neither"], answer: 0 },
    { question: "Q14: Thermodynamic free energy F = U - TS is minimized at equilibrium for constant T: true or false?", options: ["True", "False", "only adiabatic", "only isolated"], answer: 0 },
    { question: "Q15: If vector field has non-zero divergence, there exist sources or sinks: true or false?", options: ["True", "False", "only rotations", "only conservative"], answer: 0 },
    { question: "Q16: For a rigid body, center of mass motion decouples from rotational motion if external force acts at center: true or false?", options: ["True", "False", "only if torque zero", "only if massless"], answer: 0 },
    { question: "Q17: In time-independent perturbation theory, first-order energy correction is:", options: ["⟨ψ0|H'|ψ0⟩", "0", "second order", "⟨ψ1|H'|ψ1⟩"], answer: 0 },
    { question: "Q18: If electrical conductivity σ increases, resistivity ρ:", options: ["decreases", "increases", "unchanged", "zero"], answer: 0 },
    { question: "Q19: If two events are space-like separated they cannot be causally connected: true or false?", options: ["True", "False", "depends", "sometimes"], answer: 0 },
    { question: "Q20: The Riemann zeta function ζ(s) analytic continuation has zeros on critical strip Re(s) between 0 and 1: Riemann hypothesis states nontrivial zeros lie on:", options: ["Re(s)=1/2", "Re(s)=1", "Re(s)=0", "varies"], answer: 0 }
  ]
];

// ---------- State ----------
let currentLevel = 0;
let totalScore = 0;
let levelScore = 0;
let answered = [];
let timer = null;
let timeLeft = TIME_PER_LEVEL_SECONDS;

// ---------- Helpers ----------
function ensureLevelLength(levelQuestions) {
  // Ensures exactly QUESTIONS_PER_LEVEL items (but our bank already does)
  const copy = Array.isArray(levelQuestions) ? [...levelQuestions] : [];
  if (copy.length > QUESTIONS_PER_LEVEL) copy.length = QUESTIONS_PER_LEVEL;
  while (copy.length < QUESTIONS_PER_LEVEL) {
    copy.push({
      question: `Placeholder Q${copy.length+1}`,
      options: ["A","B","C","D"],
      answer: -1
    });
  }
  return copy;
}

// ---------- Rendering ----------
function startLevel(levelIndex) {
  levelScore = 0;
  answered = new Array(QUESTIONS_PER_LEVEL).fill(false);

  document.getElementById("next-level-btn").classList.add("hidden");
  document.getElementById("submit-btn").classList.remove("hidden");
  document.getElementById("level-info").innerText = `Level ${levelIndex + 1} of ${allLevels.length}`;

  const rawQuestions = allLevels[levelIndex] || [];
  const questions = ensureLevelLength(rawQuestions);

  const oldSummary = document.getElementById("level-summary");
  if (oldSummary) oldSummary.remove();

  renderQuestions(questions);
  resetTimer();
}

function renderQuestions(questions) {
  const container = document.getElementById("questions");
  container.innerHTML = "";

  questions.forEach((q, idx) => {
    const qBlock = document.createElement("div");
    qBlock.className = "question-block";

    const qText = document.createElement("div");
    qText.className = "question";
    qText.innerHTML = `<strong>Q${idx + 1}.</strong> ${q.question || ""}`;

    const optionsDiv = document.createElement("div");
    optionsDiv.className = "options";

    // Render exactly 4 option buttons
    for (let optIndex = 0; optIndex < 4; optIndex++) {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "option";
      btn.innerText = (q.options && q.options[optIndex]) ? q.options[optIndex] : `Option ${optIndex+1}`;
      btn.addEventListener("click", () => handleAnswer(idx, optIndex, q.answer, optionsDiv));
      optionsDiv.appendChild(btn);
    }

    qBlock.appendChild(qText);
    qBlock.appendChild(optionsDiv);
    container.appendChild(qBlock);
  });
}

// ---------- Answering & Scoring ----------
function handleAnswer(questionIndex, chosenIndex, correctIndex, optionsDiv) {
  if (answered[questionIndex]) return; // prevent double scoring
  answered[questionIndex] = true;

  const opts = optionsDiv.querySelectorAll(".option");
  opts.forEach((optEl, i) => {
    optEl.disabled = true;
    if (i === correctIndex && correctIndex >= 0) optEl.classList.add("correct");
  });

  if (chosenIndex !== correctIndex) {
    const chosenEl = opts[chosenIndex];
    if (chosenEl) chosenEl.classList.add("incorrect");
  } else if (chosenIndex === correctIndex && correctIndex >= 0) {
    levelScore++;
    totalScore++;
  }
}

// ---------- Timer ----------
function resetTimer() {
  clearInterval(timer);
  timeLeft = TIME_PER_LEVEL_SECONDS;
  updateTimerDisplay();
  timer = setInterval(() => {
    timeLeft--;
    updateTimerDisplay();
    if (timeLeft <= 0) {
      clearInterval(timer);
      finishLevel();
    }
  }, 1000);
}

function updateTimerDisplay() {
  const mm = String(Math.floor(timeLeft / 60)).padStart(2, "0");
  const ss = String(timeLeft % 60).padStart(2, "0");
  const timerEl = document.getElementById("timer");
  if (timerEl) timerEl.innerText = `Time Left: ${mm}:${ss}`;
}

// ---------- Finish / Submit ----------
function revealAllCorrectAnswers() {
  const questionBlocks = document.querySelectorAll("#questions .question-block");
  const levelQuestions = ensureLevelLength(allLevels[currentLevel]);
  questionBlocks.forEach((block, qIndex) => {
    const correctIndex = levelQuestions[qIndex] ? levelQuestions[qIndex].answer : -1;
    const opts = block.querySelectorAll(".option");
    opts.forEach((optEl, i) => {
      optEl.disabled = true;
      if (i === correctIndex && correctIndex >= 0) optEl.classList.add("correct");
    });
  });
}

function finishLevel() {
  clearInterval(timer);
  document.getElementById("submit-btn").classList.add("hidden");
  revealAllCorrectAnswers();

  const container = document.getElementById("questions");
  const levelTotal = QUESTIONS_PER_LEVEL;
  const cumulativeTotalPossible = (currentLevel + 1) * QUESTIONS_PER_LEVEL;

  // Remove existing summary if present
  const oldSummary = document.getElementById("level-summary");
  if (oldSummary) oldSummary.remove();

  // Insert summary after questions container
  container.insertAdjacentHTML(
    "afterend",
    `<div id="level-summary" style="text-align:center;margin-top:12px;color:#fff">
      <h2>Level Completed!</h2>
      <p style="font-weight:700">This level: ${levelScore} / ${levelTotal}</p>
      <p>Total so far: ${totalScore} / ${cumulativeTotalPossible}</p>
    </div>`
  );

  if (currentLevel < allLevels.length - 1) {
    document.getElementById("next-level-btn").classList.remove("hidden");
  } else {
    // all done
    const endMsg = document.createElement("div");
    endMsg.style.textAlign = "center";
    endMsg.style.marginTop = "10px";
    endMsg.innerHTML = `<h2>All Levels Completed!</h2><p style="font-weight:700">Final Score: ${totalScore} / ${cumulativeTotalPossible}</p>`;
    document.getElementById("quiz-container").appendChild(endMsg);
  }
}

function submitNow() {
  const oldSummary = document.getElementById("level-summary");
  if (oldSummary) oldSummary.remove();
  finishLevel();
}

// ---------- UI Hooks ----------
document.getElementById("submit-btn").addEventListener("click", () => {
  submitNow();
});

document.getElementById("next-level-btn").addEventListener("click", () => {
  const oldSummary = document.getElementById("level-summary");
  if (oldSummary) oldSummary.remove();

  currentLevel++;
  if (currentLevel >= allLevels.length) currentLevel = allLevels.length - 1;
  startLevel(currentLevel);
});

// ---------- Start ----------
startLevel(currentLevel);

