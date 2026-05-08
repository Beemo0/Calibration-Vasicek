# Vasicek Model Calibration and Simulation
 
## Overview
 
This project implements the **Vasicek short rate model** calibration using both **yield curve market data** and **historical interest rate data**. The main objectives are to calibrate the model parameters to observed market yield curves, simulate interest rate trajectories, and study the shape of the term structure under different initial conditions.
 
The implementation is written in **Python** and relies on the **Levenberg-Marquardt nonlinear least squares algorithm** for parameter calibration, as well as standard stochastic simulation techniques.
 
---
 
## Project Structure
 
```
TP3/
├── TP3.py                        # Python implementation of all computations
└── README.md                     # Project documentation
```
 
---
 
# Theoretical Background
 
## Vasicek Model
 
The Vasicek model is a **short rate model** that describes the dynamics of the instantaneous interest rate $r_t$ via the stochastic differential equation:
 
$$dr_t = (\eta - \gamma r_t)dt + \sigma dW_t$$
 
where $W_t$ is a standard Brownian motion.
 
| Parameter | Meaning |
|-----------|---------|
| $\eta$ | Long-term mean level (times $\gamma$) |
| $\gamma$ | Speed of mean reversion |
| $\sigma$ | Volatility of the short rate |
| $r_0$ | Initial short rate |
 
The **exact solution** is:
 
$$r_t = r_0 e^{-\gamma t} + \frac{\eta}{\gamma}(1 - e^{-\gamma t}) + \sigma \int_0^t e^{-\gamma(t-s)} dW_s$$
 
with mean $\mathbb{E}[r_t] = r_0 e^{-\gamma t} + \frac{\eta}{\gamma}(1 - e^{-\gamma t})$ and variance $\text{Var}[r_t] = \frac{\sigma^2}{2\gamma}(1 - e^{-2\gamma t})$.
 
---
 
## Zero-Coupon Bond Pricing
 
The price of a zero-coupon bond under the Vasicek model has the closed-form solution:
 
$$P(r, t; T) = e^{A(t,T) - r \cdot B(t,T)}$$
 
with:
 
$$B(t, T) = \frac{1 - e^{-\gamma \tau}}{\gamma}, \quad \tau = T - t$$
 
$$A(t, T) = (B(t,T) - \tau) \frac{\eta\gamma - \sigma^2/2}{\gamma^2} - \frac{\sigma^2 B^2(t,T)}{4\gamma}$$
 
---
 
## Yield Curve
 
The **yield curve** (or spot rate) is derived from bond prices:
 
$$Y(t, T) = -\frac{\log P(r, t; T)}{T - t} = -\frac{A(0,T) - r_0 \cdot B(0,T)}{T}$$
 
The Vasicek model produces three possible curve shapes depending on $r_0$:
- **Upward sloping** (typical): $r_0 < \frac{\eta}{\gamma}$
- **Humped**: intermediate $r_0$
- **Inverted**: $r_0 > \frac{\eta}{\gamma}$
Limiting behaviors:
 
$$\lim_{T \to 0} Y(0, T) = r_0, \qquad \lim_{T \to \infty} Y(0, T) = \frac{\eta}{\gamma} - \frac{1}{2}\left(\frac{\sigma}{\gamma}\right)^2$$
 
---
 
# Implemented Functions
 
### Bond Pricing
- `bond_price(r0, t, T, eta, gamma, sigma)` — Computes $P(r, t; T)$ using the closed-form Vasicek formula.
---
 
### Yield Curve
- `yield_curve(r0, T, eta, gamma, sigma)` — Computes $Y(0, T)$ for a given set of parameters.
---
 
### Jacobian
- `jacobian(r0, T, eta, sigma2, gamma)` — Computes the $10 \times 3$ Jacobian matrix used in the Levenberg-Marquardt algorithm:
$$J_{pj} = \frac{\partial(\text{Res})_p}{\partial \beta_j}$$
 
---
 
### Levenberg-Marquardt Calibration
- `levenberg_marquardt(Y_market, T, r0)` — Minimizes the sum of squared residuals:
$$\min_{\beta} \Phi = \sum_{p=1}^{10} \left(Y^{\text{market}}_p - Y(\beta, T_p)\right)^2$$
 
---
 
### Simulation
- `simulate_vasicek(r0, T, N, eta, gamma, sigma)` — Simulates a path of the short rate using the discretization:
$$r_{i+1} = r_i e^{-\gamma \Delta t} + \frac{\eta}{\gamma}(1 - e^{-\gamma \Delta t}) + \sigma \sqrt{\frac{1 - e^{-2\gamma \Delta t}}{2\gamma}} \cdot \mathcal{N}(0,1)$$
 
---
 
# Part I — Yield Curve Shape Analysis
 
The yield curve $T \mapsto Y(0, T)$ is plotted for fixed parameters:
 
$$\gamma = 0.25, \quad \eta = 0.25 \times 0.03, \quad \sigma = 0.02$$
 
Three shapes are illustrated by varying $r_0$:
 
| Case | $r_0$ | Shape |
|------|--------|-------|
| a | 0.01 | Upward sloping |
| b | 0.027 | Slightly humped |
| c | 0.05 | Inverted |
 
---
 
# Part II — Calibration to Market Yield Curve
 
## Market Data (t = 0)
 
| Maturity $T_i$ (months) | 3 | 6 | 9 | 12 | 15 | 18 | 21 | 24 | 27 | 30 |
|--------------------------|-------|-------|--------|-------|--------|--------|--------|--------|------|--------|
| $Y^{\text{market}}_i$ | 0.035 | 0.041 | 0.0439 | 0.046 | 0.0484 | 0.0494 | 0.0507 | 0.0514 | 0.052 | 0.0523 |
 
## Algorithm: Levenberg-Marquardt
 
**Parameters to calibrate:**
 
$$\beta_1 \equiv \eta, \quad \beta_2 \equiv \sigma^2, \quad \beta_3 \equiv \gamma$$
 
**Initialization:**
 
$$\eta = 1, \quad \sigma^2 = 1, \quad \gamma = 1, \quad r_0 = 0.023, \quad \varepsilon = 10^{-9}, \quad \lambda = 0.01$$
 
**Update rule:**
 
$$d_k = -(J^T J + \lambda I)^{-1} \cdot J^T \cdot \text{Res}$$
$$\beta_{k+1} = \beta_k + d_k$$
 
The algorithm converges when $||\beta_{k+1} - \beta_k|| < \varepsilon$.
 
---
 
## Recalibration at t = 1
 
The yield curve is recalibrated one year later with $r_1 = 0.04$ and new market data:
 
| Maturity $T_i$ (months) | 3 | 6 | 9 | 12 | 15 | 18 | 21 | 24 | 27 | 30 |
|--------------------------|-------|-------|------|------|------|------|------|------|--------|------|
| $Y^{\text{market}}_i$ | 0.056 | 0.064 | 0.074 | 0.081 | 0.082 | 0.09 | 0.087 | 0.092 | 0.0895 | 0.091 |
 
Note: the time to maturity is now $T - t = T - 1$.
 
---
 
# Part III — Calibration to Historical Data
 
## Simulation
 
A path of the Vasicek process is simulated with:
 
$$T = 5, \quad \eta = 0.6, \quad \gamma = 4, \quad \sigma = 0.08$$
 
The pairs $(x_i, y_i) = (r_i, r_{i+1})$ are plotted and fitted to the linear model $y = ax + b$.
 
## Parameter Recovery
 
From the regression coefficients $a$, $b$ and residual variance $D^2$, the Vasicek parameters are recovered as:
 
$$\gamma = -\frac{\ln a}{\Delta t}, \qquad \eta = \gamma \frac{b}{1 - a}, \qquad \sigma = D \cdot \sqrt{\frac{-2 \ln a}{\Delta t (1 - a^2)}}$$
 
The recovered values are compared to the original simulation parameters to validate the calibration procedure.
 
---
 
# Libraries Used
 
```
numpy
matplotlib
pandas
scipy
scikit-learn
```
 
Install dependencies with:
 
```bash
pip install numpy matplotlib pandas scipy scikit-learn
```
 
---
 
# Usage
 
Run the main script:
 
```bash
python TP3.py
```
 
The script will:
1. Plot yield curve shapes (Part I)
2. Calibrate model parameters to market yield data at $t = 0$ (Part II)
3. Recalibrate at $t = 1$ with updated market data
4. Calibrate to LIBOR and Indonesian Government Securities data
5. Simulate historical paths and recover parameters (Part IV)
---
 
# Key Concepts Covered
 
- Vasicek short rate model
- Zero-coupon bond pricing (closed-form)
- Term structure of interest rates
- Yield curve shapes (normal, humped, inverted)
- Levenberg-Marquardt nonlinear least squares
- Jacobian computation for gradient-based optimization
- Stochastic simulation of interest rate paths
- Historical calibration via linear regression
---
 
# Disclaimer
 
This repository contains the **code and methodology** used in an academic laboratory session at CY Tech (Université Paris-Cergy). The original written report has been removed due to **academic and data usage restrictions**.
