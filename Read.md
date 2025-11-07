# Optimal Transport-Based Variable Recoding

This repository contains the code and simulation framework used in the paper:

> **Optimal Transport-Based Variable Recoding: A Comparative Study with MICE and kNN Approaches**  

---

## 📘 Overview

This project investigates the problem of **variable recoding** in heterogeneous observational datasets where categorical or ordinal outcomes are defined differently across sources.  
We compare a **Optimal Transport (OT)** method with two standard machine learning techniques:

- **MICE** (Multiple Imputation by Chained Equations)  
- **kNN** (k-Nearest Neighbours Imputation)

The study evaluates how each method performs under varying sample sizes, correlations, and missingness mechanisms.  
Performance is assessed through **precision scores**  between reconstructed and reference distributions.

---
We also  consider a real dataset, the National Child Development Styudy(NCDS) and comapre distributions of results between the methods using the **Wasserstein distances**. 
## ⚙️ Simulation Design

The simulations are structured into **10 scenarios**, summarized below:

| Case | Scenario Type | Missingness | Imputation | Coefficient Equality |
|------|----------------|--------------|-------------|----------------------|
| 1–4  | Complete data (varying sample size, imbalance, R²) | No | No | Yes/No |
| 5–7  | Incomplete data (MCAR, MAR, MNAR) | Yes | No | No |
| 8–10 | Imputed data (MCAR, MAR, MNAR) | Yes | Yes | No |

Each scenario tests the robustness of OT, MICE, and kNN in reconstructing consistent variable codings under different data conditions.

---

## 🧩 Repository Structure


