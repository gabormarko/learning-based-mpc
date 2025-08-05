# Learning-based Rigid Tube Model Predictive Control
Code implementation for the Learning-based Rigid Tube Model Predictive
Control project as part of the Practical Course Motion planning for the autonomous vehicle EDGAR SS2025.

## Abstract
In this report we address model predictive control (MPC) for discrete-time linear systems subject to bounded additive disturbances, where the true disturbance set is unknown. Building on the learning-based rigid tube MPC framework, which integrates conservative prior knowledge with online disturbance learning, we extend its implementation to support higher-dimensional systems. The disturbance set is parametrized as a homothetic transformation of a conservative prior, enabling efficient construction and update of the rigid tube that bounds disturbance propagation.

We modularize all core components of the learning-based rigid tube MPC framework to support higher-order systems. To evaluate performance in higher dimensions, we conduct a scalability analysis of the core components, including disturbance learning, set calculations, tube updates, and optimization problems. We analyze the scalability of the proposed method through complexity calculations and simulations on 3D and 4D systems. Comparative results with baseline rigid tube MPC and scenario MPC demonstrate that our generalized approach maintains reduced conservatism, robust constraint satisfaction, and computational feasibility for real-time control in higher-order systems.

## General information
This repository is an extension for the codebase for the following article:
```
@inproceedings{gao2024learning,
  title={Learning-based Rigid Tube Model Predictive Control},
  author={Yulong Gao, Shuhao Yan, Jian Zhou, Mark Cannon, Alessandro Abate, and Karl H. Johansson},
  booktitle={Learning for Dynamics and Control Conference},
  year={2023},
  pages={}
} 
```

Although the theoretical framework presented in the article supports arbitrary state and input dimensions, the original codebase (https://github.com/Safe-Autonomy-and-Intelligence-Lab/learning-based-rigid-tube-rmpc) was limited to 2D systems. In this repository, we refactored the control pipeline into modular, simension-agnostic components to ensure general applicability.


## Packages for running the code
To run the code you need to install:

**CasADi**: https://web.casadi.org/;

**MPT**: https://www.mpt3.org/Main/Installation; (The installation will automatically install Yalmip, which is also necessary for running the code.)

## Introduction to the files

## offline_parameters_computation_nxn.m
Generalized version of the `offline_parameters_computation.m` file.
Calculate and define all parameters and save the results in **parameters_{nx}x{nx}.mat**. Since the parameters have been provided in the repository, if you do not want to update those values you can avoid running this file.

This file is modularized, meaning the dimensions of the state space representation can be set also to higher dimensions. Tested for nx={2,3,4}  dimensional systems.

## Functions_general
Generalized version of all of the functions in the original `Functions` folder.
### MRPISet.m:

Compute the minimum robust positive invariant set $\mathbb{S}$.

### ComputeFeasibleRegion.m:

Compute the feasible region corresponding to different disturbance sets, e.g., $\mathbb{W}$, $\mathbb{W}_{\rm true}$, and $\hat{\mathbb{W}}_k^{\star}$.

### InitialSetComputation.m:

Learn the initial uncertainty set $\hat{\mathbb{W}}^{\star}_0$ based on $\mathbb{W}$, using initial information set
$\mathcal{I}_0^w$.

### ModelingCar.m:

Modeling of the ego vehicle (EV) and the leading vehicle (LV).

### NominalRobustMPC.m:

Conventional Robust MPC controller.

### UQRobustMPC.m:

The proposed uncertainty quantification-based Robust MPC controller.

### ScenarioMPC.m:
The scenario MPC controller.

## Cases_general
Generalized version of all of the scripts in the original `Cases` folder.

### Case_1_Feasible_Region.m
Compute the feasible region of UQ-RMPC and nominal RMPC, results are saved in **Results/Results_1.mat**, and figures are produced by **Results/Fig_Case_1.m**.

### Case_2_MC_Different_Initial_InformationSet.m
Monte-Carlo simulation learning the set $\mathbb{W}_{\rm true}$ with different sizes of initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_2.mat**, and figures are produced by **Results/Fig_Case_2.m**.

### Case_3_Online_UQRMPC_Different_Initial_InformationSet.m
Online evaluation of UQ-RMPC with different initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_3_large.mat** ($|\mathcal{I}_0^w| = 20000$) and **Results/Results_3_small.mat** ($|\mathcal{I}_0^w| = 100$), and figures are produced by **Results/Fig_Case_3.m**.

### Case_4_Online_UQRMPC_Long_Simulation_Step.m
Monte-Carlo simulation of UQ-RMPC when the simulation time is long enough, results are saved in **Results/Results_4.mat**, and figures are produced by **Results/Fig_Case_4.m**.

### Case_5_Feasibility_Evaluation_UQRPC.m
Monte-Carlo simulation to evaluate the feasibility of UQ-RMPC with different initial information set $\mathcal{I}_0^w$, results are saved in **Results/Results_5.mat**, e.g., **Results_5_10.mat** indicates the results with $|\mathcal{I}_0^w| = 10$, and so on.

### Case_6_Compare_With_SCMPC.m
Comparing the computation time with Scenario MPC.

## Visualizing 3D results
To visualize 3D disturbance sets and feasible regions, run the `visualize_3d_results.m` script with the corresponding `parameters_3x3.mat` and `Results_1_3x3.mat` files.

## Some implementation details
(1) It is not necessary to run **offline_parameters_computation_nxn.m** if you do not want to update the parameter values, but you need to run **run_first.m** to add the path of folders.

(2) In the article, the horizon $\nu_k$ is updated according to Algorithm 1, but this will change the number of constraints of MPC. In our code, we implemented Algorithm 1 and found that $\nu_k$ is almost equal to $\nu_s$. Therefore, we use $\nu_s$ to replace $\nu_k$. In practical applications, we can make $\nu_k$ long enough such that the condition on $\nu_k$ will be satisfied. For example, a suggestion can be $\nu_k = 2\nu_s$.
