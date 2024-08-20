# DynaVac: Modeling the Evolution of Immunogenicity Dynamics Provides Insight into SARS-CoV-2 Immunization Strategies

## Project Description

DynaVac is an innovative computational model developed to quantitatively describe and predict the immune response dynamics induced by SARS-CoV-2 vaccination. The model uses ordinary differential equations (ODEs) based on advanced immunological mechanisms to abstract key cellular and molecular processes of antibody immune response, including mRNA translation to produce antigens, naive B-cell affinity maturation, differentiation of memory B cells and plasma cells, as well as antigen-antibody neutralization. The model also accounts for the activation, proliferation, competition, and differentiation of memory B cells generated during primary immunization when reactivated in subsequent immunizations. Specifically, the equations describe the dynamic changes of key components such as antigens ($Ag$), naïve B cells ($N$), memory B cells ($M_{on}, M_{off}$), and antibodies ($Ab$). The general form of the equations is as follows:
```math
\begin{gather}
\frac{dAg_i}{dt} = -\sum_{j}^{n} \gamma_{\text{neu}_i} c_{i,j}  Ag_i  Ab_j - \gamma_{\text{Ag}}  Ag_i \\
\frac{dN_i}{dt} = a_N(t) s_N \frac{Ag_i}{Ag_i + K} N_i (1-N_i) - d_N N_i \\
\frac{dM_{\text{off}_i}}{dt} = k_{N2M} N_i \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{off}_i} \\
\frac{dM_{\text{on}_i}}{dt} = a_M(t) s_M \left(\sum_{k}^{n} \frac{Ag_k}{Ag_k + K} \frac{c_{i,k} M_{\text{on}_i}}{\sum_{j}^{n}c_{j,k} M_{\text{on}_j} + m_0} \right) M_{\text{on}_i} \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{on}_i} \\
\frac{dAb_i}{dt} = a_N(t) p_N N_i + a_M(t) \frac{p_M M_{\text{on}_i} c_{i,k}}{c_{i,k} + c_0} - \gamma_{\text{neu}_i} \sum_{j}^{n} c_{j,i} Ag_j Ab_i - \gamma_{\text{Ab}} Ab_i
\end{gather}
```
The subscript $i$ represents the variant index. For detailed derivation of the equations and parameter meanings, please refer to the supplementary notes in the related publication [link to the paper](paper).

By fitting antibody-pseudovirus titration experimental data from various vaccination combinations, DynaVac can comprehensively and quantitatively characterize and predict the humoral immune response dynamics induced by SARS-CoV-2 vaccines during primary and booster immunizations. Our simulation results align with and expand upon the current understanding of immune imprinting. DynaVac provides a powerful quantitative tool for optimizing vaccine composition design and immunization strategies.

This repository includes the code to reproduce the main results of the paper (Figures 3 to 6) and offers methods for simulating personalized vaccination strategies. For quick online interactive simulations, visit the [DynaVac Online Interface](interface).

## Usage

### 1. Reproduction of Main Results

#### 1.1 Figure 3: Model Validation in Mice

To reproduce Figure 3 from the main text:

- Run `Fig3_validation_mice/plot_fig3.m`.
- The files `yd_parameters.mat`, `cyl_parameters.mat`, and `mix_parameters.mat` contain the parameter sets that were obtained by parameterizing the model using data from this study, the Yisimayi et al. dataset, and a mixed dataset, respectively.
- Run `plot_validation_this_study.m` or `plot_validation_Yisimayi.m` to visualize the model's performance on different datasets. Change the loaded parameter set at line 4 of the script as needed.

#### 1.2 Figure 4: Model Validation in Humans

To reproduce Figure 4:

- Run `Fig4_validation_human/plot_fig4.m`.
- The parameter set trained on human clinical data is stored in `Fig4_validation_human/human_parameter.mat`.

#### 1.3 Figure 5: Immune Imprinting vs. Antigenic Distance 

To reproduce Figure 5:

- Run `Fig5/plot_fig5.m`.
- To reproduce Figure 5f only, run `Fig5/plot_fig5f.m`.
- Simulation data required for Figure 5 is stored in `Fig5/fig5_data.mat`.
- To regenerate and store simulation data, run `Fig5/simu_fig5.m`. This simulation fixes the primary and booster vaccination schedules and primary vaccination dosage, and generates multiple simulation results by varying antigenic distance and booster dosage. Users can adjust fixed variables such as `vaccine_time` or `primary_dosage` to modify the simulation scenario.

#### 1.4 Figure 6a: Immune Imprinting vs. Primary Vaccination Dosage

To reproduce Figure 6a:

- Run `Fig6a/plot_fig6a.m`.
- Simulation data required for Figure 6a is stored in `Fig6a/fig6a_data.mat`.
- To regenerate and store simulation data, run `Fig6a/simu_fig6a.m`. This simulation fixes the primary and booster vaccination schedules and booster vaccination dosage, and generates multiple simulation results by varying antigenic distance and primary vaccination dosage. Users can adjust fixed variables such as `vaccine_time` or `booster_dosage` to modify the simulation scenario.

#### 1.5 Figure 6b: Immune Imprinting vs. Primary-Booster Interval

To reproduce Figure 6b:

- Run `Fig6b/plot_fig6b.m`.
- Simulation data required for Figure 6b is stored in `Fig6b/fig6b_data.mat`.
- To regenerate and store simulation data, run `Fig6b/simu_fig6b.m`. This simulation fixes the primary and booster vaccination dosage, and generates multiple simulation results by varying antigenic distance and the interval between primary and booster vaccinations. Users can adjust fixed variables such as `vaccine_amonut` to modify the simulation scenario.

### 2. Personalized Vaccination Strategy Simulation

`Personalized_vac/main.m` enables simulation of the molecular and cellular dynamics under personalized vaccination strategies. Consider a simulation involving $n$ variants and $m$ consecutive vaccinations. The vaccination strategy can be described by four $m$-dimensional vectors: $T,V,D,P$, where:

- $T = [t_{01},t_{02} ,..., t_{0m}]$ with $0 ≤ t_{01} < t_{02} < ... < t_{0m}$, representing the time of the `k-th` vaccination in days.
- $V = [v_1 ,v_2 ,..., v_m]$ with $vk ∈ {1,2,...,6}$, representing the variant used in the `k-th` vaccination (`1: WT`, `2: Alpha/Beta`, `3: Delta`, `4: BA.1`, `5: BA.2/BA.5`, `6: XBB.1.5`).
- $D = [d_1 ,d_2 ,... ,d_m]$ with $d_k > 0$, representing the dosage of the `k-th` vaccination in μg.
- $P = [p_1 ,p_2, ... ,p_m]$ with $p_k ∈ {1,2,3}$, representing the vaccine type (`1: protein vaccine`, `2: mRNA vaccine`, `3: inactivated vaccine`).

For multivalent vaccines, extend the monovalent vaccine simulation method by replacing the `m`-dimensional vector $V$ with a $q×m$ matrix:

```math

V = \begin{pmatrix}
v_{1,1} & v_{1,2} & \dots & v_{1,m} \\
v_{2,1} & v_{2,2} & \dots & v_{2,m} \\
\vdots & \vdots & \ddots & \vdots \\
v_{q,1} & v_{q,2} & \dots & v_{q,m}
\end{pmatrix}

```

Where $q$ represents the highest vaccine valency in the strategy, and $v_{jk} \in \{0,1,2,...,n\}$ represents the variant of the `j-th` antigenic component in the vaccine used in the `k-th` vaccination. If the valency of the vaccine used in the `k-th` vaccination is denoted as $q_k$ and is lower than the highest vaccine valency $q$, then $v_{j,k} = 0$ for all $j$ such that $q_k < j ≤ q$.

Running `main.m` generates area plots showing the dynamic changes in antigen, naïve B-cell maturity, memory B-cells, and antibody levels specific to each variant over time.

Users can customize the vaccination strategy for the 6 SARS-CoV-2 variants in `main.m` lines 13-18:

```matlab
% -----------personalized vaccination strategy---------
vaccine_variant = [1 4 3 2 5;
                   0 0 5 0 6]; % V: vaccination variant
% (1: WT, 2: Alpha/Beta, 3: Delta, 4: BA.1, 5: BA.2/BA.5, 6: XBB.1.5) 
vaccine_amount = [10 30 30 20 30]; % D: vaccination dosage
vaccine_time = [0 200 250 420 440]; % T: vaccination time (days)
vaccine_type = [2 2 2 2 2]; % P: vaccination type (1: protein, 2: mRNA, 3: inactivated)
```

(Ensure that the columns of the four input arrays are of the same length.)

The antigen-antibody cross-neutralization matrix for the six SARS-CoV-2 variants, derived from pseudovirus titration experiments in this study, does not require user customization. 
The legend for the variant-specific colors is found in `Personalized_vac/variant_legend.png`.
![legend](https://github.com/Jinkaichun/DynaVac/blob/main/Personalized_vac/variant_legend.png)
## License

This project is licensed under the XXX License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work was supported by [funding sources].

## Contact

For questions or issues, please contact jinkaichun@stu.pku.edu.cn.

