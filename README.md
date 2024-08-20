# DynaVac: Modeling the Evolution of Immunogenicity Dynamics Provides Insight into SARS-CoV-2 Immunization Strategies

## Project Description

DynaVac is an innovative computational model developed to quantitatively describe and predict the immune response dynamics induced by SARS-CoV-2 vaccination. The model uses ordinary differential equations (ODEs) based on advanced immunological mechanisms to abstract key cellular and molecular processes of antibody immune response, including mRNA translation to produce antigens, naive B-cell affinity maturation, differentiation of memory B-cells and plasma cells, as well as antigen-antibody neutralization. The model also accounts for the activation, proliferation, competition, and differentiation of memory B-cells generated during primary immunization when reactivated in subsequent immunizations. Specifically, the equations describe the dynamic changes of key components such as antigens ($Ag$), naïve B-cells ($N$), memory B-cells ($Mon, Moff$), and antibodies ($Ab$). The general form of the equations is as follows:
$$
\begin{gather}
\frac{dAg_i}{dt} = -\sum_{j}^{n} \gamma_{\text{neu}_i} c_{i,j}  Ag_i  Ab_j - \gamma_{\text{Ag}}  Ag_i \\
\frac{dN_i}{dt} = a_N(t) s_N \frac{Ag_i}{Ag_i + K} N_i (1-N_i) - d_N N_i \\
\frac{dM_{\text{off}_i}}{dt} = k_{N2M} N_i \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{off}_i} \\
\frac{dM_{\text{on}_i}}{dt} = a_M(t) s_M \left(\sum_{k}^{n} \frac{Ag_k}{Ag_k + K} \frac{c_{i,k} M_{\text{on}_i}}{c_{i,k} M_{\text{on}_i} + m_0} \right) M_{\text{on}_i} \left[ 1 - \sum_{j}^{n} (M_{\text{off}_j} + M_{\text{on}_j}) \right] - d_M M_{\text{on}_i} \\
\frac{dAb_i}{dt} = a_N(t) p_N N_i + a_M(t) \frac{p_M M_{\text{on}_i} c_{i,k}}{c_{i,k} + c_0} - \gamma_{\text{neu}_i} \sum_{j}^{n} c_{j,i} Ag_j Ab_i - \gamma_{\text{Ab}} Ab_i
\end{gather}
$$
