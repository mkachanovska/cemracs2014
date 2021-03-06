\documentclass[proc]{edpsmath}


\begin{document}
\section{Discretization}
\subsection{Discretization of the Frequency Domain Problem}

We look for a ratio $h(\nu)$ that would ensure that 
\begin{align*}
 \|E^{\nu}_{1}-E^{\nu,h}_{1}\|_{L_{2}(\Omega)}<\epsilon,
\end{align*}
for a fixed $\epsilon>0$. 

To conjecture the dependence of $h(\nu)$, we suggest the following (not rigorous argument). 
Recall that for the problem (\ref{}) the C\'ea's lemma gives
\begin{align*}
 \|\mathbf{E}^{\nu}-\mathbf{E}^{h,\nu}\|\leq \frac{C_c}{C_i}\min_{\mathbf{v}\in V_h}\|\mathbf{E}-\mathbf{v}\|_{V},
\end{align*}
where $C_c$ is the continuity and $C_i$ is the coercivity constants.
\begin{lemma}{}
Let $\alpha(x), \; \delta(x)$ satisfy
\begin{align*}
\end{align*}
Then, for all $\nu_*>0$ and for all $\nu>\nu_*$, the error between the numerical solution $\mathbf{E}^{h,\nu}$ and the 
solution $\mathbf{E}^{\nu}$
\begin{align*}
\|\mathbf{E}^{\nu}-\mathbf{E}^{h,\nu}\|_{V}\leq 
\end{align*}
\end{lemma}
\begin{proof}
The C\'ea's lemma states 
\begin{align*}
 \|\mathbf{E}^{\nu}-\mathbf{E}^{h,\nu}\|\leq \frac{C_c}{C_i}\min_{\mathbf{v}\in V_h}\|\mathbf{E}-\mathbf{v}\|_{V},
\end{align*}
where $C_c$ is the continuity and $C_i$ is the coercivity constants.
Indeed,
\begin{align*}
 \min_{\mathbf{v}\in V_h}\|\mathbf{E}-\mathbf{v}\|_{V}\leq \min_{\mathbf{v}\in V_h}\|\mathbf{E}-I^{h}\mathbf{v}\|,
\end{align*}
where $I^{h}\mathbf{v}$ is an interpolation operator. The estimates from \cite[Chapter 0]{Brenner_Scott} give us
\begin{align*}
 \|E_1-I^{h}E_1\|_{L_{2}(\Omega)}\leq Ch^2|E_1''|_{L_{2}(\Omega)},\\
 \|E_2-I^{h}E_2\|_{H_{1}(\Omega)}\leq Ch|E_2''|_{H_{1}(\Omega)}.
\end{align*}
To obtain these bounds, we require an explicit bound on $E_1, \; E_2$. In \cite{Despers et al Hybrid} 
it was shown that 
\begin{align*}
 
\end{align*}


\end{proof}

\subsection{Discretization of the Time Domain Problem}
\section{Numerical Experiments}
\subsection{Frequency Domain}
\subsection{Time Domain}
\end{document}