::: {.callout-warning}
All text in this post is AI generated until otherwise indicated.
:::


# A Trajectory-Based Framework for Long-Term Holdout Evaluation

## 1. Notation and Exposure Space

To evaluate the cumulative and long-term impact of multiple distinct treatments (innovations) introduced sequentially over time, we define a multi-dimensional, longitudinal exposure space. Let $k$ denote the total number of distinct treatments under consideration. We represent the exposure state of subject $i$ at a discrete time block $t$ as a $k$-dimensional binary vector:

$$\mathbf{a}_{t} = [a_{t,1}, a_{t,2}, \dots, a_{t,k}]^\top \in \{0, 1\}^k$$

where $a_{t,j} = 1$ if subject $i$ is exposed to treatment $j$ at time $t$, and $0$ otherwise. The special cases $\mathbf{1}$ and $\mathbf{0}$ denote the all-ones and all-zeros vectors, respectively.

Because treatments roll out progressively across historical cycles, we track a subject's longitudinal path as an exposure trajectory. Let $\mathbf{\bar a}_t = \{\mathbf{a}_\tau\}_{\tau=1}^t$ be the sequence of exposure vectors from the initiation of the program $t=1$ up to time $t$. The complete trajectory over an experimental rollout horizon $T$ is represented by the $T \times k$ matrix object $\mathbf{\bar a} = \{\mathbf{a}_t\}_{t=1}^T$. We write $\mathbf{\bar 1}$ for the counterfactual trajectory in which every treatment is active at all times, and $\mathbf{\bar 0}$ for the trajectory in which no treatment is ever active.

Following the potential outcomes framework extended to longitudinal treatments, let $Y_i(\mathbf{\bar a})$ denote the potential outcome for subject $i$ under the complete exposure trajectory $\mathbf{\bar a}$. Crucially, $Y_i(\mathbf{\bar a}) \in \mathbb{R}$ is a scalar-valued metric mapped from the matrix-valued trajectory argument, representing a cumulative or post-rollout behavioral summary.

---

## 2. Structural Assumptions

To bound the data-generating process and isolate the causal estimands of interest, we formalize the following core scoping and structural assumptions.

### Assumption 1: No Interference (SUTVA-Style Isolation)
The potential outcome for any subject $i$ under a joint treatment assignment vector across the entire population depends strictly on their own assigned trajectory:

$$Y_i(\mathbf{\bar a}_1, \mathbf{\bar a}_2, \dots, \mathbf{\bar a}_N) = Y_i(\mathbf{\bar a}_i) \quad \forall i \in \{1, \dots, N\}$$

### Assumption 2: Path Dependence
Let $\mathbf{\bar a}^\alpha$ and $\mathbf{\bar a}^\beta$ be two distinct exposure trajectories over the rollout window $1 \dots T$. If they terminate in the identical cross-sectional exposure configuration at time $T$ ($\mathbf{a}_T^\alpha = \mathbf{a}_T^\beta$), the potential outcomes are not assumed to be equal:

$$\mathbb{E}[Y_i(\mathbf{\bar a}^\alpha)] \neq \mathbb{E}[Y_i(\mathbf{\bar a}^\beta)]$$

This non-Markovian behavior is driven by two distinct structural mechanisms:
1. *Behavioral Memory and Hysteresis:* The underlying psychological or operational outcome process possesses memory. Phenomena such as user adaptation, learning curves, novelty decay, or UI fatigue imply that the effect of an exposure depends heavily on its duration and historical onset.
2. *Time-Aggregated Outcomes:* When the macro-metric $Y_i$ is computed as a sum or average over the experimental timeline, the trajectory matters mechanically. Even under a memoryless, instantaneous treatment effect, an exposure that was active longer contributes greater mass to the aggregate scalar than a late-stage deployment.

---

## 3. Cohort Dynamics and the Steady-State Horizon

In practice, a population is randomized at $t=1$ into three operational arms: the *Status Quo* ($SQ$), the *Winners Only* ($WO$), and the *Active Testing Pool* ($R$). Over the sequential rollout period $t = 1, \dots, T$, winning treatments selected from group $R$ are cascaded to the $WO$ cohort. This enforces a strict structural monotonicity on the $WO$ trajectory as coordinates of $\mathbf{a}_t$ shift irreversibly from $0$ to $1$. 

Eventually, the rollout phase terminates at time $T$, at which point all $k$ successful treatments are permanently activated for the $WO$ cohort, yielding the terminal vector $\mathbf{a}_T = \mathbf{1}$. To evaluate the true cumulative impact while accounting for the lag effects implied by Assumption 2, we introduce a post-rollout evaluation window spanning from $T+1$ to $T+K$. 

During this steady-state horizon, exposures are frozen:
$$\mathbf{a}_\tau^{WO} = \mathbf{1} \quad \text{and} \quad \mathbf{a}_\tau^{SQ} = \mathbf{0} \quad \forall \tau \in \{T+1, \dots, T+K\}$$

Let $\mathbf{\bar a}_{T+K}^{WO}$ and $\mathbf{\bar a}_{T+K}^{SQ}$ represent these full, extended trajectories. We introduce a final stabilizing assumption:

### Assumption 3: Equilibration (Absence of Distant Carryover)
For a sufficiently large evaluation buffer $K$, the marginal memory of the historical rollout path dissipates, allowing the potential outcome process to reach a steady-state equilibrium:

$$\mathbb{E}[Y_{i, T+K}(\mathbf{\bar a}_{T+K}^{WO})] \approx \mathbb{E}[Y_{i, T+K}(\mathbf{\bar 1}_{T+K})] \quad \text{and} \quad \mathbb{E}[Y_{i, T+K}(\mathbf{\bar a}_{T+K}^{SQ})] = \mathbb{E}[Y_{i, T+K}(\mathbf{\bar 0}_{T+K})]$$

---

## 4. Causal Estimands and Identification

Under the steady-state horizon, the principal objective of the holdout framework is to recover the structural contrast between the fully optimized product and the absolute baseline. The terminal causal estimand evaluated over the window $T+1$ to $T+K$ is defined as:

$$\tau_{\text{prog}}(T+K) = \mathbb{E}[Y_{i, T+K}(\mathbf{\bar a}_{T+K}^{WO})] - \mathbb{E}[Y_{i, T+K}(\mathbf{\bar a}_{T+K}^{SQ})]$$

By virtue of Assumption 3, the trajectory-bound historical artifact yields to a clean cross-sectional proxy:

$$\tau_{\text{prog}}(T+K) \approx \mathbb{E}[Y_{i, T+K}(\mathbf{\bar 1})] - \mathbb{E}[Y_{i, T+K}(\mathbf{\bar 0})]$$

This estimand effectively maps to the classical average treatment effect (ATE) of the fully unified treatment suite against the status quo, isolated from the temporal noise of the deployment schedule. 

Given initial randomization to the $WO$ and $SQ$ paths at $t=1$, the trajectories are orthogonal to individual potential outcomes ($\mathbf{\bar a} \perp \{Y_i(\mathbf{\bar a}^{WO}), Y_i(\mathbf{\bar 0})\}$). Consequently, the structural estimand is non-parametrically identified directly via the contrast of the observed empirical cohort means during the steady-state window:

$$\hat{\tau}_{\text{prog}}(T+K) = \frac{1}{N_{WO}} \sum_{i \in WO} Y_{i, T+K} - \frac{1}{N_{SQ}} \sum_{i \in SQ} Y_{i, T+K}$$