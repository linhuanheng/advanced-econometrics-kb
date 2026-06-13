---
title: 线性非观测效应面板数据模型
type: framework
domain: econometrics
source: Chapter V — Basic Linear Unobserved Effects Panel Data Models
created: 2026-05-19
updated: 2026-05-24
---

# Ch5: 线性非观测效应面板数据模型

## 核心问题

面板数据的核心动机：**解决遗漏变量问题**。当存在不可观测的个体异质性 $c_i$（如能力、动机、管理质量），且 $c_i$ 与解释变量相关时，单次横截面OLS产生有偏不一致的估计。面板数据通过追踪相同个体在不同时间点上的观测，允许我们**消除** $c_i$。

## 1. 基本框架

### 1.1 非观测效应模型 (Unobserved Effects Model, UEM)

对于随机抽取的横截面个体 $i$，在时期 $t$：

$$y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + c_i + u_{it}, \quad t = 1,2,...,T$$

| 符号 | 含义 | 关键性质 |
|------|------|---------|
| $y_{it}$ | 个体 $i$ 在时期 $t$ 的结果变量 | 可观测 |
| $\mathbf{x}_{it}$ | $1 \times K$ 时变解释变量向量 | 可观测 |
| $\boldsymbol{\beta}$ | $K \times 1$ 参数向量（我们关心的） | 待估计 |
| $c_i$ | **非观测效应** (unobserved effect) | 时间不变；代表个体 $i$ 的所有不可观测的持久特征 |
| $u_{it}$ | **特异性误差** (idiosyncratic error) | 随 $i$ 和 $t$ 都变化 |

**$c_i$ 的别名**: unobserved component, latent variable, unobserved heterogeneity, individual effect (当 $i$ 索引个人时), firm effect (当 $i$ 索引公司时)。

**平衡面板 (Balanced Panel)**: 每个横截面单位都有相同的时间段 $t=1,...,T$。$T$ 固定，$N \to \infty$。

### 1.2 三类解释变量

| 类型 | 变动维度 | 例子 | FE可否估计 |
|------|---------|------|-----------|
| 仅随 $i$ 变 | 跨 $i$，不跨 $t$ | 性别、种族、是否临河 | **否** |
| 仅随 $t$ 变 | 跨 $t$，不跨 $i$ | 年度宏观变量、时间虚拟变量 | 需与时间虚拟变量合并 |
| 随 $i$ 和 $t$ 都变 | 跨 $i$ 和 $t$ | 收入、工作经验、R&D支出 | **是** |

### 1.3 随机效应 vs 固定效应

| | 随机效应 (Random Effects) | 固定效应 (Fixed Effects) |
|------|------|------|
| $c_i$ 的性质 | 随机变量 | 待估参数（或任意相关） |
| 核心假设 | $\text{Cov}(\mathbf{x}_{it}, c_i) = 0$ | 允许任意 $\text{Cov}(\mathbf{x}_{it}, c_i)$ |
| 对 $c_i$ 的处理 | 放入复合误差项 | 通过变换**消除** |
| 时间不变变量 | 可以包含 | 无法包含 |

---

## 2. 严格外生性假设

### 2.1 条件期望形式 (最揭示性的表述)

$$E(y_{it} | \mathbf{x}_{i1}, \mathbf{x}_{i2}, ..., \mathbf{x}_{iT}, c_i) = E(y_{it} | \mathbf{x}_{it}, c_i) = \mathbf{x}_{it}\boldsymbol{\beta} + c_i \tag{10.12}$$

**含义**: 一旦控制了 $\mathbf{x}_{it}$ 和 $c_i$，其他时期的 $\mathbf{x}_{is}$ ($s \neq t$) 对 $y_{it}$ 无偏效应。

### 2.2 特异性误差形式

$$E(u_{it} | \mathbf{x}_{i1}, ..., \mathbf{x}_{iT}, c_i) = 0, \quad t = 1,2,...,T$$

**强度**: $E(\mathbf{x}_{is}' u_{it}) = 0$ 对所有 $s, t$ 成立——远强于仅要求同期不相关 $E(\mathbf{x}_{it}' u_{it}) = 0$。

### 2.3 为什么需要严格外生性？

严格外生性确保组内变换或差分变换后的OLS正交条件成立。若仅假设同期外生 $E(\mathbf{x}_{it}' u_{it}) = 0$，则对于 $s \neq t$：
- $E[(\mathbf{x}_{it} - \bar{\mathbf{x}}_i)' u_{it}]$ 中的 $\bar{\mathbf{x}}_i$ 与 $u_{it}$ 可能相关 → FE不一致
- $E[(\mathbf{x}_{it} - \mathbf{x}_{i,t-1})' (u_{it} - u_{i,t-1})]$ 中的 $\mathbf{x}_{i,t-1}$ 与 $u_{it}$ 可能相关 → FD不一致

### 2.4 严格外生性合理的例子

**大豆产出模型**: $y_{it}$ = 大豆产出, $\mathbf{x}_{it}$ = (资本, 劳动, 材料, 降雨量), $c_i$ = (土地质量, 管理能力)

→ 一旦控制当期投入和 $c_i$，其他年份的投入不影响当期产出：**合理**。

### 2.5 严格外生性不合理的例子

**含滞后因变量的动态模型**:

$$\log(\text{wage}_{it}) = \rho \log(\text{wage}_{i,t-1}) + c_i + u_{it}$$

因为 $\mathbf{x}_{i,t+1} = \log(\text{wage}_{it})$ 与 $u_{it}$ 相关（$E(\text{wage}_{it} \cdot u_{it}) = E(u_{it}^2) > 0$），严格外生性**必然不成立**。

**项目评估模型**: $u_{it}$（过去的工资冲击）可能影响未来的项目参与决策 $(prog_{i,t+1})$ → 严格外生性可能不成立。

---

## 3. 混合OLS (Pooled OLS)

### 3.1 估计

将模型写为 $y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + v_{it}$，其中 $v_{it} = c_i + u_{it}$，然后进行OLS。

### 3.2 一致性条件

需要 $E(\mathbf{x}_{it}' v_{it}) = 0, \; t = 1,...,T$。在 $E(\mathbf{x}_{it}' u_{it}) = 0$ 下，等价于：

$$E(\mathbf{x}_{it}' c_i) = 0, \quad t = 1,2,...,T$$

即**非观测效应与解释变量不相关**。

### 3.3 偏误后果

若 $\text{Cov}(\mathbf{x}_{it}, c_i) \neq 0$ → **混合OLS有偏且不一致**。

### 3.4 效率问题

即使一致，混合OLS也非有效。因为 $v_{it} = c_i + u_{it}$ 存在序列相关：

$$\text{Corr}(v_{it}, v_{is}) = \frac{\sigma_c^2}{\sigma_c^2 + \sigma_u^2} > 0, \quad s \neq t$$

且相关不随时间间隔 $|t-s|$ 增大而衰减。

---

## 4. 随机效应 (Random Effects, RE) 估计

### 4.1 RE 的三个核心假设

**假设 RE.1 (严格外生性 + 正交性)**:

$$\text{RE.1a: } E(u_{it} | \mathbf{x}_i, c_i) = 0, \quad t = 1,...,T$$

$$\text{RE.1b: } E(c_i | \mathbf{x}_i) = E(c_i) = 0$$

RE.1b 是**关键**——要求 $c_i$ 与所有时期的 $\mathbf{x}_{it}$ 均值独立。这是RE与FE的根本区别。

**假设 RE.2 (秩条件)**:

$$\text{rank } E(\mathbf{X}_i' \boldsymbol{\Omega}^{-1} \mathbf{X}_i) = K$$

**假设 RE.3 (方差结构)**:

(a) $E(\mathbf{u}_i \mathbf{u}_i' | \mathbf{x}_i, c_i) = \sigma_u^2 \mathbf{I}_T$ ——特异性误差条件同方差且序列无关

(b) $E(c_i^2 | \mathbf{x}_i) = \sigma_c^2$ ——$c_i$ 的条件方差为常数

### 4.2 RE的方差-协方差矩阵

在RE.1-RE.3下，复合误差 $\mathbf{v}_i = c_i \mathbf{j}_T + \mathbf{u}_i$ 的方差矩阵具有**随机效应结构**：

$$\boldsymbol{\Omega} \equiv E(\mathbf{v}_i \mathbf{v}_i') = \sigma_u^2 \mathbf{I}_T + \sigma_c^2 \mathbf{j}_T \mathbf{j}_T' \tag{10.31}$$

$$\boldsymbol{\Omega} = \begin{pmatrix}
\sigma_c^2 + \sigma_u^2 & \sigma_c^2 & \cdots & \sigma_c^2 \\
\sigma_c^2 & \sigma_c^2 + \sigma_u^2 & \cdots & \sigma_c^2 \\
\vdots & \vdots & \ddots & \vdots \\
\sigma_c^2 & \sigma_c^2 & \cdots & \sigma_c^2 + \sigma_u^2
\end{pmatrix}_{T \times T}$$

**关键性质**: $\boldsymbol{\Omega}$ 仅依赖于**两个参数** ($\sigma_c^2, \sigma_u^2$)，而非 $T(T+1)/2$ 个自由参数 → 对任何 $T$ 都简洁。

### 4.3 方差成分的估计

**Step 1**: 混合OLS得 $\hat{\boldsymbol{\beta}}$，计算 $\hat{v}_{it} = y_{it} - \mathbf{x}_{it}\hat{\boldsymbol{\beta}}$

**Step 2**: 估计 $\sigma_v^2 = \sigma_c^2 + \sigma_u^2$：

$$\hat{\sigma}_v^2 = \frac{1}{NT - K} \sum_{i=1}^{N} \sum_{t=1}^{T} \hat{v}_{it}^2$$

**Step 3**: 估计 $\sigma_c^2 = E(v_{it} v_{is})$ ($t \neq s$)：

$$\hat{\sigma}_c^2 = \frac{1}{NT(T-1)/2 - K} \sum_{i=1}^{N} \sum_{t=1}^{T-1} \sum_{s=t+1}^{T} \hat{v}_{it} \hat{v}_{is}$$

然后 $\hat{\sigma}_u^2 = \hat{\sigma}_v^2 - \hat{\sigma}_c^2$。

### 4.4 RE估计量 (FGLS)

$$\hat{\boldsymbol{\beta}}_{RE} = \left( \sum_{i=1}^{N} \mathbf{X}_i' \hat{\boldsymbol{\Omega}}^{-1} \mathbf{X}_i \right)^{-1} \left( \sum_{i=1}^{N} \mathbf{X}_i' \hat{\boldsymbol{\Omega}}^{-1} \mathbf{y}_i \right) \tag{10.34}$$

### 4.5 拟去均值 (Quasi-Demeaning) 变换

RE有一个等价的实现方式。定义：

$$\lambda = 1 - \sqrt{\frac{\sigma_u^2}{\sigma_u^2 + T\sigma_c^2}}, \quad 0 \leq \lambda \leq 1$$

**变换后的变量**: $y_{it}^* = y_{it} - \lambda \bar{y}_i$, $\mathbf{x}_{it}^* = \mathbf{x}_{it} - \lambda \bar{\mathbf{x}}_i$

RE估计量等价于 $y_{it}^*$ 对 $\mathbf{x}_{it}^*$ 的混合OLS。

| $\lambda = 0$ | $\lambda = 1$ |
|--------------|--------------|
| 混合OLS | 固定效应 (FE) |

当 $\hat{\lambda} \to 1$（$T$ 大或 $\sigma_c^2 / \sigma_u^2$ 大）时，RE和FE估计值趋于一致。

### 4.6 RE的性质

- 在RE.1和RE.2下：**一致** ($N \to \infty$)
- 在RE.3下：**渐近有效**（在 $E(\mathbf{v}_i | \mathbf{x}_i) = 0$ 的所有一致估计量中）
- 若RE.3不成立（异方差/序列相关），RE仍一致（若RE.1成立），但需用**稳健方差矩阵**

### 4.7 检验非观测效应的存在 (Breusch-Pagan LM检验)

$H_0: \sigma_c^2 = 0$（不存在非观测效应 → 混合OLS即可）

检验统计量基于 $\hat{v}_{it}\hat{v}_{is}$ ($t \neq s$) 的均值。拒绝 $H_0$ → 存在序列相关（可能来自 $c_i$，也可能来自其他来源）。

**注意**: 拒绝 $H_0$ 不应直接解释为"随机效应结构必须正确"——可能只是存在一般的序列相关。

---

## 5. 固定效应 (Fixed Effects, FE) 估计

### 5.1 假设

**假设 FE.1 (严格外生性)**:

$$E(u_{it} | \mathbf{x}_{i1},...,\mathbf{x}_{iT}, c_i) = 0, \quad t = 1,...,T$$

**注意**: FE**不要求** $E(c_i | \mathbf{x}_i) = 0$——$c_i$ 与 $\mathbf{x}_{it}$ 可以任意相关。这是FE比RE更稳健的关键。

**假设 FE.2 (秩条件)**:

$$\text{rank} \sum_{t=1}^{T} E(\ddot{\mathbf{x}}_{it}' \ddot{\mathbf{x}}_{it}) = \text{rank } E(\ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i) = K$$

要求去均值后的解释变量满秩 → 排除时间不变变量 → 要求每个 $\mathbf{x}_{itj}$ 在**至少部分个体**上有跨时变异。

**假设 FE.3 (效率条件)**:

$$E(\mathbf{u}_i \mathbf{u}_i' | \mathbf{x}_i, c_i) = \sigma_u^2 \mathbf{I}_T$$

即 $u_{it}$ 同方差且序列无关。

### 5.2 组内变换 (Within Transformation)

**Step 1 — 时间平均**:

$$\bar{y}_i = \bar{\mathbf{x}}_i \boldsymbol{\beta} + c_i + \bar{u}_i, \quad \bar{y}_i = \frac{1}{T} \sum_{t=1}^{T} y_{it}$$

**Step 2 — 去均值**:

$$y_{it} - \bar{y}_i = (\mathbf{x}_{it} - \bar{\mathbf{x}}_i)\boldsymbol{\beta} + (u_{it} - \bar{u}_i)$$

即 $\ddot{y}_{it} = \ddot{\mathbf{x}}_{it}\boldsymbol{\beta} + \ddot{u}_{it}$。

**$c_i$ 被完全消除**。

### 5.3 矩阵形式

定义**时间去均值矩阵**: $\mathbf{Q}_T = \mathbf{I}_T - \mathbf{j}_T(\mathbf{j}_T'\mathbf{j}_T)^{-1}\mathbf{j}_T' = \mathbf{I}_T - \frac{1}{T}\mathbf{j}_T\mathbf{j}_T'$

- $\mathbf{Q}_T \mathbf{y}_i = \ddot{\mathbf{y}}_i$
- $\mathbf{Q}_T \mathbf{X}_i = \ddot{\mathbf{X}}_i$
- $\mathbf{Q}_T(c_i \mathbf{j}_T) = \mathbf{0}$（$c_i$ 被消除）

### 5.4 FE估计量 (组内估计量)

$$\hat{\boldsymbol{\beta}}_{FE} = \left( \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i \right)^{-1} \left( \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \ddot{\mathbf{y}}_i \right)$$
$$= \left( \sum_{i=1}^{N} \sum_{t=1}^{T} \ddot{\mathbf{x}}_{it}' \ddot{\mathbf{x}}_{it} \right)^{-1} \left( \sum_{i=1}^{N} \sum_{t=1}^{T} \ddot{\mathbf{x}}_{it}' \ddot{y}_{it} \right) \tag{10.48}$$

### 5.5 虚拟变量回归

FE等价于为每个个体引入虚拟变量进行OLS（**虚拟变量估计量**）：

$$y_{it} = \alpha_1 d_i^1 + \alpha_2 d_i^2 + \cdots + \alpha_N d_i^N + \mathbf{x}_{it}\boldsymbol{\beta} + u_{it}$$

$\hat{\boldsymbol{\beta}}$ 与 $\hat{\boldsymbol{\beta}}_{FE}$ 相同。但 $\hat{\alpha}_i$ **不一致**（因为每增加一个个体，就增加一个参数 → 信息不累积）。

### 5.6 去均值误差的方差结构

即使 $u_{it}$ 满足FE.3（同方差+序列无关），$\ddot{u}_{it} = u_{it} - \bar{u}_i$ 也**不是同方差的且存在负的序列相关**：

$$E(\ddot{u}_{it}^2) = \sigma_u^2 \left(1 - \frac{1}{T}\right) \tag{10.51}$$

$$E(\ddot{u}_{it} \ddot{u}_{is}) = -\frac{\sigma_u^2}{T}, \quad t \neq s$$

$$\text{Corr}(\ddot{u}_{it}, \ddot{u}_{is}) = -\frac{1}{T-1}$$

### 5.7 渐近方差估计

在FE.3下：

$$\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}_{FE}) = \hat{\sigma}_u^2 \left( \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i \right)^{-1}$$

$$\hat{\sigma}_u^2 = \frac{1}{N(T-1) - K} \sum_{i=1}^{N} \sum_{t=1}^{T} \hat{\ddot{u}}_{it}^2, \quad \hat{\ddot{u}}_{it} = \ddot{y}_{it} - \ddot{\mathbf{x}}_{it}\hat{\boldsymbol{\beta}}_{FE}$$

**为什么自由度是 $N(T-1)-K$？** 组内变换消耗了 $N$ 个自由度（$\bar{y}_i$ 的计算），有效观测数为 $N(T-1)$。

### 5.8 稳健方差矩阵 (Arellano, 1987)

当 $u_{it}$ 存在异方差或序列相关时：

$$\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}_{FE}) = (\ddot{\mathbf{X}}'\ddot{\mathbf{X}})^{-1} \left[ \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \hat{\ddot{\mathbf{u}}}_i \hat{\ddot{\mathbf{u}}}_i' \ddot{\mathbf{X}}_i \right] (\ddot{\mathbf{X}}'\ddot{\mathbf{X}})^{-1} \tag{10.59}$$

允许任意形式的异方差和序列相关（$T$ 相对于 $N$ 较小时有效）。

### 5.9 时间不变变量与时间虚拟变量的交互

虽然无法直接估计时间不变变量（如性别）的主效应，但可以估计其与时间虚拟变量的交互效应：

$$y_{it} = \mu_1 + \mu_2 d2_t + \cdots + \mu_T dT_t + \mathbf{z}_i \boldsymbol{\gamma}_1 + d2_t \cdot \mathbf{z}_i \boldsymbol{\gamma}_2 + \cdots + dT_t \cdot \mathbf{z}_i \boldsymbol{\gamma}_T + \mathbf{w}_{it}\boldsymbol{\delta} + c_i + u_{it}$$

- $\boldsymbol{\gamma}_1$ **不可识别**（与 $c_i$ 无法区分）
- $\boldsymbol{\gamma}_2, ..., \boldsymbol{\gamma}_T$ **可识别**——捕捉了时间不变变量偏效应随时间的变化

**例子**: $\mathbf{z}_i = \text{female}$（女性虚拟变量）→ 检验性别工资差距是否随时间变化。

---

## 6. 一阶差分 (First Differencing, FD) 估计

### 6.1 假设

**假设 FD.1 (严格外生性)**: 与FE.1相同

$$\Delta y_{it} = \Delta \mathbf{x}_{it} \boldsymbol{\beta} + \Delta u_{it}, \quad t = 2, 3, ..., T \tag{10.63}$$

$c_i$ 被差分消除。损失第一个时期，剩 $T-1$ 个时期。

### 6.2 FD估计量

$$\hat{\boldsymbol{\beta}}_{FD} = \left( \sum_{i=1}^{N} \sum_{t=2}^{T} \Delta \mathbf{x}_{it}' \Delta \mathbf{x}_{it} \right)^{-1} \left( \sum_{i=1}^{N} \sum_{t=2}^{T} \Delta \mathbf{x}_{it}' \Delta y_{it} \right) \tag{10.65}$$

### 6.3 效率条件与检验

**假设 FD.3**: $E(\mathbf{e}_i \mathbf{e}_i' | \mathbf{x}_{i1},...,\mathbf{x}_{iT}, c_i) = \sigma_e^2 \mathbf{I}_{T-1}$，其中 $\mathbf{e}_i$ 为 $(T-1) \times 1$ 向量，$e_{it} = \Delta u_{it}$。

FD.3 要求差分后的误差 $e_{it}$ 同方差且**序列无关**。

#### 为什么 FE.3 成立时 FD 误差是序列相关的？

**FE.3**: $u_{it}$ 序列无关，$\text{Var}(u_{it}) = \sigma_u^2$。计算 $\text{Cov}(e_{it}, e_{i,t-1})$:

$$\begin{aligned}
e_{it} &= \Delta u_{it} = u_{it} - u_{i,t-1} \\
e_{i,t-1} &= u_{i,t-1} - u_{i,t-2}
\end{aligned}$$

$$\begin{aligned}
\text{Cov}(e_{it}, e_{i,t-1}) &= \text{Cov}(u_{it} - u_{i,t-1},\; u_{i,t-1} - u_{i,t-2}) \\
&= \underbrace{\text{Cov}(u_{it}, u_{i,t-1})}_{0} - \underbrace{\text{Cov}(u_{it}, u_{i,t-2})}_{0} - \text{Cov}(u_{i,t-1}, u_{i,t-1}) + \underbrace{\text{Cov}(u_{i,t-1}, u_{i,t-2})}_{0} \\
&= -\text{Var}(u_{i,t-1}) = -\sigma_u^2
\end{aligned}$$

$\text{Var}(e_{it}) = \text{Var}(u_{it} - u_{i,t-1}) = \sigma_u^2 + \sigma_u^2 = 2\sigma_u^2$

$$\rho_1 \equiv \text{Corr}(e_{it}, e_{i,t-1}) = \frac{-\sigma_u^2}{2\sigma_u^2} = -0.5$$

**直觉**: 相邻两期的差分共享同一个 $u_{i,t-1}$——$e_{it}$ 含 $-u_{i,t-1}$，$e_{i,t-1}$ 含 $+u_{i,t-1}$ → 符号相反 → 负相关，且恰好为 $-0.5$。

#### 为什么 FD.3 成立时原误差是序列相关的？

FD.3 要求 $\text{Cov}(e_{it}, e_{i,t-1}) = 0$。反推 $u_{it}$ 的性质:

$$\text{Cov}(e_{it}, e_{i,t-1}) = \text{Cov}(u_{it} - u_{i,t-1},\; u_{i,t-1} - u_{i,t-2}) = -\text{Cov}(u_{i,t-1}, u_{i,t-1}) + \text{Cov}(u_{it}, u_{i,t-2}) + \cdots = 0$$

该条件满足当 $u_{it}$ 遵循**随机游走**: $u_{it} = u_{i,t-1} + \varepsilon_{it}$，$\varepsilon_{it}$ 为白噪声。此时 $e_{it} = u_{it} - u_{i,t-1} = \varepsilon_{it}$，自然序列无关。

**FD.3 和 FE.3 是互斥的**（除非 $\sigma_u^2 = 0$）:

| 假设 | 对 $u_{it}$ 的要求 | 对 $e_{it} = \Delta u_{it}$ 的要求 | $\rho_1 = \text{Corr}(e_{it}, e_{i,t-1})$ |
|------|-------------------|----------------------------------|------------------------------------------|
| **FE.3** | 序列无关, $\text{Var} = \sigma_u^2$ | 必然是序列相关的 | $-0.5$ |
| **FD.3** | 随机游走（必然序列相关） | 序列无关, $\text{Var} = \sigma_e^2$ | $0$ |
| 两者都成立 | $\sigma_u^2 = 0$（平凡情况） | $e_{it} = 0$ | 无定义 |

**检验方法**: 回归 $\hat{e}_{it}$ 对 $\hat{e}_{i,t-1}$ ($t = 3,...,T$)，检验 $H_0: \rho_1 = 0$（需 $T \geq 3$）。通常的t统计量在FD.3下渐近有效。

**实践含义**: 若 $\hat{\rho}_1 \neq 0$ 且 $\neq -0.5$ → FE.3和FD.3都可能被拒绝 → 应使用**稳健方差矩阵**，因为它对 $u_{it}$ 的序列相关结构不做任何假设。

---

## 7. FE、FD 和 RE 的比较

### 7.1 FE vs FD 综合对比

| 维度 | FE | FD |
|------|-----|-----|
| 变换 | $y_{it} - \bar{y}_i$ | $y_{it} - y_{i,t-1}$ |
| $T=2$ 时 | 等价 | 等价 |
| $T>2$ 且 $u_{it}$ 序列无关 | **更有效** | 次有效 |
| $T>2$ 且 $u_{it}$ 随机游走 | 次有效 | **更有效** |
| 严格外生性不成立时的不一致性 | $O(T^{-1})$（随 $T$ 增大而减小） | $O(1)$（不随 $T$ 减小） |
| 序贯外生下可用IV？ | 否（$\bar{\mathbf{x}}_i$ 含未来值） | **是**（滞后 $\mathbf{x}$ 为有效IV） |
| 实现难度 | 需专用软件 | 普通OLS即可（手动差分） |

### 7.2 Hausman检验 (RE vs FE)

**检验问题**: 是否存在 $\text{Cov}(\mathbf{x}_{it}, c_i) \neq 0$？

$$H_0: E(c_i | \mathbf{x}_i) = 0$$

**检验统计量**:

$$H = (\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE})' [\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}_{FE}) - \widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}_{RE})]^{-1} (\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE}) \xrightarrow{d} \chi^2_M$$

其中 $M$ 为时变解释变量的个数。在 $H_0$ 和 RE.3 下成立。

**单参数版本**: $t = (\hat{\beta}_{j,FE} - \hat{\beta}_{j,RE}) / \sqrt{\text{se}(\hat{\beta}_{j,FE})^2 - \text{se}(\hat{\beta}_{j,RE})^2}$

**关键**: $\text{Avar}(\hat{\boldsymbol{\beta}}_{FE}) - \text{Avar}(\hat{\boldsymbol{\beta}}_{RE})$ 为正定矩阵——FE比RE有更大的渐近方差（RE在$H_0$下更有效）。

**重要警告**: Hausman检验维持**严格外生性**假设（RE.1a）。若严格外生性不成立，检验结果不可靠。检验对RE.3不成立无系统性检验力。

### 7.3 何时偏好FD？

1. **严格外生性可能不成立** → FD与滞后$\mathbf{x}$为IV兼容（详见Ch6）
2. **手动实现简单**（无需专用面板软件）
3. **$u_{it}$ 近似随机游走**时FD更有效

---

## 8. 政策分析中的FE/FD

### 8.1 FD的政策分析模型

$$y_{it} = \mathbf{z}_{it}\boldsymbol{\gamma} + \delta w_{it} + c_i + u_{it}$$

差分后: $\Delta y_{it} = \Delta \mathbf{z}_{it}\boldsymbol{\gamma} + \delta \Delta w_{it} + \Delta u_{it}$

$w_{it}$ 为政策变量（如是否参与项目）。

### 8.2 项目评估中的关键问题

**Example 10.1 (职业培训对工资的影响)**:

$$\log(\text{wage}_{it}) = \mu_t + \mathbf{z}_{it}\boldsymbol{\gamma} + \delta_1 \text{prog}_{it} + c_i + u_{it}$$

- 在 $t=1$ 时所有人都未参与 ($\text{prog}_{i1}=0$)
- 在 $t=2$ 时一些人参与了
- $c_i$ 代表能力——可能同时影响工资和参与决策（**自选择问题**）
- 严格外生性：$u_{it}$ 与过去/当前/未来的 $\text{prog}_{is}$ 无关？→ 若过去的工资冲击影响未来参与决策，**严格外生性不成立** → 需Ch6方法

---

## 9. 检验严格外生性

**T = 2 时**: 在FD方程中添加 $\mathbf{x}_{i2}$，检验其是否显著。

**T > 2 时**: 在FD方程中添加 $\mathbf{w}_{i,t+1}$（未来 $\mathbf{x}$ 的子集），检验 $H_0: \boldsymbol{\delta} = \mathbf{0}$。

**使用FE**: 估计 $y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + \mathbf{w}_{i,t+1}\boldsymbol{\delta} + c_i + u_{it}$ ($t=1,...,T-1$)，检验 $\boldsymbol{\delta} = \mathbf{0}$。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + c_i + u_{it}$ | 基准非观测效应模型 (10.11) |
| $\boldsymbol{\Omega} = \sigma_u^2\mathbf{I}_T + \sigma_c^2\mathbf{j}_T\mathbf{j}_T'$ | RE方差结构 (10.31) |
| $\ddot{y}_{it} = y_{it} - \bar{y}_i$ | FE组内变换 |
| $\hat{\boldsymbol{\beta}}_{FE} = (\sum \ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i)^{-1}(\sum \ddot{\mathbf{X}}_i' \ddot{\mathbf{y}}_i)$ | FE估计量 (10.48) |
| $E(\ddot{u}_{it}^2) = \sigma_u^2(1 - 1/T)$ | 去均值误差方差 (10.51) |
| $\Delta y_{it} = y_{it} - y_{i,t-1}$ | FD变换 (10.63) |
| $\lambda = 1 - \sqrt{\sigma_u^2/(\sigma_u^2 + T\sigma_c^2)}$ | RE拟去均值参数 |
| $H = (\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE})'[\widehat{\text{Avar}}_{FE} - \widehat{\text{Avar}}_{RE}]^{-1}(\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE})$ | Hausman统计量 |

---

## 附录 A: RE方差-协方差矩阵 $\boldsymbol{\Omega}$ 的详细推导与讲解

### A.1 从复合误差出发

模型为 $y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + c_i + u_{it}$。对于个体 $i$，将其 $T$ 个时期的方程堆叠：

$$\mathbf{y}_i = \mathbf{X}_i \boldsymbol{\beta} + \mathbf{v}_i, \quad \mathbf{v}_i = c_i \mathbf{j}_T + \mathbf{u}_i$$

其中 $\mathbf{j}_T = (1, 1, ..., 1)'$ 为 $T \times 1$ 全1向量。

$\mathbf{v}_i$ 是 $T \times 1$ 复合误差向量，$\mathbf{u}_i = (u_{i1}, u_{i2}, ..., u_{iT})'$。

### A.2 逐元素推导

$$\boldsymbol{\Omega} = E[\mathbf{v}_i \mathbf{v}_i' | \mathbf{x}_i] = E[(c_i \mathbf{j}_T + \mathbf{u}_i)(c_i \mathbf{j}_T + \mathbf{u}_i)' | \mathbf{x}_i]$$

展开期望：

$$\boldsymbol{\Omega} = E[c_i^2 | \mathbf{x}_i] \cdot \mathbf{j}_T \mathbf{j}_T' + E[c_i \mathbf{u}_i' | \mathbf{x}_i] \cdot \mathbf{j}_T + E[\mathbf{u}_i c_i | \mathbf{x}_i] \cdot \mathbf{j}_T' + E[\mathbf{u}_i \mathbf{u}_i' | \mathbf{x}_i]$$

**逐个分析这四项**:

| 项 | 表达式 | RE.3下的值 | 原因 |
|----|--------|-----------|------|
| ① | $E[c_i^2 \mid \mathbf{x}_i] \cdot \mathbf{j}_T \mathbf{j}_T'$ | $\sigma_c^2 \cdot \mathbf{j}_T \mathbf{j}_T'$ | RE.3(b): $c_i$ 条件方差为常数 |
| ② | $E[c_i \mathbf{u}_i' \mid \mathbf{x}_i] \cdot \mathbf{j}_T$ | $\mathbf{0}$ | RE.1a: $E(\mathbf{u}_i \mid \mathbf{x}_i, c_i) = \mathbf{0}$ → $E(c_i \mathbf{u}_i' \mid \mathbf{x}_i) = E[c_i \cdot E(\mathbf{u}_i' \mid \mathbf{x}_i, c_i) \mid \mathbf{x}_i] = \mathbf{0}$ |
| ③ | $E[\mathbf{u}_i c_i \mid \mathbf{x}_i] \cdot \mathbf{j}_T'$ | $\mathbf{0}$ | 同上，对称 |
| ④ | $E[\mathbf{u}_i \mathbf{u}_i' \mid \mathbf{x}_i]$ | $\sigma_u^2 \cdot \mathbf{I}_T$ | RE.3(a): $u_{it}$ 同方差且序列无关 |

因此：

$$\boxed{\boldsymbol{\Omega} = \sigma_u^2 \mathbf{I}_T + \sigma_c^2 \mathbf{j}_T \mathbf{j}_T'} \tag{A.1}$$

### A.3 矩阵展开：每个元素的含义

$\mathbf{j}_T \mathbf{j}_T'$ 是一个 **$T \times T$ 全1矩阵**（每个元素都等于1）。所以：

$$\boldsymbol{\Omega} = \sigma_u^2 \begin{pmatrix}
1 & 0 & \cdots & 0 \\
0 & 1 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 1
\end{pmatrix} + \sigma_c^2 \begin{pmatrix}
1 & 1 & \cdots & 1 \\
1 & 1 & \cdots & 1 \\
\vdots & \vdots & \ddots & \vdots \\
1 & 1 & \cdots & 1
\end{pmatrix} = \begin{pmatrix}
\sigma_c^2 + \sigma_u^2 & \sigma_c^2 & \cdots & \sigma_c^2 \\
\sigma_c^2 & \sigma_c^2 + \sigma_u^2 & \cdots & \sigma_c^2 \\
\vdots & \vdots & \ddots & \vdots \\
\sigma_c^2 & \sigma_c^2 & \cdots & \sigma_c^2 + \sigma_u^2
\end{pmatrix}$$

**对角元素** ($t = s$):

$$\omega_{tt} = E(v_{it}^2) = E[(c_i + u_{it})^2] = E(c_i^2) + 2E(c_i u_{it}) + E(u_{it}^2) = \sigma_c^2 + \sigma_u^2$$

| 成分 | 来源 | 含义 |
|------|------|------|
| $\sigma_c^2$ | $c_i$ 的方差 | 个体间的持久差异对总波动的贡献 |
| $\sigma_u^2$ | $u_{it}$ 的方差 | 特异性暂时冲击对总波动的贡献 |

**非对角元素** ($t \neq s$):

$$\omega_{ts} = E(v_{it} v_{is}) = E[(c_i + u_{it})(c_i + u_{is})] = E(c_i^2) + E(c_i u_{it}) + E(c_i u_{is}) + E(u_{it} u_{is}) = \sigma_c^2$$

| 成分 | 来源 | 含义 |
|------|------|------|
| $E(c_i^2) = \sigma_c^2$ | $c_i$ 的持久影响 | 同一 $c_i$ 出现在所有时期 → 产生跨期相关 |
| $E(u_{it} u_{is}) = 0$ | $u_{it}$ 序列无关 (RE.3a) | 特异性冲击不产生序列相关 |
| $E(c_i u_{it}) = 0$ | 交叉项 (RE.1a) | $c_i$ 和 $u_{it}$ 不相关 |

### A.4 核心性质

#### (1) 复合对称 (Compound Symmetry)

$\boldsymbol{\Omega}$ 的所有对角元素相等，所有非对角元素也相等。这种结构称为**复合对称**或**可交换相关结构**。

#### (2) 仅含两个自由参数

无论 $T$ 多大，$\boldsymbol{\Omega}$ 仅由 $\sigma_c^2$ 和 $\sigma_u^2$ 两个参数完全确定——而非无约束的 $T(T+1)/2$ 个参数。这对于 $T$ 较大时GLS估计至关重要。

#### (3) 复合误差的相关系数

$$\text{Corr}(v_{it}, v_{is}) = \frac{\sigma_c^2}{\sigma_c^2 + \sigma_u^2} \equiv \rho, \quad \forall t \neq s$$

- $\rho$ 是**任何两个不同时期**复合误差之间的相关系数
- $\rho$ 也是 $c_i$ 对总误差方差的**贡献率**
- $\rho > 0$（只要 $\sigma_c^2 > 0$）
- $\rho$ **不随** $|t-s|$ 增大而衰减——这是 $c_i$ 持久性导致的（与AR(1)序列相关不同）

#### (4) 特征值与谱分解

$\boldsymbol{\Omega}$ 的特征值：

$$\begin{aligned}
\lambda_1 &= \sigma_u^2 + T\sigma_c^2 \quad &\text{(重数 1, 特征向量 } \mathbf{j}_T/\sqrt{T}\text{)} \\
\lambda_2 &= \sigma_u^2 \quad &\text{(重数 } T-1\text{, 任意与 } \mathbf{j}_T \text{ 正交的向量)}
\end{aligned}$$

**含义**: 最大的特征值 $\lambda_1$ 对应"均值方向"（组间变异），而 $\lambda_2 = \sigma_u^2$ 对应所有与均值正交的方向（组内变异）。

### A.5 $\boldsymbol{\Omega}$ 的逆

$$\boldsymbol{\Omega}^{-1} = \frac{1}{\sigma_u^2} \left( \mathbf{I}_T - \frac{\sigma_c^2}{\sigma_u^2 + T\sigma_c^2} \mathbf{j}_T \mathbf{j}_T' \right)$$

**推导思路**: 利用 $\boldsymbol{\Omega} = \sigma_u^2 \mathbf{I}_T + \sigma_c^2 \mathbf{j}_T \mathbf{j}_T'$ 是一个"恒等矩阵 + 秩1修正"的结构，通过 Sherman-Morrison 公式求逆。

令 $\psi = \frac{\sigma_u^2}{\sigma_u^2 + T\sigma_c^2}$，则 $\boldsymbol{\Omega}^{-1}$ 可写为更简洁的形式：

$$\boldsymbol{\Omega}^{-1} = \frac{1}{\sigma_u^2} \left[ \mathbf{I}_T - (1-\psi) \cdot \frac{1}{T} \cdot \mathbf{j}_T \mathbf{j}_T' \right]$$

由于 $\mathbf{P}_T = \frac{1}{T}\mathbf{j}_T \mathbf{j}_T'$ 是时间平均投影矩阵，$\mathbf{Q}_T = \mathbf{I}_T - \mathbf{P}_T$ 是去均值矩阵：

$$\boldsymbol{\Omega}^{-1} = \frac{1}{\sigma_u^2} \left( \mathbf{Q}_T + \psi \cdot \mathbf{P}_T \right)$$

### A.6 $\boldsymbol{\Omega}^{-1/2}$ 与拟去均值参数 $\lambda$

GLS需要的是 $\boldsymbol{\Omega}^{-1/2}$（而非 $\boldsymbol{\Omega}^{-1}$）。由 $\boldsymbol{\Omega}$ 的谱分解可得：

$$\boldsymbol{\Omega}^{-1/2} = \frac{1}{\sigma_u} \left( \mathbf{I}_T - \frac{\lambda}{T} \mathbf{j}_T \mathbf{j}_T' \right) = \frac{1}{\sigma_u} (\mathbf{I}_T - \lambda \mathbf{P}_T)$$

其中：

$$\lambda = 1 - \sqrt{\psi} = 1 - \sqrt{\frac{\sigma_u^2}{\sigma_u^2 + T\sigma_c^2}}, \quad 0 \leq \lambda \leq 1$$

$\boldsymbol{\Omega}^{-1/2} \mathbf{y}_i$ 的第 $t$ 个元素为：

$$y_{it} - \lambda \bar{y}_i$$

这就是**拟去均值变换**的由来。

| 参数关系 | 公式 | 值域 |
|---------|------|------|
| 非对角相关 | $\rho = \frac{\sigma_c^2}{\sigma_c^2 + \sigma_u^2}$ | $[0, 1)$ |
| 方差比 | $\psi = \frac{\sigma_u^2}{\sigma_u^2 + T\sigma_c^2}$ | $(0, 1]$ |
| 拟去均值参数 | $\lambda = 1 - \sqrt{\psi}$ | $[0, 1]$ |

**极端情况**:

| $\lambda \to 0$ ($\sigma_c^2 \to 0$ 或 $T \to 1$) | $\text{RE} \to \text{混合OLS}$ |
| $\lambda \to 1$ ($\sigma_c^2 / \sigma_u^2 \to \infty$ 或 $T \to \infty$) | $\text{RE} \to \text{FE}$ |

### A.7 数值示例

设 $\sigma_c^2 = 4$, $\sigma_u^2 = 1$, $T = 3$：

$$\boldsymbol{\Omega} = \begin{pmatrix}
5 & 4 & 4 \\
4 & 5 & 4 \\
4 & 4 & 5
\end{pmatrix}$$

$$\rho = \frac{4}{4+1} = 0.8, \quad \psi = \frac{1}{1+3\cdot 4} = \frac{1}{13}, \quad \lambda = 1 - \sqrt{\frac{1}{13}} \approx 0.723$$

复合误差之间的相关系数为 $0.8$——非常高，意味着 $c_i$ 占误差波动的主导地位。$\lambda \approx 0.723$ 意味着 RE 估计量消除了个体均值的约 72%，比混合 OLS ($\lambda=0$) 更接近 FE ($\lambda=1$)。

### A.8 与OLS和FE方差结构的关系

| 估计方法 | 隐含的误差方差结构 | $\boldsymbol{\Omega}$ 的参数约束 |
|---------|-------------------|------------------------------|
| 混合OLS | $\boldsymbol{\Omega} \propto \mathbf{I}_T$ | $\sigma_c^2 = 0$ (无序列相关) |
| RE | $\boldsymbol{\Omega} = \sigma_u^2 \mathbf{I}_T + \sigma_c^2 \mathbf{j}_T \mathbf{j}_T'$ | 复合对称（两个参数） |
| FE（组内变换后） | $\ddot{\boldsymbol{\Omega}} = \sigma_u^2 \mathbf{Q}_T$ | 仅 $\sigma_u^2$（$\sigma_c^2$ 被消除） |
| 通用FGLS | $\boldsymbol{\Omega}$ 无约束 (任意 $T(T+1)/2$ 参数) | 每个方差和协方差自由估计 |

RE是混合OLS和通用FGLS之间的折中——既不像OLS那样忽略序列相关，也不像通用FGLS那样需要估计 $T(T+1)/2$ 个参数。
