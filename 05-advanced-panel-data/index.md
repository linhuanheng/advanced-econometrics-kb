---
title: 高级面板数据专题
type: framework
domain: econometrics
source: Chapter VI — More Topics in Linear Unobserved Effects Models
created: 2026-05-19
updated: 2026-05-24
---

# Ch6: 高级面板数据专题

## 概述

Ch5假设严格外生性（$u_{it}$ 与**所有时期**的 $\mathbf{x}_{is}$ 无关）。本章放宽该假设，处理仅满足**序贯外生性**的模型，并涵盖个体特定斜率、Chamberlain方法、Hausman-Taylor模型等。

---

## 1. 序贯外生性与严格外生性

### 1.1 严格外生性 (Strict Exogeneity) 回顾

$$E(u_{it} | \mathbf{x}_{i1}, ..., \mathbf{x}_{iT}, c_i) = 0, \quad t = 1,...,T$$

→ 所有时期的 $\mathbf{x}$ 与所有时期的 $u$ 都不相关。在动态面板、反馈模型中**常不成立**。

### 1.2 序贯外生性 (Sequential Exogeneity)

**Chamberlain (1992b)** 提出序贯矩条件：

$$E(u_{it} | \mathbf{x}_{it}, \mathbf{x}_{i,t-1}, ..., \mathbf{x}_{i1}, c_i) = 0, \quad t = 1,2,...,T \tag{11.2}$$

| 对比维度 | 严格外生性 | 序贯外生性 |
|---------|-----------|-----------|
| 条件集 | 所有时期 $\mathbf{x}_{i1},...,\mathbf{x}_{iT}$ | 仅当期和过去 $\mathbf{x}_{it},...,\mathbf{x}_{i1}$ |
| 未来 $\mathbf{x}$ 与当期 $u$ | 必须无关 | **允许相关** |
| 对动态模型的适用性 | 不适用 | 适用 |
| 对理性预期检验 | 可能不满足 | 通常满足 (Keane & Runkle 1992) |

**关键含义**: $E(y_{it} | \mathbf{x}_{it}, \mathbf{x}_{i,t-1}, ..., \mathbf{x}_{i1}, c_i) = E(y_{it} | \mathbf{x}_{it}, c_i) = \mathbf{x}_{it}\boldsymbol{\beta} + c_i$

→ 一旦控制当期 $\mathbf{x}$ 和 $c_i$，过去的 $\mathbf{x}$ 不再提供额外信息。**未来 $\mathbf{x}$ 与当期 $y$ 的关系则不受限制**。

### 1.3 序贯外生下的正交条件

$$E(\mathbf{x}_{is}' u_{it}) = 0, \quad s = 1, 2, ..., t$$

$$\text{但对 } s > t: \; E(\mathbf{x}_{is}' u_{it}) \neq 0 \text{（可能）}$$

对差分方程，关键的正交条件为：

$$E(\mathbf{x}_{is}' \Delta u_{it}) = 0, \quad s = 1, 2, ..., t-1$$

---

## 2. 严格外生性不成立时 FE 和 FD 的渐近偏误

### 2.1 FE 的不一致性

$$\text{plim } \hat{\boldsymbol{\beta}}_{FE} = \boldsymbol{\beta} + \left[ \frac{1}{T} \sum_{t=1}^{T} E(\ddot{\mathbf{x}}_{it}' \ddot{\mathbf{x}}_{it}) \right]^{-1} \cdot \frac{1}{T} \sum_{t=1}^{T} E(\ddot{\mathbf{x}}_{it}' u_{it})$$

在序贯外生下，$E(\ddot{\mathbf{x}}_{it}' u_{it}) = -E(\bar{\mathbf{x}}_i' u_{it}) \neq 0$（因为 $\bar{\mathbf{x}}_i$ 包含未来的 $\mathbf{x}_{is}$）。

**偏误量级**: 在"矩有界"和"弱相依"假设下：

$$\left| E(\bar{x}_{ij} \bar{u}_i) \right| \leq [\text{Var}(\bar{x}_{ij}) \cdot \text{Var}(\bar{u}_i)]^{1/2} = O(T^{-1})$$

$$\text{plim } \hat{\boldsymbol{\beta}}_{FE} = \boldsymbol{\beta} + O(T^{-1})$$

→ FE的不一致性**随 $T$ 增大而减小**。

**Hsiao的结果**: 对AR(1)面板，$|\rho| < 1$ 时偏误为 $O(T^{-1})$。但**若 $\rho$ 接近1**（$\mathbf{x}_{it}$ 高度持久），即使 $T$ 较大，FE偏误也**可能相当可观**。

### 2.2 FD 的不一致性

$$\text{plim } \hat{\boldsymbol{\beta}}_{FD} = \boldsymbol{\beta} + \left[ \frac{1}{T-1} \sum_{t=2}^{T} E(\Delta \mathbf{x}_{it}' \Delta \mathbf{x}_{it}) \right]^{-1} \cdot \frac{1}{T-1} \sum_{t=2}^{T} E(\Delta \mathbf{x}_{it}' \Delta u_{it})$$

**关键推导**:

$$\begin{aligned}
E(\Delta \mathbf{x}_{it}' \Delta u_{it}) &= E(\mathbf{x}_{it}' u_{it}) - E(\mathbf{x}_{it}' u_{i,t-1}) - E(\mathbf{x}_{i,t-1}' u_{it}) + E(\mathbf{x}_{i,t-1}' u_{i,t-1}) \\
&= -E(\mathbf{x}_{i,t-1}' u_{it}) \quad \text{（序贯外生下仅此项非零）}
\end{aligned}$$

在**平稳性**假设下，$E(\mathbf{x}_{i,t-1}' u_{it})$ 不依赖于 $t$ 或 $T$ → 是一个**常数**。

$$\text{plim } \hat{\boldsymbol{\beta}}_{FD} = \boldsymbol{\beta} + O(1)$$

→ FD的不一致性**不随 $T$ 增大而消失**。

### 2.3 对实践的指导

| 情况 | 偏好 |
|------|------|
| 严格外生可能不成立，$T > 2$ | FE（偏误 $O(T^{-1})$ < FD的 $O(1)$） |
| 序贯外生成立，可用IV | **FD**（滞后 $\mathbf{x}$ 是有效工具变量） |

---

## 3. 动态面板模型

### 3.1 AR(1) 含非观测效应

**Example 11.1**:

$$y_{it} = \mathbf{z}_{it}\boldsymbol{\gamma} + \rho y_{i,t-1} + c_i + u_{it} \tag{11.4}$$

令 $\mathbf{x}_{it} = (\mathbf{z}_{it}, y_{i,t-1})$。在序贯外生下：

$$E(y_{it} | \mathbf{z}_{it}, y_{i,t-1}, \mathbf{z}_{i,t-1}, ..., \mathbf{z}_{i1}, y_{i0}, c_i) = E(y_{it} | \mathbf{z}_{it}, y_{i,t-1}, c_i) = \mathbf{z}_{it}\boldsymbol{\gamma} + \rho y_{i,t-1} + c_i$$

这是**以 $c_i$ 为条件的动态完备性**——一期滞后足以捕捉动态。

**检验 $H_0: \rho = 0$**: 控制 $c_i$ 和当前及过去 $\mathbf{z}$ 后，$y_{i,t-1}$ 是否仍预测 $y_{it}$？

- $\rho \neq 0$ → **状态依赖 (State Dependence)**: 当前状态受上一期状态影响，即使在控制了 $c_i$ 后

### 3.2 为什么严格外生性必然不成立？

$$E(y_{it} \cdot u_{it}) = E(u_{it}^2) > 0$$

而 $\mathbf{x}_{i,t+1}$ 包含 $y_{it}$（作为滞后因变量）→ $E(\mathbf{x}_{i,t+1}' u_{it}) \neq 0$

→ 严格外生性**在含滞后因变量的非观测效应模型中永远不成立**。

### 3.3 标准序贯外生假设

$$E(u_{it} | y_{i,t-1}, y_{i,t-2}, ..., y_{i0}, c_i) = 0 \tag{10.19}$$

→ 仅当前和过去的 $y$ 在条件集中——所有动态已被一阶滞后捕捉。

### 3.4 $y_{i,t-1}$ 与 $c_i$ 必然相关

因为 $y_{i,t-1}$ 在 $t-1$ 时期是左边变量 → $y_{i,t-1}$ 包含 $c_i$ 的成分 → $\text{Cov}(y_{i,t-1}, c_i) \neq 0$

→ **混合OLS和RE都不一致**。

---

## 4. 差分方程的IV估计

### 4.1 FD后的矩条件

对 (11.4) 取一阶差分：

$$\Delta y_{it} = \Delta \mathbf{z}_{it}\boldsymbol{\gamma} + \rho \Delta y_{i,t-1} + \Delta u_{it}$$

$\Delta y_{i,t-1}$ 与 $\Delta u_{it}$ 相关（因为 $y_{i,t-1}$ 与 $u_{i,t-1}$ 相关）。

但在序贯外生下，$y_{i,t-2}$ 与 $\Delta u_{it}$ **不相关** → $y_{i,t-2}$ 可作为 $\Delta y_{i,t-1}$ 的工具变量。

### 4.2 一般模型：差分方程的IV

对 $\Delta y_{it} = \Delta \mathbf{x}_{it}\boldsymbol{\beta} + \Delta u_{it}$ ($t=3,...,T$)：

**可用工具变量**: $\mathbf{X}_{i,t-1}^0 = (\mathbf{x}_{i1}, \mathbf{x}_{i2}, ..., \mathbf{x}_{i,t-1})$

**为什么有效？** $E(\mathbf{x}_{is}' \Delta u_{it}) = 0$ 对 $s = 1,...,t-1$（序贯外生下的矩条件）。

**确切识别**: 使用 $\Delta \mathbf{x}_{i,t-1}$（或等价地 $\mathbf{x}_{i,t-1}$）为工具变量。需满足：

1. $E(\Delta \mathbf{x}_{i,t-1}' \Delta u_{it}) = 0$（序贯外生保证）
2. $\text{rank } E(\Delta \mathbf{x}_{i,t-1}' \Delta \mathbf{x}_{it}) = K$（工具变量与内生变量相关）

**$\mathbf{T=3}$ 时**: 等式变为对截面数据的2SLS——$(\mathbf{x}_{i2} - \mathbf{x}_{i1})$ 为 $(\mathbf{x}_{i3} - \mathbf{x}_{i2})$ 的工具变量。

**$\mathbf{T=2}$ 时**: $\Delta y_{i2} = \Delta \mathbf{x}_{i2}\boldsymbol{\beta} + \Delta u_{i2}$，$\mathbf{x}_{i1}$ 与 $\Delta u_{i2}$ 不相关 → $\mathbf{x}_{i1}$ 为 $\Delta \mathbf{x}_{i2}$ 的IV。但 $\mathbf{x}_{i1}$ 与 $(\mathbf{x}_{i2} - \mathbf{x}_{i1})$ 的相关可能很弱 → **大渐近方差**。

### 4.3 用滞后水平替代滞后差分

使用 $(\mathbf{x}_{i,t-1}, \mathbf{x}_{i,t-2})$ 而非 $\Delta \mathbf{x}_{i,t-1}$ 作为IV **不降低效率**——因为 $\Delta \mathbf{x}_{i,t-1}$ 是前两者的线性组合。

---

## 5. GMM工具变量矩阵 (Arellano-Bond型)

为了在每个时点使用所有可用的工具变量：

$$\mathbf{Z}_i = \begin{pmatrix}
\mathbf{X}_{i1}^0 & \mathbf{0} & \cdots & \mathbf{0} \\
\mathbf{0} & \mathbf{X}_{i2}^0 & \cdots & \mathbf{0} \\
\vdots & \vdots & \ddots & \vdots \\
\mathbf{0} & \mathbf{0} & \cdots & \mathbf{X}_{i,T-1}^0
\end{pmatrix}$$

其中 $\mathbf{X}_{i,t-1}^0 = (\mathbf{x}_{i1}, \mathbf{x}_{i2}, ..., \mathbf{x}_{i,t-1})$。

**关键特性**:
- $\mathbf{Z}_i$ 有 $T-1$ 行（对应 $T-1$ 个差分方程）
- 每行使用**不同的工具变量**——越往后时期，可用工具变量越多
- 这是经典Arellano-Bond差分GMM的结构

---

## 6. 同时含严格外生和序贯外生变量的模型

### 6.1 模型设定

$$y_{it} = \mathbf{z}_{it}\boldsymbol{\gamma} + \mathbf{w}_{it}\boldsymbol{\delta} + c_i + u_{it} \tag{11.18}$$

| 变量 | 外生性类型 | 条件 |
|------|-----------|------|
| $\mathbf{z}_{it}$ | **严格外生** | $E(\mathbf{z}_{is}' u_{it}) = 0$，所有 $s,t$ |
| $\mathbf{w}_{it}$ | **序贯外生** | $E(u_{it} \mid \mathbf{z}_i, \mathbf{w}_{it}, \mathbf{w}_{i,t-1},...,\mathbf{w}_{i1}) = 0$ |

### 6.2 差分与IV

差分后: $\Delta y_{it} = \Delta \mathbf{z}_{it}\boldsymbol{\gamma} + \Delta \mathbf{w}_{it}\boldsymbol{\delta} + \Delta u_{it}$

在时刻 $t$ 可用的IV: $(\mathbf{z}_i, \mathbf{w}_{i,t-1}, ..., \mathbf{w}_{i1})$——整个 $\mathbf{z}_i$ 向量（所有时期）+ 过去的 $\mathbf{w}$

用**混合2SLS**估计。

---

## 7. 反馈模型 (Feedback Models)

### 7.1 Example 11.2

**主方程**:

$$y_{it} = \mathbf{z}_{it}\boldsymbol{\gamma} + \delta w_{it} + c_i + u_{it} \tag{11.6}$$

**反馈方程**:

$$w_{it} = \mathbf{z}_{it}\boldsymbol{\psi} + \rho_1 y_{i,t-1} + \psi c_i + \varepsilon_{it} \tag{11.8}$$

**经济解释**: $y_{it}$ = 城市 $i$ 第 $t$ 年的人均安全套销量，$w_{it}$ = HIV感染率。

- 方程 (11.6): 安全套使用是否受HIV传播影响？
- 方程 (11.8): HIV传播是否受过去的安全套使用影响？

### 7.2 为什么严格外生性失败？

$$E(w_{i,t+1} \cdot u_{it}) = \rho_1 E(y_{it} \cdot u_{it}) = \rho_1 E(u_{it}^2) > 0$$

→ 当期的 $u_{it}$ 通过影响 $y_{it}$ → 影响 $w_{i,t+1}$ → $w$ 与未来 $u$ 相关。

**除非 $\rho_1 = 0$**（无反馈），否则严格外生性必然不成立。

---

## 8. 面板数据中的测量误差

### 8.1 模型

$$y_{it} = \beta x_{it}^* + c_i + u_{it}$$

$x_{it} = x_{it}^* + r_{it}$（$r_{it}$ 为测量误差）。

### 8.2 混合OLS的偏误

$$\text{plim } \hat{\beta}_{POLS} = \beta + \frac{\text{Cov}(x_{it}, c_i) - \beta\sigma_r^2}{\text{Var}(x_{it})}$$

**两个偏误来源**:
1. 遗漏变量偏误: $\text{Cov}(x_{it}, c_i)$
2. 衰减偏误: $-\beta\sigma_r^2$

若 $x_{it}$ 与 $c_i$ 正相关且 $\beta > 0$，两个偏误方向相反 → **可能相互抵消**。

### 8.3 差分加剧测量误差偏误 (Solon 1985)

对 $T=2$:

$$\text{plim } \hat{\beta}_{FD} = \beta - \frac{\sigma_r^2(1-\rho_r)}{\sigma_{x^*}^2(1-\rho_{x^*}) + \sigma_r^2(1-\rho_r)}$$

- 当 $\rho_{x^*} \to 1$（$x^*$ 高度持久），偏误 $\to -\beta$（**完全衰减**）
- 差分大幅减少 $x$ 的变异（分母中 $\text{Var}(\Delta x^*)$ 小） → **放大**信噪比问题

**一般方法**: (1) 变换消除 $c_i$，(2) 为变换后方程中的内生变量选择工具变量。

---

## 9. 个体特定斜率模型

### 9.1 随机趋势模型 (Random Trend Model)

$$y_{it} = c_i + g_i t + \mathbf{x}_{it}\boldsymbol{\beta} + u_{it}, \quad t = 1,...,T \tag{11.38}$$

每个个体有自己的时间趋势 $g_i t$。若 $y_{it}$ 是对数形式，$g_i$ 近似为**平均增长率**。

**估计**: 差分消除 $c_i$ → $\Delta y_{it} = g_i + \Delta \mathbf{x}_{it}\boldsymbol{\beta} + \Delta u_{it}$ ($t=2,...,T$)

→ 变成标准非观测效应模型（$g_i$ 为新"$c_i$"）。可用FE或FD估计 $\boldsymbol{\beta}$。**需要 $T \geq 3$**。

### 9.2 一般个体特定斜率模型

$$y_{it} = \mathbf{z}_{it}\mathbf{a}_i + \mathbf{x}_{it}\boldsymbol{\beta} + u_{it} \tag{11.43}$$

- $\mathbf{z}_{it}$: $1 \times J$（可与 $\mathbf{x}_{it}$ 有重叠）
- $\mathbf{a}_i$: $J \times 1$ 个体特定斜率（允许个体异质性与 $\mathbf{z}_{it}$ 交互）

**特例**:
- $J=1$, $\mathbf{z}_{it}=1$ → 标准FE模型 ($T \geq 2$)
- $J=2$, $\mathbf{z}_{it}=(1,t)$ → 随机趋势模型 ($T \geq 3$)
- $\mathbf{z}_{it} = (1, \text{prog}_{it})$ → 程序参与效应存在个体异质性

### 9.3 $\mathbf{M}_i$ 变换（扩展的组内变换）

令 $\mathbf{Z}_i$ 为 $T \times J$ 矩阵（第 $t$ 行为 $\mathbf{z}_{it}$），$\mathbf{X}_i$ 为 $T \times K$ 矩阵。

定义投影矩阵:

$$\mathbf{M}_i = \mathbf{I}_T - \mathbf{Z}_i(\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'$$

$\mathbf{M}_i \mathbf{Z}_i = \mathbf{0}$ → **$\mathbf{a}_i$ 被消除**。

$$\ddot{\mathbf{y}}_i = \mathbf{M}_i \mathbf{y}_i, \quad \ddot{\mathbf{X}}_i = \mathbf{M}_i \mathbf{X}_i, \quad \ddot{\mathbf{u}}_i = \mathbf{M}_i \mathbf{u}_i$$

$$\ddot{\mathbf{y}}_i = \ddot{\mathbf{X}}_i \boldsymbol{\beta} + \ddot{\mathbf{u}}_i$$

### 9.4 FE估计量与渐近方差

$$\hat{\boldsymbol{\beta}}_{FE} = \left( \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i \right)^{-1} \left( \sum_{i=1}^{N} \ddot{\mathbf{X}}_i' \ddot{\mathbf{y}}_i \right)$$

在 $\text{FE.3}'$ ($E(\mathbf{u}_i\mathbf{u}_i' | \mathbf{Z}_i, \mathbf{X}_i, \mathbf{a}_i) = \sigma_u^2 \mathbf{I}_T$) 下：

$$\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}_{FE}) = \hat{\sigma}_u^2 \left( \sum \ddot{\mathbf{X}}_i' \ddot{\mathbf{X}}_i \right)^{-1}, \quad \hat{\sigma}_u^2 = \frac{1}{N(T-1-J)-K} \sum_{i,t} \hat{u}_{it}^2 \tag{11.50}$$

### 9.5 平均偏效应 (APE) 的估计

关注 $\boldsymbol{\alpha} = E(\mathbf{a}_i)$（$\mathbf{z}_{it}$ 的平均偏效应）：

$$\mathbf{a}_i = (\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\boldsymbol{\beta}) - (\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'\mathbf{u}_i$$

$$\hat{\boldsymbol{\alpha}} = \frac{1}{N} \sum_{i=1}^{N} (\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\hat{\boldsymbol{\beta}}_{FE}) \tag{11.54}$$

- $\hat{\boldsymbol{\alpha}}$ 是 $\sqrt{N}$-一致且渐近正态的
- Chamberlain (1992a) 提供了更有效的联合估计 $\boldsymbol{\alpha}$ 和 $\boldsymbol{\beta}$ 的方法

---

## 10. Chamberlain 方法

### 10.1 解决什么问题？

标准FE的组内变换虽然允许 $c_i$ 与 $\mathbf{x}_{it}$ 任意相关，但有两个代价：

1. **所有时间不变变量被一并消除** — 无法估计 $\mathbf{z}_i$ 的系数
2. **无法对 $c_i$ 做推断** — $c_i$ 被消除后，无法估计 $E(c_i)$ 或检验 $c_i$ 与哪些 $\mathbf{x}$ 相关

Chamberlain 方法的替代思路：**与其消除 $c_i$，不如将其替换为 $\mathbf{x}_i$（所有时期）的线性投影**，从而在保留 $\boldsymbol{\beta}$ 一致性的同时获得更大的灵活性。

| | 组内变换（FE） | Chamberlain 方法 |
|---|---|---|
| 处理 $c_i$ 的方式 | **消除**（差分/去均值） | **替代**（用 $\mathbf{x}_i$ 的线性投影） |
| $\boldsymbol{\beta}$ 的一致性 | ✓ | ✓ |
| 时间不变变量的系数 | 无法估计 | 可通过 $\boldsymbol{\lambda}_t$ 的限制形式纳入 |
| 对 $c_i$ 的推断 | 通常不做 | 可估计 $E(c_i)$ 等矩 |
| 假设强度 | 仅需严格外生 | 额外要求 $E(\mathbf{X}_i' v_{it}) = \mathbf{0}$（若仅需正交条件，该条件自动成立） |
| 效率 | FE.3下有效 | 系统OLS非有效（$v_{it}$ 存在序列相关和异方差） |

### 10.2 核心思想

在FE环境下：$\mathbf{x}_{it}$ 仅包含时变变量，严格外生性成立，允许 $c_i$ 与 $\mathbf{x}_{it}$ 任意相关。

假设 $c_i$ 和 $\mathbf{x}_i$ 的所有元素具有有限二阶矩。**无论真实的 $E(c_i|\mathbf{x}_i)$ 是什么函数形式**，我们总是可以写：

$$c_i = \psi + \mathbf{x}_{i1}\boldsymbol{\lambda}_1 + \mathbf{x}_{i2}\boldsymbol{\lambda}_2 + \cdots + \mathbf{x}_{iT}\boldsymbol{\lambda}_T + a_i \tag{11.59}$$

- $\psi$ 为标量
- $\boldsymbol{\lambda}_1,...,\boldsymbol{\lambda}_T$ 为 $1 \times K$ 向量（每个时期一组参数）
- $a_i$: 投影误差，满足 $E(a_i)=0$，$\text{Cov}(a_i, \mathbf{x}_{it}) = \mathbf{0}$（对所有 $t$）

**关键**：(11.59) 不是一个"假设"——它是**总成立**的线性投影（只要二阶矩存在）。它仅假设 $c_i$ 的条件**期望**是线性的，而不是说 $c_i$ 本身一定是 $\mathbf{x}_i$ 的线性函数。

### 10.3 代入后的方程

将 (11.59) 代入 $y_{it} = \mathbf{x}_{it}\boldsymbol{\beta} + c_i + u_{it}$：

$$y_{it} = \psi + \mathbf{x}_{i1}\boldsymbol{\lambda}_1 + \cdots + \mathbf{x}_{it}(\boldsymbol{\beta} + \boldsymbol{\lambda}_t) + \cdots + \mathbf{x}_{iT}\boldsymbol{\lambda}_T + v_{it}$$

其中 $v_{it} = a_i + u_{it}$。

**正交条件**：在严格外生性FE.1下：

$$E(v_{it}) = 0, \quad E(\mathbf{X}_i' v_{it}) = \mathbf{0}, \quad t = 1,...,T \tag{11.61}$$

注意 (11.61) 说的是 $\mathbf{X}_i$（**所有时期**的 $\mathbf{x}$）与 $v_{it}$ 不相关——这来自两个事实的叠加：
- $a_i$ 与所有 $\mathbf{x}_{it}$ 不相关（投影误差的定义）
- $u_{it}$ 与所有 $\mathbf{x}_{is}$ 不相关（严格外生性 FE.1）

但需警惕：除非额外假设 $E(c_i|\mathbf{x}_i)$ 是线性的，否则 $E(v_{it} | \mathbf{x}_i) \neq 0$——也就是说，(11.61) 保证的是**正交性**（一阶矩），而非**条件均值独立**（零条件期望）。

### 10.4 系统方程

堆叠所有 $T$ 个时期，写成系统形式：

$$\begin{pmatrix} y_{i1} \\ y_{i2} \\ \vdots \\ y_{iT} \end{pmatrix} = \begin{pmatrix}
1 & \mathbf{X}_{i1} & \mathbf{X}_{i2} & \cdots & \mathbf{X}_{iT} & \mathbf{X}_{i1} \\
1 & \mathbf{X}_{i1} & \mathbf{X}_{i2} & \cdots & \mathbf{X}_{iT} & \mathbf{X}_{i2} \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
1 & \mathbf{X}_{i1} & \mathbf{X}_{i2} & \cdots & \mathbf{X}_{iT} & \mathbf{X}_{iT}
\end{pmatrix} \begin{pmatrix} \psi \\ \boldsymbol{\lambda}_1 \\ \boldsymbol{\lambda}_2 \\ \vdots \\ \boldsymbol{\lambda}_T \\ \boldsymbol{\beta} \end{pmatrix} + \mathbf{v}_i$$

即 $\mathbf{y}_i = \mathbf{W}_i \boldsymbol{\mu} + \mathbf{v}_i$，其中 $\mathbf{W}_i$ 为 $T \times (1 + TK + K)$ 设计矩阵，$\boldsymbol{\mu}$ 为 $(1 + TK + K) \times 1$ 参数向量。

注意设计矩阵的结构：对第 $t$ 行，$\boldsymbol{\beta}$ 总是与 $\mathbf{x}_{it}$ 相乘，但 $\boldsymbol{\lambda}_t$ 也同时与 $\mathbf{x}_{it}$ 相乘 → $\mathbf{x}_{it}$ 一共出现在两个位置：一次乘以 $\boldsymbol{\beta}$，一次乘以 $\boldsymbol{\lambda}_t$ → 其总系数为 $\boldsymbol{\beta} + \boldsymbol{\lambda}_t$。

### 10.5 估计与性质

**系统OLS**：由 $E(\mathbf{W}_i' \mathbf{v}_i) = \mathbf{0}$，系统OLS一致（$N \to \infty$，$T$ 固定）。秩条件要求：

$$\text{rank } E(\mathbf{W}_i' \mathbf{W}_i) = 1 + TK + K$$

**效率问题**：系统OLS虽然一致，但**不是有效的**，因为：
1. $E(\mathbf{v}_i \mathbf{v}_i' | \mathbf{x}_i) \neq \sigma_v^2 \mathbf{I}_T$——$v_{it}$ 通过 $a_i$ 存在序列相关
2. $E(\mathbf{v}_i \mathbf{v}_i' | \mathbf{x}_i) \neq E(\mathbf{v}_i \mathbf{v}_i')$——可能存在条件异方差

→ 更高效的估计需使用GLS/GMM。Chamberlain (1992a) 提出了有效矩方法（Efficient Method of Moments）来联合估计 $\boldsymbol{\beta}$ 和 $\boldsymbol{\alpha} = E(\mathbf{a}_i)$。

**与CRE（相关随机效应）的联系**：Chamberlain方法是Chamberlain相关随机效应Probit/Logit模型的线性版本，也是Mundlak方法和CRE方法的理论基础。

---

## 11. Hausman-Taylor 模型

### 11.1 问题

当主要关心的变量是**时间不变**的（如教育年限），但 $c_i$ 可能与部分变量相关 → FE消除了所有时间不变变量。

### 11.2 模型与假设

$$y_{it} = \mathbf{z}_i\boldsymbol{\gamma} + \mathbf{x}_{it}\boldsymbol{\beta} + c_i + u_{it} \tag{11.65}$$

将 $\mathbf{z}_i$ 和 $\mathbf{x}_{it}$ 分别划分为两组：

$$\mathbf{z}_i = (\mathbf{z}_{i1}, \mathbf{z}_{i2}), \quad \mathbf{x}_{it} = (\mathbf{x}_{it1}, \mathbf{x}_{it2})$$

| 变量 | 维度 | 假设 |
|------|------|------|
| $\mathbf{z}_{i1}$ | $1 \times J_1$ | $E(\mathbf{z}_{i1}' c_i) = 0$（外生时间不变变量） |
| $\mathbf{z}_{i2}$ | $1 \times J_2$ | 可能与 $c_i$ 相关（内生时间不变变量） |
| $\mathbf{x}_{it1}$ | $1 \times K_1$ | $E(\mathbf{x}_{it1}' c_i) = 0, \forall t$（外生时变变量） |
| $\mathbf{x}_{it2}$ | $1 \times K_2$ | 可能与 $c_i$ 相关（内生时变变量） |

**总外生性条件**: RE.1a: $E(u_{it} | \mathbf{z}_i, \mathbf{x}_{i1},...,\mathbf{x}_{iT}, c_i) = 0, \; t = 1,...,T$ (11.66)

### 11.3 识别与工具变量

**堆叠方程**: $\mathbf{y}_i = \mathbf{z}_i\boldsymbol{\gamma} + \mathbf{X}_i\boldsymbol{\beta} + \mathbf{V}_i, \quad \mathbf{V}_i = c_i\mathbf{j}_T + \mathbf{u}_i$ (11.69)

**工具变量集合**:

1. $\mathbf{Q}_T\mathbf{X}_i$（组内变换后的 $\mathbf{X}_i$）——因为 $\mathbf{Q}_T\mathbf{V}_i = \mathbf{Q}_T\mathbf{u}_i$
2. $\mathbf{j}_T \otimes (\mathbf{z}_{i1}, \bar{\mathbf{x}}_{i1}^0)$——来自 (11.68) 的正交条件

**时刻 $t$ 的IV向量**: $(\ddot{\mathbf{x}}_{it}, \mathbf{z}_{i1}, \bar{\mathbf{x}}_{i1}^0)$

### 11.4 估计效率

在Hausman-Taylor假设下，**广义IV (GIV) 估计量**是有效的GMM估计量。所有来自拟去均值数据的混合2SLS统计量都渐近有效。

### 11.5 应用示例 (HT 1981)

估计教育回报：
- **内生时间不变变量** ($\mathbf{z}_{i2}$): 教育年限
- **外生时间不变变量** ($\mathbf{z}_{i1}$): 种族、工会身份
- **外生时变变量** ($\mathbf{x}_{it1}$): 经验、健康不良指标、前一年失业指标
- **内生时变变量** ($\mathbf{x}_{it2}$): （此处可能为空）

---

## 12. 匹配对与聚类样本

### 12.1 匹配对 (兄弟姐妹数据)

$$y_{i1} = \mathbf{x}_{i1}\boldsymbol{\beta} + f_i + u_{i1}, \quad y_{i2} = \mathbf{x}_{i2}\boldsymbol{\beta} + f_i + u_{i2}$$

$f_i$ = 不可观测的**家庭效应**。

- **RE**: $f_i$ 与 $(\mathbf{x}_{i1}, \mathbf{x}_{i2})$ 不相关 → 混合OLS一致
- **FE**: 允许任意相关 → 跨兄弟姐妹差分消除 $f_i$

**严格外生性**: 每个兄弟姐妹方程中的 $u_{is}$ 与**两个**方程中的解释变量都不相关

### 12.2 聚类样本

$$y_{is} = \mathbf{x}_{is}\boldsymbol{\beta} + c_i + u_{is}$$

$i$ = 聚类（学校/公司/家庭），$s$ = 聚类内个体。

**一般模型**: $\mathbf{y}_i = \mathbf{X}_i\boldsymbol{\beta} + c_i\mathbf{j}_{G_i} + \mathbf{u}_i$ (11.74)

$G_i$ = 聚类 $i$ 中的个体数（**允许不同聚类大小不同**）。

- **FE**: 聚类内去均值消除 $c_i$。聚类大小不同**不造成问题**——去均值在每个聚类内单独进行。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $E(u_{it} \mid \mathbf{x}_{it},...,\mathbf{x}_{i1}, c_i) = 0$ | 序贯矩条件 (11.2) |
| $\Delta y_{it} = \Delta \mathbf{x}_{it}\boldsymbol{\beta} + \Delta u_{it}$ | FD方程 (11.14) |
| $\text{plim } \hat{\boldsymbol{\beta}}_{FE} = \boldsymbol{\beta} + O(T^{-1})$ | FE偏误量级 |
| $\text{plim } \hat{\boldsymbol{\beta}}_{FD} = \boldsymbol{\beta} + O(1)$ | FD偏误量级 |
| $c_i = \psi + \sum_{t=1}^{T} \mathbf{x}_{it}\boldsymbol{\lambda}_t + a_i$ | Chamberlain投影 (11.59) |
| $\mathbf{y}_i = \mathbf{W}_i\boldsymbol{\mu} + \mathbf{v}_i$ | Chamberlain系统 (11.63) |
| $\mathbf{M}_i = \mathbf{I}_T - \mathbf{Z}_i(\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'$ | 个体斜率消除矩阵 |
| $\hat{\boldsymbol{\alpha}} = \frac{1}{N}\sum (\mathbf{Z}_i'\mathbf{Z}_i)^{-1}\mathbf{Z}_i'(\mathbf{y}_i - \mathbf{X}_i\hat{\boldsymbol{\beta}}_{FE})$ | APE估计 (11.54) |
