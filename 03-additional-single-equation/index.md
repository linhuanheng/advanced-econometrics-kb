---
title: 单方程扩展专题
type: framework
domain: econometrics
source: Chapter IV — Additional Single-Equation Topics
created: 2026-05-19
updated: 2026-05-24
---

# Ch4: 单方程扩展专题

## 概述

本章涵盖OLS和2SLS估计中的高级专题：生成回归元/工具变量的推断、控制函数方法、设定检验（内生性、过度识别、函数形式、异方差）、相关随机系数模型及双重差分估计。

---

## 1. 生成回归元与工具变量

### 1.1 OLS with Generated Regressors (Pagan, 1984)

**模型**: $y = \beta_0 + \beta_1 x_1 + \cdots + \beta_K x_K + \gamma q + u$ (6.1)

其中 $q = f(\mathbf{w}, \boldsymbol{\delta})$，$\mathbf{w}$ 为可观测变量向量，$\hat{q} = f(\mathbf{w}, \hat{\boldsymbol{\delta}})$。

**第二阶段OLS**: $y_i$ 对 $1, x_{i1}, ..., x_{iK}, \hat{q}_i$ (6.2)

**一致性**: 只要 $u$ 与 $(x_1,...,x_K, q)$ 不相关即可。

**推断问题**: (6.2) 的标准误和t统计量**一般无效**——忽略了 $\hat{\boldsymbol{\delta}}$ 的抽样变异。

**忽略 $\hat{\boldsymbol{\delta}}$ 抽样变异 的条件**:

$$E[\nabla_{\boldsymbol{\delta}} f(\mathbf{w}, \boldsymbol{\delta})' u] = \mathbf{0} \tag{6.3}$$

$$\gamma = 0 \tag{6.4}$$

条件 (6.3) 由 $E(u | \mathbf{x}, \mathbf{w}) = 0$ 保证。

### 1.2 2SLS with Generated Instruments

$y = \mathbf{x}\boldsymbol{\beta} + u$, $E(\mathbf{z}'u) = \mathbf{0}$，但 $\mathbf{z} = g(\mathbf{w}, \boldsymbol{\lambda})$。

生成工具变量: $\hat{\mathbf{z}}_i = g(\mathbf{w}_i, \hat{\boldsymbol{\lambda}})$

**一致性充分条件**: (1) $\hat{\boldsymbol{\lambda}}$ 为 $\sqrt{N}$-一致的，(2) $E[\nabla_{\boldsymbol{\lambda}} g(\mathbf{w}, \boldsymbol{\lambda})' u] = \mathbf{0}$ (6.8)

条件 (6.8) 在 $E(u | \mathbf{w}) = 0$ 下成立。$\sqrt{N}$-渐近分布与使用真实 $\boldsymbol{\lambda}$ 一致。

### 1.3 同时含生成回归元和生成工具变量

$y = \mathbf{x}\boldsymbol{\beta} + \gamma f(\mathbf{w}, \boldsymbol{\delta}) + u$, $E(u | \mathbf{z}, \mathbf{w}) = 0$

若 $\gamma = 0$，2SLS估计 $(\boldsymbol{\beta}, \gamma)'$（IV为 $(\mathbf{z}_i, f_i)$）的渐近分布不受 $\hat{\boldsymbol{\delta}}$ 影响——可用通常的2SLS t统计量（或其异方差稳健版本）检验 $H_0: \gamma = 0$。

---

## 2. 控制函数方法 (Control Function Approach)

### 2.1 核心思想

通过加入约简型残差来"切断"内生变量与不可观测因素之间的相关性。

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + u_1, \quad E(\mathbf{z}'u_1) = \mathbf{0}$$

**约简型**: $y_2 = \mathbf{z}\boldsymbol{\pi}_2 + v_2$, $E(\mathbf{z}'v_2) = 0$

**内生性的充要条件**: $u_1$ 与 $v_2$ 相关。

### 2.2 控制函数的推导

将 $u_1$ 对 $v_2$ 线性投影:

$$u_1 = \rho_1 v_2 + e_1, \quad \rho_1 = \frac{E(u_1 v_2)}{E(v_2^2)}$$

其中 $E(v_2 e_1) = 0$, $E(\mathbf{z}'e_1) = \mathbf{0}$。

代入结构方程:

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + \rho_1 v_2 + e_1$$

$e_1$ 与 $\mathbf{z}_1, y_2, v_2$ 都不相关 → OLS一致。

### 2.3 CF估计方程

用 $\hat{v}_2 = y_2 - \mathbf{z}\hat{\boldsymbol{\pi}}_2$ 替代 $v_2$:

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + \rho_1 \hat{v}_2 + \text{error} \tag{6.15}$$

OLS估计即为**控制函数 (CF) 估计量**。

### 2.4 CF vs 2SLS

**线性模型中**: CF与2SLS**完全等价**——系数估计值相同。

**CF的两个优势**:
1. 提供简单的内生性检验: $H_0: \rho_1 = 0$（对 $\hat{v}_2$ 的t检验）
2. 可适应**非线性模型**（含 $y_2^2$ 等）而2SLS无法直接推广

**CF的排除约束**: 需要 $\mathbf{z}$ 中至少有一个元素不在 $\mathbf{z}_1$ 中（与2SLS相同）。若 $\mathbf{z} = \mathbf{z}_1$，存在完全共线性。

### 2.5 CF与2SLS的分歧：含非线性函数的内生变量

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + \gamma_1 y_2^2 + u_1, \quad E(u_1 | \mathbf{z}) = 0$$

**2SLS**: 解释变量 $(\mathbf{z}_1, y_2, y_2^2)$，IV $(\mathbf{z}_1, \mathbf{z}_2, \mathbf{z}_2^2)$

**CF**: 假设 $E(u_1 | \mathbf{z}, y_2) = E(u_1 | \mathbf{z}, v_2) = E(u_1 | v_2) = \rho_1 v_2$

OLS回归 $y_1$ 对 $\mathbf{z}_1, y_2, y_2^2, \hat{v}_2$ → CF与2SLS**不同**。

**权衡**: CF更精细（需线性条件期望假设），但更有效；2SLS更稳健但效率较低。

### 2.6 多内生变量的CF

$y_2$ 为 $1 \times G_1$ 内生变量向量:

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \mathbf{y}_2\boldsymbol{\alpha}_1 + u_1$$

约简型: $\mathbf{y}_2 = \mathbf{z}\boldsymbol{\Pi}_2 + \mathbf{v}_2$，$\hat{\mathbf{v}}_2$ 为 $1 \times G_1$ 残差向量

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \mathbf{y}_2\boldsymbol{\alpha}_1 + \hat{\mathbf{v}}_2\boldsymbol{\rho}_1 + \text{error}$$

$H_0: \boldsymbol{\rho}_1 = \mathbf{0}$ → 标准F检验（$G_1$ 个约束）。或LM检验: $\hat{u}_1$（$y_1$ 对 $\mathbf{z}_1, \mathbf{y}_2$ 的OLS残差）对 $\mathbf{z}_1, \mathbf{y}_2, \hat{\mathbf{v}}_2$ 回归 → $NR_u^2 \xrightarrow{a} \chi^2_{G_1}$。

---

## 3. 设定检验

### 3.1 内生性检验 (Durbin-Wu-Hausman, DWH)

**原始DWH检验**（基于 $\hat{\boldsymbol{\beta}}_{2SLS} - \hat{\boldsymbol{\beta}}_{OLS}$ 的差）:

在 $H_0: E(\mathbf{x}'u) = 0$（无内生性）和同方差下:

$$\text{Avar}[\sqrt{N}(\hat{\boldsymbol{\beta}}_{2SLS} - \hat{\boldsymbol{\beta}}_{OLS})] = \sigma^2[E(\mathbf{x}^{*'}\mathbf{x}^*)]^{-1} - \sigma^2E(\mathbf{x}'\mathbf{x})^{-1}$$

其中 $\mathbf{x}^* = \mathbf{z}\boldsymbol{\pi}$（第一阶段拟合值）。

$$DWH = (\hat{\boldsymbol{\beta}}_{2SLS} - \hat{\boldsymbol{\beta}}_{OLS})' [(\hat{\mathbf{x}}'\hat{\mathbf{x}})^{-1} - (\mathbf{x}'\mathbf{x})^{-1}]^{-1} (\hat{\boldsymbol{\beta}}_{2SLS} - \hat{\boldsymbol{\beta}}_{OLS}) / \hat{\sigma}^2 \xrightarrow{d} \chi^2$$

**回归型Hausman检验** (CF方法——更方便):

1. $y_2$ 对 $\mathbf{z}$ 回归 → $\hat{v}_2$
2. $y_1$ 对 $\mathbf{z}_1, y_2, \hat{v}_2$ 回归
3. $H_0: \rho_1 = 0$ → t检验（或异方差稳健t检验）

**Hausman t统计量**（同方差下）：

$$\frac{\hat{\alpha}_{1,2SLS} - \hat{\alpha}_{1,OLS}}{\sqrt{[se(\hat{\alpha}_{1,2SLS})]^2 - [se(\hat{\alpha}_{1,OLS})]^2}} \xrightarrow{d} N(0,1)$$

对多个潜在内生变量的扩展: 使用F检验或LM检验。

**Example 6.2 (Card 1995 教育内生性检验)**:

$$\log(\text{wage}) = \alpha_1 \text{educ} + \alpha_2 \text{black} \cdot \text{educ} + \mathbf{z}_1\boldsymbol{\delta}_1 + u_1 \tag{6.20}$$

IV: motheduc, fatheduc, huseduc。$\hat{v}_2$ 上的t统计量检验 $H_0: educ$ 外生。

### 3.2 过度识别约束检验 (Sargan, 1958)

$$y_1 = \mathbf{z}_1\boldsymbol{\delta}_1 + \mathbf{y}_2\boldsymbol{\alpha}_1 + u_1 \tag{6.30}$$

$\mathbf{z} = (\mathbf{z}_1, \mathbf{z}_2)$，$\mathbf{z}_2$ 为 $1 \times L_2$，$L_2 > G_1$ → 过度识别。

**步骤**:

1. 使用**所有** $\mathbf{z}$ 进行2SLS → 得残差 $\hat{u}_1$
2. $\hat{u}_1$ 对 $\mathbf{z}$ 回归 → 得 $R_u^2$
3. $NR_u^2 \xrightarrow{a} \chi^2_{Q_1}$，$Q_1 = L_2 - G_1$（过度识别约束数）

$H_0$: 所有 $\mathbf{z}_2$ 中的IV与 $u_1$ 不相关（即**额外IV有效**）。

**Example 6.3**: motheduc, fatheduc, huseduc 三个IV → $Q_1 = 3-1 = 2$ 个过度识别约束。

测试 $\hat{u}_1$ 对 $1, \text{exper}, \text{exper}^2, \text{motheduc}, \text{fatheduc}, \text{huseduc}$ 的回归。

**注意**: 拒绝 $H_0$ 时不指示**哪个**IV有问题。

### 3.3 函数形式检验 (2SLS的RESET)

$y = \mathbf{x}\boldsymbol{\beta} + u$, $E(u | \mathbf{x}) = 0$ (6.23, 6.24)

在 $E(u|\mathbf{x}) = 0$ 下，$\mathbf{x}$ 的任何函数与 $u$ 无关 → $(\mathbf{x}\boldsymbol{\beta})^p$ 与 $u$ 无关。

检验 $\hat{u}_i$（OLS或2SLS残差）是否与 $\hat{y}_i^2, \hat{y}_i^3, \hat{y}_i^4$ 相关:

- (a) 加入这些项到回归中做F检验 ($F_{3, N-K-3}$)
- (b) LM: $\hat{u}$ 对 $(\mathbf{x}, \hat{y}^2, \hat{y}^3, \hat{y}^4)$ 回归 → $NR_u^2 \sim \chi^2_3$

### 3.4 2SLS框架下的异方差检验

$y = \beta_0 + \mathbf{x}\boldsymbol{\beta} + u$, $E(u | \mathbf{x}) = 0$

$H_0: E(u^2 | \mathbf{x}) = \sigma^2$; 备择: $E(u^2 | \mathbf{x})$ 依赖 $\mathbf{x}$

$$\hat{u}_i^2 = \delta_0 + \mathbf{h}_i\boldsymbol{\delta} + v_i, \quad \mathbf{h}_i \equiv h(\mathbf{x}_i)$$

$H_0: \boldsymbol{\delta} = \mathbf{0}$

需homokurtosis条件: $E(u_i^4 | \mathbf{x}_i) = \kappa^2$ (6.37)

$N R_c^2 \xrightarrow{a} \chi^2_Q$（$Q$ 为 $\mathbf{h}_i$ 的维度）

**三种常用设定**:
1. **Breusch-Pagan**: $\mathbf{h}_i = \mathbf{x}_i$, $Q=K$
2. **White**: $\mathbf{h}_i$ 含 $\mathbf{x}_i$ 的所有非恒定、唯一元素（水平+平方+交叉）
3. **简化版**: $\mathbf{h}_i = (\hat{y}_i, \hat{y}_i^2)$，仅需 $Q=2$

---

## 4. 相关随机系数模型 (Correlated Random Coefficient, CRC)

### 4.1 CRC模型

$$y_1 = \eta_1 + \mathbf{z}_1\boldsymbol{\delta}_1 + a_1 y_2 + u_1 \tag{6.41}$$

$a_1$ 是 $y_2$ 上的"随机系数"——可能与 $y_2$ 相关 (Heckman & Vytlacil 1998)。

令 $a_1 = \alpha_1 + v_1$（$\alpha_1 = E(a_1)$ 为目标参数——平均偏效应）:

$$y_1 = \eta_1 + \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + v_1 y_2 + u_1 \tag{6.42}$$

问题: $e_1 = v_1 y_2 + u_1$ 不一定与 $\mathbf{z}$ 无关——即使 $E(u_1|\mathbf{z}) = E(v_1|\mathbf{z}) = 0$，$v_1 y_2$ 可能与 $\mathbf{z}$ 相关。

### 4.2 IV一致性的关键条件

**Wooldridge (2003b)**: IV一致性的**充分条件**为:

$$\text{Cov}(v_1, y_2 | \mathbf{z}) = \text{Cov}(v_1, y_2) \tag{6.44}$$

即 $v_1$ 与 $y_2$ 的条件协方差不依赖于 $\mathbf{z}$——等于无条件协方差。

**何时成立？** 若 $(v_1, v_2)$ 独立于 $\mathbf{z}$（其中 $y_2 = m_2(\mathbf{z}) + v_2$），则 $\text{Cov}(v_1, y_2 | \mathbf{z}) = \text{Cov}(v_1, v_2 | \mathbf{z}) = \text{Cov}(v_1, v_2)$（常数）→ 成立。

$v_1$ 不可观测 → 不可直接检验。

### 4.3 含交互项的CRC扩展

允许 $\mathbf{z}_1$ 与 $y_2$ 交互:

$$y_1 = \eta_1 + \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + (\mathbf{z}_1 - \boldsymbol{\mu}_1) y_2 \boldsymbol{\gamma}_1 + v_1 y_2 + u_1$$

其中 $\boldsymbol{\mu}_1 = E(\mathbf{z}_1)$，确保 $\alpha_1$ 为平均偏效应。

用样本均值 $\bar{\mathbf{z}}_1$ 替代 $\boldsymbol{\mu}_1$，IV估计。工具变量选择: $[1, \mathbf{z}_{i1}, (\mathbf{z}_{i1} - \bar{\mathbf{z}}_1)\hat{y}_{i2}]$（$\hat{y}_{i2} = \mathbf{z}_i\hat{\boldsymbol{\pi}}_2$）。

### 4.4 Garen (1984) 的控制函数方法

假设 $(u_1, v_1, v_2)$ 独立于 $\mathbf{z}$ 且联合正态（比实际需要的强）：

$$E(u_1 | \mathbf{z}, v_2) = \rho_1 v_2, \quad E(v_1 | \mathbf{z}, v_2) = \theta_1 v_2$$

则:

$$E(y_1 | \mathbf{z}, y_2) = \eta_1 + \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + \theta_1 v_2 y_2 + \rho_1 v_2$$

**Garen的两步法**:
1. $y_2$ 对 $\mathbf{z}$ 回归 → $\hat{v}_2$
2. OLS回归 $y_1$ 对 $1, \mathbf{z}_1, y_2, \hat{v}_2, \hat{v}_2 y_2$

$\hat{v}_2 y_2$ 这一项捕捉了随机斜率 $v_1$ 的内生性。

**外生性检验**: $H_0: \theta_1 = 0, \rho_1 = 0$（F检验，2个约束）

---

## 5. 双重差分 (Difference-in-Differences, DID)

### 5.1 混合横截面的DID

**设定**: 两个时期（$t=1,2$），两组（A=对照, B=处理）。

$$y = \beta_0 + \beta_1 dB + \delta_0 d2 + \delta_1 d2 \cdot dB + u$$

| 符号 | 含义 |
|------|------|
| $dB$ | 处理组虚拟变量（=1 若为B组） |
| $d2$ | 第二期虚拟变量（=1 若为 $t=2$） |
| $\delta_1$ | **DID估计量**——关心的因果效应 |

### 5.2 DID估计量的期望

$$\bar{y}_{A,1} = \beta_0, \quad \bar{y}_{A,2} = \beta_0 + \delta_0$$
$$\bar{y}_{B,1} = \beta_0 + \beta_1, \quad \bar{y}_{B,2} = \beta_0 + \beta_1 + \delta_0 + \delta_1$$

$$\hat{\delta}_1 = (\bar{y}_{B,2} - \bar{y}_{B,1}) - (\bar{y}_{A,2} - \bar{y}_{A,1})$$

### 5.3 为什么DID优于其他方法

| 方法 | 公式 | 问题 |
|------|------|------|
| 仅处理组前后对比 | $\bar{y}_{B,2} - \bar{y}_{B,1}$ | 均值可能因政策无关原因变化 |
| 仅第二期横截面对比 | $\bar{y}_{B,2} - \bar{y}_{A,2}$ | 处理组和对照组间存在系统性差异 |
| **DID** | $(\bar{y}_{B,2} - \bar{y}_{B,1}) - (\bar{y}_{A,2} - \bar{y}_{A,1})$ | 同时控制组别效应和时间效应 |

### 5.4 关键假设

- **平行趋势 (Parallel Trends)**: 无处理时，两组的变化趋势相同
- 政策变化不与影响 $y$ 的其他因素（隐藏在 $u$ 中）系统相关

### 5.5 实例：工伤保险金 (Meyer, Viscusi, & Durbin 1995)

研究Kentucky 1980年提高工伤赔偿上限对受伤工人停工时间的影响。

$$\log(\text{durat}) = 1.126 + 0.0077 \text{afchnge} + 0.256 \text{highearn} + 0.191 \text{afchnge} \cdot \text{highearn}$$

$\hat{\delta}_1 = 0.191$ ($t=2.77$) → 高收入工人的平均工伤停工时间增加了约19%。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $E[\nabla_{\delta}f(\mathbf{w},\delta)'u] = 0$ | 忽略生成回归元抽样变异条件 (6.3) |
| $u_1 = \rho_1 v_2 + e_1$ | CF线性投影分解 |
| $y_1 = \mathbf{z}_1\delta_1 + \alpha_1 y_2 + \rho_1 \hat{v}_2 + error$ | CF估计方程 (6.15) |
| $DWH = (\hat{\beta}_{2SLS}-\hat{\beta}_{OLS})'[...]^{-1}(...)/\hat{\sigma}^2$ | DWH内生性检验 |
| $NR_u^2 \sim \chi^2_{L_2-G_1}$ | Sargan过度识别检验 |
| $NR_c^2 \sim \chi^2_Q$ | 异方差检验 (2SLS框架) |
| $y_1 = \eta_1 + \mathbf{z}_1\delta_1 + \alpha_1 y_2 + v_1 y_2 + u_1$ | CRC模型 (6.42) |
| $\text{Cov}(v_1, y_2 \mid \mathbf{z}) = \text{Cov}(v_1, y_2)$ | CRC下IV一致性条件 (6.44) |
| $y_1 \sim 1, \mathbf{z}_1, y_2, \hat{v}_2, \hat{v}_2 y_2$ | Garen CF回归 |
| $\hat{\delta}_1 = (\bar{y}_{B,2}-\bar{y}_{B,1}) - (\bar{y}_{A,2}-\bar{y}_{A,1})$ | DID估计量 |
