---
title: 工具变量估计与两阶段最小二乘
type: framework
domain: econometrics
source: Chapter III — Instrumental Variables Estimation of Single-Equation Linear Models
created: 2026-05-19
updated: 2026-05-24
---

# Ch3: 工具变量估计与2SLS

## 核心问题

当解释变量与误差项相关（内生性），OLS估计量不一致。工具变量（IV）和两阶段最小二乘（2SLS）提供了内生性问题的通用解决方案。

内生性来源：遗漏变量（如能力）、测量误差（CEV下的衰减偏误）、联立性。

---

## 1. IV的基本原理

### 1.1 模型与假设

$$y = \beta_0 + \beta_1 x_1 + \cdots + \beta_K x_K + u \tag{5.1}$$

- $x_1, ..., x_{K-1}$ 为外生变量，$x_K$ 为**潜在内生变量**
- $E(u) = 0$, $\text{Cov}(x_j, u) = 0$ 对 $j = 1,...,K-1$（外生变量与误差无关）

### 1.2 工具变量的两个条件

工具变量 $z_1$ 需满足：

**条件一（外生性）**: $\text{Cov}(z_1, u) = 0$ (5.3)

$z_1$ 与结构误差不相关——排除约束（exclusion restriction）。

**条件二（相关性）**: 在约简型中 $\theta_1 \neq 0$

$$x_K = \delta_0 + \delta_1 x_1 + \cdots + \delta_{K-1} x_{K-1} + \theta_1 z_1 + r_K \tag{5.4}$$

$z_1$ 必须与内生变量 $x_K$ 相关。$\theta_1 \neq 0$ 是**关键识别条件**。

### 1.3 约简型与结构参数

$x_K$ 的约简型: 将内生变量写为所有外生变量的线性投影。

$y$ 的约简型: 将 (5.4) 代入 (5.1)：

$$y = \alpha_0 + \alpha_1 x_1 + \cdots + \alpha_{K-1} x_{K-1} + \lambda_1 z_1 + v$$

其中 $\alpha_j = \beta_j + \beta_K \delta_j$, $\lambda_1 = \beta_K \theta_1$。

**例子**: $y$ = 工人平均生产率, $z_1$ = 培训补贴（二元）, $x_K$ = 人均培训小时数。$\lambda_1$ 是"简化型效应"，$\beta_K$（培训小时数的结构效应）更有政策价值。

### 1.4 识别

将模型写为矩阵形式: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{u}$ (5.7)

外生变量向量: $\mathbf{Z} = (1, x_1, ..., x_{K-1}, z_1)$（$1 \times K$）

**总体矩条件**: $E(\mathbf{Z}'u) = 0$ (5.8)

代入 $u = y - \mathbf{X}\boldsymbol{\beta}$ 并取期望：

$$[E(\mathbf{Z}'\mathbf{X})]\boldsymbol{\beta} = E(\mathbf{Z}'y) \tag{5.9}$$

$E(\mathbf{Z}'\mathbf{X})$ 是 $K \times K$ 矩阵。唯一解存在的充要条件是：

$$\text{rank } E(\mathbf{Z}'\mathbf{X}) = K \tag{5.10}$$

**识别出的结构参数**:

$$\boldsymbol{\beta} = [E(\mathbf{Z}'\mathbf{X})]^{-1} E(\mathbf{Z}'y) \tag{5.11}$$

### 1.5 IV估计量

$$\hat{\boldsymbol{\beta}}_{IV} = \left(\frac{1}{N}\sum \mathbf{Z}_i'\mathbf{X}_i\right)^{-1} \left(\frac{1}{N}\sum \mathbf{Z}_i'\mathbf{y}_i\right) = (\mathbf{Z}'\mathbf{X})^{-1}\mathbf{Z}'\mathbf{Y}$$

$\mathbf{Z}$ 和 $\mathbf{X}$ 为 $N \times K$ 矩阵。

### 1.6 代理变量 vs 工具变量

| | 代理变量 (Proxy) | 工具变量 (IV) |
|------|------|------|
| 与遗漏变量的关系 | 应**高度相关** | 应**不相关** |
| 在结构模型中的角色 | 冗余 | 冗余 |
| 目的 | 替代不可观测变量 | 提供外生变异 |

### 1.7 实例：教育回报率的IV估计

**Example 5.1 (工资方程)**:

$$\log(\text{wage}) = \beta_0 + \beta_1 \text{exper} + \beta_2 \text{exper}^2 + \beta_3 \text{educ} + u \tag{5.12}$$

内生性来源: (1) 遗漏能力, (2) 教育质量差异, (3) 家庭背景。

**母亲教育 (motheduc) 作为IV的有效性条件**:
1. $\text{Cov}(\text{motheduc}, u) = 0$（母亲的受教育程度与能力无关）
2. $\theta_1 \neq 0$（母亲教育预测子女教育——可检验）

**Angrist & Krueger (1991)**: 使用**出生季度**作为教育的IV——(5.3)可能成立但(5.5)存疑（需检验第一阶段F统计量）。

**Card (1995)**: 使用**是否在四年制大学附近长大**作为教育IV。OLS估计回报率为7.5%，IV估计为13.2%。

**自然实验与项目评估**: 当个体被随机选为项目**合格者**时，"合格性"可作为"实际参与"的IV——合格性与参与相关，且合格性通常外生。

---

## 2. 两阶段最小二乘 (2SLS)

### 2.1 多工具变量的情况

当有多个工具变量 $z_1, ..., z_M$（$M \geq 1$）：

$$\text{Cov}(z_h, u) = 0, \quad h = 1,2,...,M \tag{5.13}$$

外生变量向量: $\mathbf{z} = (1, x_1, ..., x_{K-1}, z_1, ..., z_M)$，$1 \times L$，$L = K + M$。

### 2.2 2SLS的核心思想

在所有 $\mathbf{z}$ 的线性组合中，2SLS选择与 $x_K$ **最相关的**作为工具变量：即 $x_K$ 对 $\mathbf{z}$ 回归的拟合值。

### 2.3 第一阶段 (First Stage)

$$x_{iK} = \delta_0 + \delta_1 x_{i1} + \cdots + \delta_{K-1} x_{i,K-1} + \theta_1 z_{i1} + \cdots + \theta_M z_{iM} + r_{iK}$$

$$\hat{x}_{iK} = \hat{\delta}_0 + \hat{\delta}_1 x_{i1} + \cdots + \hat{\delta}_{K-1} x_{i,K-1} + \hat{\theta}_1 z_{i1} + \cdots + \hat{\theta}_M z_{iM}$$

### 2.4 第二阶段 (Second Stage)

$$y_i \text{ 对 } 1, x_{i1}, ..., x_{i,K-1}, \hat{x}_{iK} \text{ 进行OLS回归} \tag{5.19}$$

### 2.5 2SLS的矩阵表述

定义 $\hat{\mathbf{X}}_i = (1, x_{i1}, ..., x_{i,K-1}, \hat{x}_{iK})$。

$$\hat{\boldsymbol{\beta}}_{2SLS} = \left(\sum_{i=1}^{N} \hat{\mathbf{X}}_i'\hat{\mathbf{X}}_i\right)^{-1} \left(\sum_{i=1}^{N} \hat{\mathbf{X}}_i' y_i\right) = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}\hat{\mathbf{X}}'\mathbf{Y} \tag{5.17}$$

**等价表达**:

$$\hat{\boldsymbol{\beta}}_{2SLS} = (\hat{\mathbf{X}}'\mathbf{X})^{-1}\hat{\mathbf{X}}'\mathbf{Y} = (\mathbf{X}'\mathbf{P}_Z\mathbf{X})^{-1}\mathbf{X}'\mathbf{P}_Z\mathbf{y}$$

其中 $\mathbf{P}_Z = \mathbf{Z}(\mathbf{Z}'\mathbf{Z})^{-1}\mathbf{Z}'$，$\hat{\mathbf{X}} = \mathbf{P}_Z\mathbf{X}$。

**关键恒等式**: $\hat{\mathbf{X}}'\hat{\mathbf{X}} = \hat{\mathbf{X}}'\mathbf{X}$（因为 $\mathbf{P}_Z$ 幂等且对称）

### 2.6 识别状态

| 状态 | 条件 | 含义 |
|------|------|------|
| 恰好识别 | $L = K$ ($M = 1$) | IV唯一确定 |
| **过度识别** | $L > K$ ($M > 1$) | $M-1$ 个过度识别约束 |
| 不可识别 | $L < K$ | 参数无法识别 |

---

## 3. 2SLS的一般理论

### 3.1 一致性

**假设 2SLS.1** (正交条件): 存在 $1 \times L$ 向量 $\mathbf{z}$，$E(\mathbf{z}'u) = 0$。$\mathbf{z}$ 包含 $\mathbf{x}$ 的所有外生元素。

**假设 2SLS.2** (秩条件):
- (a) $\text{rank } E(\mathbf{z}'\mathbf{z}) = L$
- (b) $\text{rank } E(\mathbf{z}'\mathbf{x}) = K$

(b) 是关键的**秩条件**——$\mathbf{z}$ 与 $\mathbf{x}$ 充分线性相关。

**阶条件**（必要条件）: $L \geq K$。

**秩条件通过线性投影理解**:
- $\mathbf{x}^* = \mathbf{z}\boldsymbol{\Pi}$，$\boldsymbol{\Pi} = [E(\mathbf{Z}'\mathbf{Z})]^{-1}E(\mathbf{Z}'\mathbf{X})$
- $\mathbf{X} = \mathbf{X}^* + \mathbf{R}$, $E(\mathbf{Z}'\mathbf{R}) = \mathbf{0}$
- $E(\mathbf{X}^{*'}\mathbf{X}) = E(\mathbf{X}'\mathbf{Z})[E(\mathbf{Z}'\mathbf{Z})]^{-1}E(\mathbf{Z}'\mathbf{X})$ 满秩 $\Leftrightarrow$ $E(\mathbf{Z}'\mathbf{X})$ 秩为 $K$

**定理 1 (2SLS的一致性)**: 在2SLS.1和2SLS.2下，$\hat{\boldsymbol{\beta}}_{2SLS} \xrightarrow{p} \boldsymbol{\beta}$。

### 3.2 渐近正态性

**假设 2SLS.3** (同方差): $E(u^2 \mathbf{z}'\mathbf{z}) = \sigma^2 E(\mathbf{z}'\mathbf{z})$，其中 $\sigma^2 = E(u^2)$。

充分条件: $E(u^2 | \mathbf{z}) = \sigma^2$ 或 $\text{Var}(u | \mathbf{z}) = \sigma^2$（在 $E(u|\mathbf{z}) = 0$ 下）。

**定理 2 (2SLS的渐近正态性)**: 在2SLS.1-2SLS.3下，

$$\sqrt{N}(\hat{\boldsymbol{\beta}}_{2SLS} - \boldsymbol{\beta}) \xrightarrow{d} N\left(\mathbf{0}, \; \sigma^2 \left\{E(\mathbf{X}'\mathbf{Z})[E(\mathbf{Z}'\mathbf{Z})]^{-1}E(\mathbf{Z}'\mathbf{X})\right\}^{-1}\right) \tag{5.24}$$

### 3.3 2SLS残差的重要注意

**2SLS残差的正确定义**:

$$\hat{u}_i = y_i - \mathbf{X}_i \hat{\boldsymbol{\beta}}_{2SLS} \tag{5.25}$$

其中 $\mathbf{X}_i$ 使用**原始**的 $x_{iK}$（而非 $\hat{x}_{iK}$）。

**常见错误**: 第二阶段回归的残差 $y_i - \hat{\mathbf{X}}_i\hat{\boldsymbol{\beta}}_{2SLS}$ 不是正确的2SLS残差——因为第二阶段OLS使用 $\hat{x}_{iK}$ 替代 $x_{iK}$。

方差估计: $\hat{\sigma}^2 = \frac{1}{N-K}\sum_{i=1}^{N} \hat{u}_i^2$

### 3.4 渐近有效性

**定理 3 (2SLS的相对有效性)**: 在2SLS.1-2SLS.3下，2SLS在所有使用 $\mathbf{z}$ 的线性组合作为工具变量的IV估计量中是**渐近最有效**的。

### 3.5 异方差稳健推断

当2SLS.3不成立时（异方差存在），用稳健方差矩阵：

$$\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}) = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1} \left( \sum_{i=1}^{N} \hat{u}_i^2 \hat{\mathbf{X}}_i'\hat{\mathbf{X}}_i \right) (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1} \tag{5.34}$$

有时乘 $\frac{N}{N-K}$ 作为自由度校正。对角元素的平方根是异方差稳健2SLS标准误。

---

## 4. 假设检验

### 4.1 2SLS的F型检验

检验 $H_0: \boldsymbol{\beta}_2 = \mathbf{0}$（$K_2$ 个约束）：

$$F = \frac{(\hat{S}\hat{S}R_r - \hat{S}\hat{S}R_{ur}) / K_2}{\hat{S}\hat{S}R_{ur} / (N-K)} \sim F_{K_2, N-K}$$

**注意**: $\hat{S}\hat{S}R$ 使用**第二阶段的残差平方和**（$y$ 对 $\hat{\mathbf{x}}$ 回归的残差），而非2SLS残差。

### 4.2 Example 5.3 (父母与配偶教育作为工具变量)

$$\log(\text{wage}) = \beta_0 + \beta_1 \text{exper} + \beta_2 \text{exper}^2 + \beta_3 \text{educ} + u$$

- IV: motheduc, fatheduc, huseduc
- 第一阶段: $\text{educ} = \delta_0 + \delta_1 \text{exper} + \delta_2 \text{exper}^2 + \theta_1 \text{motheduc} + \theta_2 \text{fatheduc} + \theta_3 \text{huseduc} + r$
- F检验 $H_0: \theta_1 = \theta_2 = \theta_3 = 0$（检验工具变量的相关性）
- 2SLS教育回报约8%，OLS估计约10.7%

---

## 5. 弱工具变量问题

### 5.1 IV不一致性的公式

$$\text{plim } \hat{\beta}_1 = \beta_1 + \frac{\text{Cov}(z_1, u)}{\text{Cov}(z_1, x_1)} \tag{5.36}$$

等价地（相关性形式）:

$$\text{plim } \hat{\beta}_1 = \beta_1 + \frac{\sigma_u}{\sigma_{x_1}} \cdot \frac{\text{Corr}(z_1, u)}{\text{Corr}(z_1, x_1)} \tag{5.37}$$

### 5.2 关键含义

即使 $\text{Corr}(z_1, u)$ 很小（工具变量几乎外生），若 $\text{Corr}(z_1, x_1)$ **也很小**（弱工具变量），IV的不一致性可以**任意大**。

$$\text{plim } \hat{\beta}_{IV} - \beta_1 = \frac{\sigma_u}{\sigma_{x_1}} \cdot \frac{\text{Corr}(z_1, u)}{\text{Corr}(z_1, x_1)} \to \infty \quad \text{当 } \text{Corr}(z_1, x_1) \to 0$$

### 5.3 诊断与应对

- **经验规则 (Stock-Yogo)**: 第一阶段F统计量 > 10
- 若 $F < 10$ → 弱IV，2SLS推断不可靠
- 过度识别时检验第一阶段F统计量: $H_0: \theta_1 = ... = \theta_M = 0$

---

## 6. 遗漏变量视角下的IV

### 6.1 模型

$$y = \beta_0 + \beta_1 x_1 + \cdots + \beta_K x_K + \gamma q + \epsilon \tag{5.45}$$

其中 $q$ 为遗漏变量，$E(\epsilon | \mathbf{X}, q) = 0$。$\mathbf{x}$ 的某些元素与 $q$ 相关 → 内生性。

### 6.2 IV作为解决方案

将 $q$ 放入误差项 $u = \gamma q + \epsilon$。工具变量需满足：

1. 在 $E(y | \mathbf{X}, q)$ 中**冗余**
2. 与遗漏变量 $q$ **不相关**
3. 与 $\mathbf{X}$ 的内生元素**充分相关**

2SLS施加于 $y = \mathbf{X}\boldsymbol{\beta} + u$（$u = \gamma q + \epsilon$）→ 一致且渐近正态。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $\text{Cov}(z_1, u) = 0$ | IV外生性条件 (5.3) |
| $\theta_1 \neq 0$ | IV相关性（约简型） (5.5) |
| $\text{rank } E(\mathbf{Z}'\mathbf{X}) = K$ | 秩条件 (5.10) |
| $\boldsymbol{\beta} = [E(\mathbf{Z}'\mathbf{X})]^{-1}E(\mathbf{Z}'\mathbf{y})$ | IV识别 (5.11) |
| $\hat{\boldsymbol{\beta}}_{2SLS} = (\mathbf{X}'\mathbf{P}_Z\mathbf{X})^{-1}\mathbf{X}'\mathbf{P}_Z\mathbf{y}$ | 2SLS估计量 |
| $\hat{u}_i = y_i - \mathbf{X}_i\hat{\boldsymbol{\beta}}_{2SLS}$ | 2SLS残差（非第二阶段残差） (5.25) |
| $E(u^2\mathbf{z}'\mathbf{z}) = \sigma^2E(\mathbf{z}'\mathbf{z})$ | 同方差假设 (2SLS.3) |
| $\text{plim } \hat{\beta}_1 = \beta_1 + \frac{\sigma_u}{\sigma_{x_1}} \cdot \frac{\text{Corr}(z_1,u)}{\text{Corr}(z_1,x_1)}$ | 弱IV偏误 (5.37) |
| $\widehat{\text{Avar}} = (\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}(\sum \hat{u}_i^2\hat{\mathbf{X}}_i'\hat{\mathbf{X}}_i)(\hat{\mathbf{X}}'\hat{\mathbf{X}})^{-1}$ | 2SLS稳健方差 (5.34) |
