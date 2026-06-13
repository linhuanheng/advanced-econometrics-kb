---
title: OLS回归与异方差
type: framework
domain: econometrics
source: Chapter I — Heteroskedasticity
created: 2026-05-19
updated: 2026-05-20
---

# Ch1: 异方差 (Heteroskedasticity)

## 核心问题

同方差假设 $\text{Var}(u|x_1,x_2,...,x_k) = \sigma^2$ 不成立时，OLS估计量的性质及补救措施。

---

## 1. 异方差的后果与偏误分析

### 1.1 对点估计的影响：无偏性与一致性不受影响

OLS估计量在MLR.1-MLR.4下（不要求MLR.5同方差）仍然是**无偏且一致**的。原因在于：

$$\hat{\beta}_1 = \beta_1 + \frac{\sum_{i=1}^n (x_i - \bar{x}) u_i}{\sum_{i=1}^n (x_i - \bar{x})^2}$$

- **有限样本**: $E(\hat{\beta}_j) = \beta_j$，只需 $E(u|x)=0$（MLR.4），不依赖于 $\text{Var}(u|x)$ 的形式
- **大样本**: $\text{plim}(\hat{\beta}_j) = \beta_j$

**关键区分**: 异方差不会像遗漏变量那样导致OLS系数估计产生偏误或不一致。这是异方差与内生性问题之间的本质区别。

### 1.2 对 $R^2$ 和 $\bar{R}^2$ 的影响：不受影响

$R^2$ 和 $\bar{R}^2$ 是总体 $R^2 = 1 - \sigma_u^2 / \sigma_y^2$ 的估计量。$\sigma_u^2$ 和 $\sigma_y^2$ 都是**无条件方差**，不受 $\text{Var}(u|x)$ 是否随 $x$ 变化的影响。且 $\text{SSR}/n \xrightarrow{p} \sigma_u^2$，$\text{SST}/n \xrightarrow{p} \sigma_y^2$，无论同方差还是异方差均成立。

### 1.3 对标准误的影响：核心偏误所在

**常规OLS方差公式**: $\widehat{\text{Var}}_{\text{OLS}}(\hat{\beta}_j) = \frac{\hat{\sigma}^2}{\text{SST}_j(1-R_j^2)}$，其中 $\hat{\sigma}^2 = \frac{\sum \hat{u}_i^2}{n-k-1}$

**问题**: 该公式假设每个观测的误差方差均等于 $\sigma^2$。在异方差下：
- $\hat{\sigma}^2$ 估计的是一个"平均"方差，掩盖了方差在不同 $x$ 值处的真实差异
- 标准误的**偏误方向不固定**——大多数情况下（横截面数据中）OLS标准误**低估**真实标准误（导致t统计量偏高，过度拒绝原假设），但也可能高估
- **大样本不能解决此问题**: 即使 $n \to \infty$，常规OLS标准误仍不一致

**偏误的直觉理解**: 异方差意味着某些观测的信息量少于其他观测。OLS给所有观测相同权重，导致对估计精度的错误判断——通常表现为过于乐观。

### 1.4 对推断的影响

| 统计量 | 偏误后果 |
|--------|---------|
| t统计量 | 分布不再是 $t_{n-k-1}$，即使在 $n\to\infty$ 下也不标准正态（若使用常规标准误） |
| F统计量 | 不再服从F分布 |
| LM统计量 | 不再有渐近 $\chi^2$ 分布 |
| 置信区间 | 覆盖率不等于名义水平 |
| 检验 size | 实际显著性水平 $\neq$ 名义显著性水平（通常偏高→过度拒绝） |

### 1.5 对效率的影响

OLS不再是BLUE（最佳线性无偏估计量）。存在其他线性估计量（如GLS/WLS）具有更小的方差：
$$\text{Var}(\hat{\beta}_{\text{GLS}}) \leq \text{Var}(\hat{\beta}_{\text{OLS}})$$

OLS也不是渐近有效的——存在其他一致估计量（如FGLS）具有更小的渐近方差。

---

## 2. 解决方案一：异方差稳健推断

### 2.1 核心思路

不改变OLS点估计（反正它是无偏一致的），而是修正方差估计量使其在任意异方差形式下都是一致的。

### 2.2 简单回归中的White稳健标准误

考虑模型 $y_i = \beta_0 + \beta_1 x_i + u_i$，$\text{Var}(u_i|x_i) = \sigma_i^2$。

#### 推导起点：OLS估计量的表达

$$\hat{\beta}_1 = \beta_1 + \frac{\sum_{i=1}^n (x_i - \bar{x}) u_i}{\sum_{i=1}^n (x_i - \bar{x})^2}$$

| 符号 | 含义 |
|------|------|
| $\hat{\beta}_1$ | $\beta_1$ 的OLS估计量 |
| $\beta_1$ | 真实斜率参数（固定未知常数） |
| $x_i$ | 第 $i$ 个观测的自变量值 |
| $\bar{x}$ | 样本均值 $\frac{1}{n}\sum_{i=1}^n x_i$ |
| $x_i - \bar{x}$ | $x_i$ 的离差——决定第 $i$ 个观测对斜率估计的"杠杆" |
| $u_i$ | 不可观测的真实误差 |
| $n$ | 样本容量 |

#### 异方差下 $\hat{\beta}_1$ 的真实方差

$$\text{Var}(\hat{\beta}_1) = \frac{\sum_{i=1}^n (x_i - \bar{x})^2 \sigma_i^2}{\text{SST}_x^2} \tag{8.2}$$

| 符号 | 含义 |
|------|------|
| $\sigma_i^2$ | $\text{Var}(u_i\mid x_i)$，异方差下随 $i$ 变化（若同方差则 $\sigma_i^2 = \sigma^2$） |
| $\text{SST}_x$ | $\sum_{i=1}^n (x_i-\bar{x})^2$，$x$ 的总平方和 |
| 分子 | 以 $(x_i-\bar{x})^2$ 为权重的误差方差加权和 |

#### White方差估计量

White (1980) 用OLS残差平方 $\hat{u}_i^2$ 替代不可观测的 $\sigma_i^2$：

$$\widehat{\text{Var}}(\hat{\beta}_1) = \frac{\sum_{i=1}^n (x_i - \bar{x})^2 \hat{u}_i^2}{\text{SST}_x^2}$$

| 符号 | 含义 |
|------|------|
| $\hat{u}_i$ | OLS残差，$\hat{u}_i = y_i - \hat{\beta}_0 - \hat{\beta}_1 x_i$ |
| $\hat{u}_i^2$ | 平方残差——作为 $\sigma_i^2$ 的代理变量 |

**一致性保证**: 由大数定律和中心极限定理，$n \cdot \frac{\sum (x_i-\bar{x})^2 \hat{u}_i^2}{\text{SST}_x^2} \xrightarrow{p} \frac{E[(x_i-\mu_x)^2 u_i^2]}{(\sigma_x^2)^2}$

**稳健标准误**: $\text{se}_{\text{robust}}(\hat{\beta}_1) = \frac{\sqrt{\sum (x_i-\bar{x})^2 \hat{u}_i^2}}{\text{SST}_x}$

**稳健t统计量**: $t_{\text{robust}} = \frac{\hat{\beta}_1 - \beta_1^0}{\text{se}_{\text{robust}}(\hat{\beta}_1)} \xrightarrow{d} N(0,1)$

### 2.3 多元回归中的White稳健标准误

对多元模型 $y_i = \beta_0 + \beta_1 x_{i1} + \cdots + \beta_k x_{ik} + u_i$：

$$\widehat{\text{Var}}(\hat{\beta}_j) = \frac{\sum_{i=1}^n \hat{r}_{ij}^2 \hat{u}_i^2}{\text{SSR}_j^2} \tag{8.4}$$

| 符号 | 含义 | 深层理解 |
|------|------|---------|
| $\hat{\beta}_j$ | 第 $j$ 个系数的OLS估计量 | $j=1,2,...,k$ |
| $\hat{u}_i$ | OLS残差 | $y_i$ 中不能被所有 $x$ 线性解释的部分 |
| $\hat{r}_{ij}$ | 辅助回归残差 | $x_{ij}$ 对其他所有自变量回归的残差——衡量 $x_j$ 的"净变异" |
| $\text{SSR}_j$ | 辅助回归残差平方和 | $\sum \hat{r}_{ij}^2 = \text{SST}_j(1-R_j^2)$，$R_j^2$ 高则 $\text{SSR}_j$ 小→方差膨胀 |
| 分子 | $\sum \hat{r}_{ij}^2 \hat{u}_i^2$ | 以 $x_j$ 的独立变异为权重，对平方残差加权 |

#### 自由度校正

$$\widehat{\text{Var}}_{\text{adj}}(\hat{\beta}_j) = \frac{n}{n-k-1} \cdot \frac{\sum \hat{r}_{ij}^2 \hat{u}_i^2}{\text{SSR}_j^2}$$

校正因子 $\frac{n}{n-k-1} \to 1$（当 $n\to\infty$）；若所有 $\hat{u}_i^2$ 相等，校正后恰好还原为OLS标准误。

### 2.4 稳健 vs 常规OLS标准误

| 对比维度 | 常规OLS标准误 | White稳健标准误 |
|----------|-------------|----------------|
| 需同方差假设 | 是 | 否 |
| 有限样本分布 | 正态误差下精确 $t$ | 仅大样本渐近有效 |
| 大样本有效性 | 仅同方差下 | 始终有效 |
| 使用建议 | 有同方差证据时使用 | 横截面应用默认推荐 |

### 2.5 稳健F统计量与稳健LM统计量

#### 2.5.1 为什么需要稳健版本的F和LM？

常规F统计量和LM统计量的有效性都建立在**同方差假设**（MLR.5）之上。在异方差条件下：

| 常规统计量 | 同方差下 | 异方差下 |
|-----------|---------|---------|
| F统计量 | $F_{q, n-k-1}$（精确） | 分布未知，实际显著性水平 $\neq$ 名义水平 |
| LM统计量 ($nR^2$) | $\chi^2_q$（渐近） | 分布未知，不再渐近卡方 |
| 两者都依赖于 | $\text{Var}(u\|x) = \sigma^2$ | 该假设不成立 → 推断失效 |

**异方差稳健Wald统计量**（亦称异方差稳健F统计量）和**异方差稳健LM统计量**解决了这一问题——它们利用White型稳健方差估计，在任何异方差形式下都渐近有效。

#### 2.5.2 异方差稳健Wald/F统计量

**核心思想**: 用异方差稳健方差-协方差矩阵 $\hat{\mathbf{V}}$ 替代同方差下的 $\hat{\sigma}^2(\mathbf{X}'\mathbf{X})^{-1}$。

**一般框架**: 考虑 $q$ 个线性约束 $H_0: \mathbf{R}\boldsymbol{\beta} = \mathbf{r}$，其中 $\mathbf{R}$ 是 $q \times (k+1)$ 约束矩阵。

**稳健Wald统计量**:

$$W = (\mathbf{R}\hat{\boldsymbol{\beta}} - \mathbf{r})' \left[ \mathbf{R} \hat{\mathbf{V}} \mathbf{R}' \right]^{-1} (\mathbf{R}\hat{\boldsymbol{\beta}} - \mathbf{r}) \xrightarrow{d} \chi^2_q$$

| 符号 | 含义 | 说明 |
|------|------|------|
| $\hat{\boldsymbol{\beta}}$ | OLS估计量向量 | $(k+1) \times 1$，含截距 |
| $\mathbf{R}$ | $q \times (k+1)$ 约束矩阵 | 每行代表一个线性约束 |
| $\mathbf{r}$ | $q \times 1$ 约束值向量 | 通常 $\mathbf{r} = \mathbf{0}$ |
| $\hat{\mathbf{V}}$ | 稳健方差-协方差矩阵 | $\hat{\mathbf{V}} = (\mathbf{X}'\mathbf{X})^{-1} \left( \sum_{i=1}^n \hat{u}_i^2 \mathbf{x}_i' \mathbf{x}_i \right) (\mathbf{X}'\mathbf{X})^{-1}$（注：这是White稳健方差矩阵的另一种等价写法） |
| $q$ | 约束条件个数 | 被检验的排除性约束数量 |

**以排除性约束为例**（最常见的应用）：

检验 $H_0: \beta_4 = 0, \beta_5 = 0$（排除 $x_4, x_5$）：

$$\mathbf{R} = \begin{pmatrix} 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 \end{pmatrix}, \quad \mathbf{r} = \begin{pmatrix} 0 \\ 0 \end{pmatrix}$$

$$\mathbf{R}\hat{\boldsymbol{\beta}} - \mathbf{r} = \begin{pmatrix} \hat{\beta}_4 \\ \hat{\beta}_5 \end{pmatrix}$$

则 $W = (\hat{\beta}_4, \hat{\beta}_5) \left[ \mathbf{R} \hat{\mathbf{V}} \mathbf{R}' \right]^{-1} \begin{pmatrix} \hat{\beta}_4 \\ \hat{\beta}_5 \end{pmatrix}$

**异方差稳健F统计量**: 对Wald统计量做简单变换即可得到更接近有限样本F分布的形式：

$$F_{\text{robust}} = \frac{W}{q}$$

在 $H_0$ 下 $F_{\text{robust}}$ 渐近服从 $F_{q, n-k-1}$ 分布（或直接报告 $W \sim \chi^2_q$）。

**计算上注意**: 许多计量软件（如Stata）直接提供异方差稳健F统计量。`reg y x1 x2 x3 x4 x5, robust` 后 `test x4 x5` 自动使用稳健协方差矩阵。

#### 2.5.3 异方差稳健LM统计量 — 回归实现方法

稳健Wald/F统计量要求软件能计算稳健协方差矩阵 $\hat{\mathbf{V}}$。不是所有回归软件包都自动支持这一功能。稳健LM统计量的一大优势是**只需OLS回归即可计算**，几乎任何软件包都能实现。

**模型设定**: 以原模型含5个自变量为例：

$$y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_3 x_3 + \beta_4 x_4 + \beta_5 x_5 + u$$

检验 $H_0: \beta_4 = 0, \beta_5 = 0$（$q = 2$）。

---

#### 2.5.4 稳健LM统计量的完整计算步骤

**Step 1: 估计受约束模型，获得受约束残差 $\tilde{u}_i$**

在 $H_0$ 约束下（排除 $x_4, x_5$），估计：

$$y_i = \tilde{\beta}_0 + \tilde{\beta}_1 x_{i1} + \tilde{\beta}_2 x_{i2} + \tilde{\beta}_3 x_{i3} + \tilde{u}_i$$

| 符号 | 含义 |
|------|------|
| $\tilde{u}_i$ | 受约束残差——包含了 $x_4$ 和 $x_5$ 可能具有的解释力（如果它们真的影响 $y$） |
| $\tilde{\beta}_j$ | 受约束OLS估计量 |

**Step 2: 将每个被排除的变量对所有包含变量回归，获得残差 $\tilde{r}_{ij}$**

分别对 $x_4$ 和 $x_5$ 进行辅助回归：

$$\begin{aligned}
x_{i4} &= \tilde{\gamma}_{40} + \tilde{\gamma}_{41} x_{i1} + \tilde{\gamma}_{42} x_{i2} + \tilde{\gamma}_{43} x_{i3} + \tilde{r}_{i4} \\
x_{i5} &= \tilde{\gamma}_{50} + \tilde{\gamma}_{51} x_{i1} + \tilde{\gamma}_{52} x_{i2} + \tilde{\gamma}_{53} x_{i3} + \tilde{r}_{i5}
\end{aligned}$$

| 符号 | 含义 | 深层理解 |
|------|------|---------|
| $\tilde{r}_{i4}$ | $x_4$ 对包含变量回归的残差 | $x_4$ 中与 $x_1, x_2, x_3$ 正交的"净变异"部分 |
| $\tilde{r}_{i5}$ | $x_5$ 对包含变量回归的残差 | 同上，$x_5$ 的"净变异" |
| 两个残差都来自 | 以**包含变量**为回归元的回归 | 与稳健方差公式(8.4)中 $\hat{r}_{ij}$ 的精神一致 |

**直觉**: $\tilde{r}_{i4}$ 和 $\tilde{r}_{i5}$ 分别提取了 $x_4$ 和 $x_5$ 中**不能被 $x_1, x_2, x_3$ 线性解释的部分**——这正是检验 $H_0: \beta_4 = \beta_5 = 0$ 所需要的信息。

**Step 3: 构造乘积项 $\tilde{r}_{ij} \cdot \tilde{u}_i$**

对每个观测 $i = 1, 2, ..., n$，计算：

$$s_{i1} = \tilde{r}_{i4} \cdot \tilde{u}_i, \quad s_{i2} = \tilde{r}_{i5} \cdot \tilde{u}_i$$

| 符号 | 含义 |
|------|------|
| $s_{i1}$ | 第 $i$ 个观测上 $x_4$ 的净变异与受约束残差的乘积 |
| $s_{i2}$ | 第 $i$ 个观测上 $x_5$ 的净变异与受约束残差的乘积 |

**直觉**: 若 $x_4$ 和 $x_5$ 确实影响 $y$（即 $\beta_4 \neq 0$ 或 $\beta_5 \neq 0$），则它们的净变异 $\tilde{r}_{i4}$ 和 $\tilde{r}_{i5}$ 应该与受约束残差 $\tilde{u}_i$ 系统相关——因为 $\tilde{u}_i$ 中包含了 $x_4$ 和 $x_5$ 的对 $y$ 的贡献。乘积项 $s_{ij}$ 的均值（或总和）反映了这种相关的强度。

**Step 4: 无截距回归——核心检验方程**

将**常数 1** 对 $s_{i1}$ 和 $s_{i2}$ 回归（**不含截距**）：

$$1 = \theta_1 s_{i1} + \theta_2 s_{i2} + e_i$$

| 符号 | 含义 | 为什么可以这样回归？ |
|------|------|---------------------|
| $1$ | 一个所有观测都等于1的"因变量" | 这是一个**构造性回归**——目的不是估计参数，而是检验 $s_{ij}$ 的均值是否为零 |
| $\theta_1, \theta_2$ | 乘积项的回归系数 | 若 $E(s_{ij}) = 0$（$H_0$ 成立），则 $\theta_j$ 应接近0 |
| $e_i$ | 此构造回归的残差 | 度量了"1"中不能被 $s_{i1}, s_{i2}$ 解释的部分 |
| 无截距 | "1"本身充当了类似截距的角色 | 加上截距会产生完全共线性 |

**为什么回归"1"而不是别的？** 这个构造源于得分检验（score test）的原理：在 $H_0$ 下，得分函数的样本均值应接近零。回归"1"对 $s_i$ 无截距等价于检验 $E(s_i) = 0$。

**Step 5: 计算稳健LM统计量**

$$\text{robust LM} = n - \text{SSR}^*$$

| 符号 | 含义 | 说明 |
|------|------|------|
| $\text{SSR}^*$ | Step 4中回归的残差平方和 | $\text{SSR}^* = \sum_{i=1}^n (1 - \hat{\theta}_1 s_{i1} - \hat{\theta}_2 s_{i2})^2$ |
| $n - \text{SSR}^*$ | 样本量减去残差平方和 | 等价于该回归的"解释平方和" (explained sum of squares) |

**等价表达**: $\text{robust LM} = n - \text{SSR}^*$ 也可写为 $n \times$（"1对 $s_i$ 的回归中1的被解释比例"），但前者数值上更直接。

在 $H_0$ 下，$\text{robust LM} \xrightarrow{d} \chi^2_q$（此处 $q=2$）。

**决策规则**: 若稳健LM统计量超过 $\chi^2_q$ 的临界值（或等价的p值小于显著性水平），则拒绝 $H_0$，认为至少有一个被排除变量与 $y$ 显著相关——即排除性约束不成立。

---

#### 2.5.5 一般模型下的稳健LM统计量

将上述步骤推广到任意模型 $y = \beta_0 + \beta_1 x_1 + \cdots + \beta_{K-1} x_{K-1} + \beta_K x_K + u$，检验 $H_0: \beta_{k-q+1} = \cdots = \beta_K = 0$：

| Step | 操作 | 产出 |
|------|------|------|
| 1 | $y$ 对 $1, x_1, ..., x_{K-q}$ 回归（受约束模型） | 受约束残差 $\tilde{u}_i$ |
| 2 | 对 $j = K-q+1, ..., K$：$x_j$ 对 $1, x_1, ..., x_{K-q}$ 回归 | $q$ 组辅助回归残差 $\tilde{r}_{ij}$ |
| 3 | 对每个 $j$ 计算 $s_{ij} = \tilde{r}_{ij} \cdot \tilde{u}_i$ | $q$ 个乘积序列 |
| 4 | 回归 $1$ 对 $s_{i, K-q+1}, ..., s_{iK}$（无截距） | 残差平方和 $\text{SSR}^*$ |
| 5 | $\text{robust LM} = n - \text{SSR}^*$ | $\xrightarrow{d} \chi^2_q$ |

---

#### 2.5.6 稳健LM统计量的推导逻辑（直觉）

| 层次 | 逻辑 |
|------|------|
| 1. 为什么用受约束残差 $\tilde{u}_i$？ | 如果 $H_0$ 成立，$\tilde{u}_i$ 应该仅包含噪声——与任何 $x$ 的函数都无关 |
| 2. 为什么用 $\tilde{r}_{ij}$（$x_j$ 的净变异）？ | 提取 $x_j$ 中"真正独立"的信息——消除与包含变量的共线性干扰 |
| 3. 为什么要乘积 $\tilde{r}_{ij} \cdot \tilde{u}_i$？ | 乘积项反映 $x_j$ 与残差之间的相关强度——本质上是得分函数的组成部分 |
| 4. 为什么要回归"1"？ | 在 $H_0$ 下，每个 $s_{ij}$ 的均值为0。回归"1"检验其均值是否显著偏离零——等价于检验矩条件 $E(\tilde{r}_j \tilde{u}) = 0$ |
| 5. 为什么 $n - \text{SSR}^*$？ | 当 $s_{ij}$ 的均值接近0时（$H_0$ 成立），它们几乎不能解释"1" → $\text{SSR}^* \approx n$ → 统计量接近0。均值偏离0时 → $\text{SSR}^* < n$ → 统计量增大 → 拒绝 $H_0$ |

---

#### 2.5.7 稳健F/Wald vs 稳健LM：方法对比

| 对比维度 | 稳健Wald/F | 稳健LM |
|---------|-----------|--------|
| 基础 | 无约束估计的稳健协方差矩阵 | 受约束模型残差 + 辅助回归 |
| 计算复杂度 | 需要 $\hat{\mathbf{V}}$（软件支持） | 仅需OLS回归（任何软件均可） |
| 适用场景 | 软件支持时首选 | 需要手动计算或软件不支持稳健Wald时 |
| 等价性 | 大样本下与稳健LM渐近等价 | 大样本下与稳健Wald渐近等价 |
| 小样本表现 | 如果 $n$ 较小，可能高估检验力 | 类似 |
| 软件实现 | Stata: `reg, robust` + `test` | 需要手动编程或专用命令 |

#### 2.5.8 一个完整的数值示例

考虑模型 $y = \beta_0 + \beta_1 \text{exper} + \beta_2 \text{exper}^2 + \beta_3 \text{educ} + u$，检验 $H_0: \beta_2 = 0$（$q = 1$，即是否需要 exper 的平方项？如果不需要平方项，是否说明模型设定有误？等等——这里只是示例）：

**Step 1**: $y$ 对 $1, \text{exper}, \text{educ}$ 回归 → $\tilde{u}_i$

**Step 2**: $\text{exper}^2$ 对 $1, \text{exper}, \text{educ}$ 回归 → $\tilde{r}_{i2}$（$x_2$ 中除去 $\text{exper}$ 和 $\text{educ}$ 线性影响后的剩余部分）

**Step 3**: $s_i = \tilde{r}_{i2} \cdot \tilde{u}_i$

**Step 4**: 回归 $1$ 对 $s_i$（无截距）→ 设 $\text{SSR}^* = 937.2$，$n = 1000$

**Step 5**: $\text{robust LM} = 1000 - 937.2 = 62.8$

查 $\chi^2_1$ 分布：临界值 $\chi^2_{1, 0.05} = 3.84$。$62.8 \gg 3.84$ → 强烈拒绝 $H_0$ → $\text{exper}^2$ 显著。

---

## 3. 解决方案二：加权最小二乘 (WLS/GLS)——寻求效率改进

### 3.1 基本假定

$$\text{Var}(u_i|x_i) = \sigma^2 h(x_i) = \sigma^2 h_i, \quad h_i > 0$$

$h(x)$ 是解释变量的某个函数，决定异方差结构。例如 $\text{Var}(\text{sav}_i | \text{inc}_i) = \sigma^2 \cdot \text{inc}_i$。

### 3.2 WLS变换

将原方程除以 $\sqrt{h_i}$：

$$\frac{y_i}{\sqrt{h_i}} = \beta_0 \frac{1}{\sqrt{h_i}} + \beta_1 \frac{x_{i1}}{\sqrt{h_i}} + \cdots + \beta_k \frac{x_{ik}}{\sqrt{h_i}} + \frac{u_i}{\sqrt{h_i}}$$

即 $y_i^* = \beta_0 x_{i0}^* + \beta_1 x_{i1}^* + \cdots + \beta_k x_{ik}^* + u_i^*$

**关键结果**: $\text{Var}(u_i^*|x_i) = \frac{1}{h_i}\text{Var}(u_i|x_i) = \frac{1}{h_i} \cdot \sigma^2 h_i = \sigma^2$

→ 变换后的方程满足所有经典线性回归假设（MLR.1-MLR.6）！

### 3.3 WLS估计量

$$\min_{b_0,...,b_k} \sum_{i=1}^n \frac{(y_i - b_0 - b_1 x_{i1} - \cdots - b_k x_{ik})^2}{h_i}$$

注意：平方残差除以 $h_i$（$\propto 1/\text{Var}$），变换变量除以 $\sqrt{h_i}$。

**直觉**: 给方差大的观测较小权重，方差小的观测较大权重。OLS给所有观测均等权重（$h_i \equiv 1$）。

### 3.4 GLS/WLS的性质

- **BLUE**: 若 $h(x)$ 正确设定，WLS是 $\beta_j$ 的BLUE
- **比OLS更有效**: $\text{Var}(\hat{\beta}_j^{\text{WLS}}) \leq \text{Var}(\hat{\beta}_j^{\text{OLS}})$
- **变换方程上的t/F检验精确有效**（若误差正态）

---

## 4. 解决方案三：可行GLS (FGLS)——$h(x)$ 未知时的做法

### 4.1 FGLS步骤

| 步骤 | 操作 | 目的 |
|------|------|------|
| 1 | OLS回归 $y$ 对 $x_1,...,x_k$ | 获得残差 $\hat{u}_i$ |
| 2 | 计算 $\log(\hat{u}_i^2)$ | 使得异方差函数的估计不受量纲影响 |
| 3 | $\log(\hat{u}_i^2)$ 对 $x_1,...,x_k$ 回归 | 估计 $\hat{g}_i = \hat{\alpha}_0 + \sum \hat{\alpha}_j x_{ij}$ |
| 4 | $\hat{h}_i = \exp(\hat{g}_i)$ | 指数函数保证 $\hat{h}_i > 0$ |
| 5 | 以 $1/\hat{h}_i$ 为权重进行WLS | FGLS估计 |

**为什么用指数函数？** $\text{Var}(u|x) = \sigma^2 \exp(\alpha_0 + \alpha_1 x_1 + \cdots)$ 自动保证方差为正；线性模型 $u^2 = \alpha_0 + \alpha_1 x_1 + e$ 可能产生负的预测值。

### 4.2 FGLS的有限样本性质（与GLS的对比）

| 性质 | GLS（$h(x)$已知） | FGLS（$\hat{h}(x)$估计） |
|------|------------------|------------------------|
| 无偏性 | 是 | **否**（$\hat{h}_i$ 使用了相同数据） |
| BLUE | 是 | **否** |
| 一致性 | 是 | **是** |
| 渐近有效 | 是 | **是**（异方差函数正确设定下优于OLS） |

### 4.3 异方差函数误设的后果

| 问题 | 影响 | 补救 |
|------|------|------|
| 一致性 | **不受影响**（$E(u\|x)=0$ 下WLS总能保持正交性） | 无需补救 |
| WLS标准误 | **无效**（即使大样本） | 对变换后的方程使用White稳健标准误 |
| 效率 | **不一定**优于OLS | 在强异方差下，即使误设的WLS通常仍比OLS好 |

### 4.4 关于OLS与WLS估计值差异的警告

若OLS和WLS产生**符号相反**或**量级差距巨大**的估计值，可能意味着：
1. $E(y|x)$ 的函数形式被误设（$\text{MLR.1}$ 不成立）
2. 零条件均值假设 $E(u|x)=0$（MLR.4）不成立

→ 可用Hausman检验正式比较OLS和WLS估计值。

### 4.5 FGLS的简化版本（White型FGLS）

仅使用 $\hat{y}$ 和 $\hat{y}^2$ 建模方差函数：

$$\log(\hat{u}_i^2) = \alpha_0 + \alpha_1 \hat{y}_i + \alpha_2 \hat{y}_i^2 + e_i$$

---

## 5. 解决方案选择：决策树

```
检测到异方差？
├── 否 → 报告常规OLS（更精确的有限样本性质）
└── 是 → 目标是什么？
    ├── 仅需有效的推断 → 报告White稳健标准误+t/F统计量
    │   （对系数估计无影响，仅修正标准误）
    └── 追求更高效率 → 是否知道异方差形式？
        ├── 是 → WLS/GLS（BLUE）
        └── 否 → FGLS
```

**实践建议**:
- 横截面数据中，默认报告异方差稳健标准误（不仅是"安全"的选择，且许多顶刊要求）
- 若样本量大（$n > 500$）且有强异方差证据，FGLS通常带来显著的效率改进
- 若有内生性怀疑（OLS与WLS估计值差异大），优先解决内生性问题

---

## 6. 异方差检验

### 6.1 Breusch-Pagan (BP) Test

$$H_0: E(u^2|x_1,...,x_k) = E(u^2) = \sigma^2$$

| 步骤 | 操作 |
|------|------|
| 1 | OLS回归 $y$ 对 $x_1,...,x_k$，得 $\hat{u}_i$ |
| 2 | 回归 $\hat{u}_i^2 = \delta_0 + \delta_1 x_{i1} + \cdots + \delta_k x_{ik} + e_i$，得 $R^2_{\hat{u}^2}$ |
| 3 | $LM = n \cdot R^2_{\hat{u}^2} \xrightarrow{d} \chi^2_k$ |

### 6.2 White Test

在BP基础上添加所有 $x_j^2$ 和 $x_j x_m$（$j \neq m$）。以3个变量为例，辅助回归含9个回归元（不含截距）。

- **优点**: 可检测任意异方差形式
- **缺点**: 自由度消耗大（$k$ 个变量 → $k(k+1)/2 + k$ 个回归元）

### 6.3 自由度节省版White Test

$$\hat{u}_i^2 = \alpha_0 + \alpha_1 \hat{y}_i + \alpha_2 \hat{y}_i^2 + e_i$$

$H_0: \alpha_1 = \alpha_2 = 0$。始终仅2个约束——因为 $\hat{y}_i^2$ 包含了所有 $x_j$ 的平方和交叉项信息。

**重要警告**: 函数形式误设可导致异方差检验显著——即使 $\text{Var}(y|x)$ 实为常数。建议**先做RESET**排除函数形式误设，再进行异方差检验。

---

## 7. 线性概率模型 (LPM) 的异方差

LPM中 $\text{Var}(y|x) = p(x)[1-p(x)]$，其中 $p(x) = P(y=1|x) = X\beta$。

**必然异方差**: $p(x)$ 随 $x$ 变化 → 方差也随 $x$ 变化。

**FGLS处理**: $\hat{h}_i = \hat{y}_i(1-\hat{y}_i)$。若 $\hat{y}_i \notin [0,1]$，需特殊调整或改用稳健OLS。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $\text{Var}(\hat{\beta}_1) = \frac{\sum (x_i-\bar{x})^2\sigma_i^2}{\text{SST}_x^2}$ | 异方差下OLS真实方差 (8.2) |
| $\widehat{\text{Var}}(\hat{\beta}_1) = \frac{\sum (x_i-\bar{x})^2 \hat{u}_i^2}{\text{SST}_x^2}$ | White稳健方差（简单回归） |
| $\widehat{\text{Var}}(\hat{\beta}_j) = \frac{\sum \hat{r}_{ij}^2 \hat{u}_i^2}{\text{SSR}_j^2}$ | White稳健方差（多元回归）(8.4) |
| $t_{\text{robust}} = \frac{\hat{\beta}_j - \beta_j^0}{\text{se}_{\text{robust}}(\hat{\beta}_j)}$ | 稳健t统计量 |
| $LM = n \cdot R^2_{\hat{u}^2} \sim \chi^2_k$ | BP检验统计量 |
| $\hat{h}_i = \exp(\hat{g}_i)$ | FGLS权重估计 |
| $\widehat{\text{Var}}(\hat{\beta}_j^{\text{WLS}}) = \widehat{\text{Var}}_{\text{robust}}[y^* \sim x^*]$ | WLS稳健标准误 |
