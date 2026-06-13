---
title: 二元响应模型
type: framework
domain: econometrics
source: Chapter VII — Binary Response Models
created: 2026-05-19
updated: 2026-05-24
---

# Ch7: 二元响应模型 (Binary Response Models)

## 核心问题

因变量为二值变量 $y \in \{0, 1\}$。关注**响应概率** $p(\mathbf{x}) = P(y=1|\mathbf{x})$ 及解释变量的偏效应。

---

## 1. 线性概率模型 (LPM)

### 1.1 模型设定

$$P(y=1|\mathbf{x}) = \beta_0 + \beta_1 x_1 + \cdots + \beta_K x_K \tag{15.4}$$

$\beta_j$ 直接解释为 $x_j$ 一单位变化带来的响应概率变化。

### 1.2 基本性质

$$\begin{aligned}
E(y|\mathbf{x}) &= 1 \cdot p(\mathbf{x}) + 0 \cdot (1-p(\mathbf{x})) = p(\mathbf{x}) = \mathbf{x}\boldsymbol{\beta} \\
\text{Var}(y|\mathbf{x}) &= p(\mathbf{x})[1-p(\mathbf{x})] = \mathbf{x}\boldsymbol{\beta}(1-\mathbf{x}\boldsymbol{\beta})
\end{aligned}$$

→ OLS一致甚至无偏，但**方差随$\mathbf{x}$变化**（必然异方差）。

### 1.3 估计策略

- 使用**异方差稳健标准误**和t统计量
- **WLS**: $\hat{\sigma}_i = [\hat{y}_i(1-\hat{y}_i)]^{1/2}$，回归 $y_i/\hat{\sigma}_i$ 对 $\mathbf{x}_i/\hat{\sigma}_i$
- 若 $\hat{y}_i \notin [0,1]$，WLS无法直接应用

### 1.4 优点与局限

- **优点**: 偏效应容易解释（直接为概率变化），近似平均偏效应好
- **局限**: 预测概率可能超出 $[0,1]$；偏效应恒定为常数

---

## 2. Logit 与 Probit 模型

### 2.1 指数模型

$$P(y=1|\mathbf{x}) = G(\mathbf{x}\boldsymbol{\beta}) \equiv p(\mathbf{x}) \tag{15.8}$$

响应概率仅通过**指数** $\mathbf{x}\boldsymbol{\beta} = \beta_1 + \beta_2 x_2 + \cdots + \beta_K x_K$ 依赖 $\mathbf{x}$。

### 2.2 潜变量推导

$$y^* = \mathbf{x}\boldsymbol{\beta} + e, \quad y = \mathbf{1}[y^* > 0]$$

$e$ 连续分布，独立于 $\mathbf{x}$，关于零对称:

$$P(y=1|\mathbf{x}) = P(y^* > 0|\mathbf{x}) = P(e > -\mathbf{x}\boldsymbol{\beta}|\mathbf{x}) = 1 - G(-\mathbf{x}\boldsymbol{\beta}) = G(\mathbf{x}\boldsymbol{\beta})$$

### 2.3 Probit

$$G(z) = \Phi(z) = \int_{-\infty}^{z} \phi(v)dv, \quad \phi(z) = (2\pi)^{-1/2}\exp(-z^2/2) \tag{15.10}$$

### 2.4 Logit

$$G(z) = \Lambda(z) = \frac{e^z}{1+e^z} \tag{15.12}$$

性质: $0 < G(z) < 1$ 对所有 $z \in \mathbb{R}$；$G$ 严格递增。

### 2.5 偏效应

**连续变量 $x_j$**:

$$\frac{\partial P(y=1|\mathbf{x})}{\partial x_j} = g(\mathbf{x}\boldsymbol{\beta}) \cdot \beta_j$$

其中 $g(z) = dG(z)/dz$:
- Logit: $g(z) = \Lambda(z)[1-\Lambda(z)]$
- Probit: $g(z) = \phi(z)$

**偏效应比**: $\frac{\partial p/\partial x_j}{\partial p/\partial x_h} = \frac{\beta_j}{\beta_h}$（比例恒定）

**二元变量 $x_K$**: $G(\beta_1 + \beta_2 x_2 + \cdots + \beta_{K-1} x_{K-1} + \beta_K) - G(\beta_1 + \beta_2 x_2 + \cdots + \beta_{K-1} x_{K-1})$

**$\beta_j$ 的符号给出偏效应的方向**（$g>0$），但 $\beta_j$ **本身不是偏效应**。

---

## 3. 最大似然估计 (MLE)

### 3.1 似然函数

给定 $\mathbf{x}_i$ 下 $y_i$ 的密度:

$$f(y_i | \mathbf{x}_i; \boldsymbol{\beta}) = [G(\mathbf{x}_i\boldsymbol{\beta})]^{y_i} [1 - G(\mathbf{x}_i\boldsymbol{\beta})]^{1-y_i}, \quad y_i = 0,1 \tag{15.16}$$

**对数似然**: $\ell_i(\boldsymbol{\beta}) = y_i \log G(\mathbf{x}_i\boldsymbol{\beta}) + (1-y_i) \log[1 - G(\mathbf{x}_i\boldsymbol{\beta})]$

$$\mathcal{L}(\boldsymbol{\beta}) = \sum_{i=1}^{N} \ell_i(\boldsymbol{\beta})$$

### 3.2 得分函数与Hessian

$$\mathbf{s}_i(\boldsymbol{\beta}) = \frac{g(\mathbf{x}_i\boldsymbol{\beta}) \mathbf{x}_i' [y_i - G(\mathbf{x}_i\boldsymbol{\beta})]}{G(\mathbf{x}_i\boldsymbol{\beta})[1 - G(\mathbf{x}_i\boldsymbol{\beta})]}$$

**期望Hessian**: $-E[\mathbf{H}_i(\boldsymbol{\beta}) | \mathbf{x}_i] = \frac{g(\mathbf{x}_i\boldsymbol{\beta})^2 \mathbf{x}_i'\mathbf{x}_i}{G(\mathbf{x}_i\boldsymbol{\beta})[1 - G(\mathbf{x}_i\boldsymbol{\beta})]} \equiv \mathbf{A}(\mathbf{x}_i, \boldsymbol{\beta})$

### 3.3 渐近方差

$$\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}) = \left[ \sum_{i=1}^{N} \frac{g(\mathbf{x}_i\hat{\boldsymbol{\beta}})^2 \mathbf{x}_i'\mathbf{x}_i}{G(\mathbf{x}_i\hat{\boldsymbol{\beta}})[1 - G(\mathbf{x}_i\hat{\boldsymbol{\beta}})]} \right]^{-1} = \hat{\mathbf{V}}$$

**异方差稳健**: $\widehat{\text{Avar}}(\hat{\boldsymbol{\beta}}) = [\sum \mathbf{H}_i]^{-1} [\sum \mathbf{s}_i \mathbf{s}_i'] [\sum \mathbf{H}_i]^{-1}$

---

## 4. 假设检验

### 4.1 三种经典检验

$$H_0: \boldsymbol{\gamma} = \mathbf{0} \; \text{ in } \; P(y=1|\mathbf{x},\mathbf{z}) = G(\mathbf{x}\boldsymbol{\beta} + \mathbf{z}\boldsymbol{\gamma})$$

| 检验 | 公式 |
|------|------|
| Wald | 软件直接计算 |
| LR | $LR = 2(\mathcal{L}_{ur} - \mathcal{L}_r) \xrightarrow{d} \chi^2_Q$ |
| LM (Score) | 受约束模型残差回归 → $NR_u^2 \sim \chi^2_Q$ |

### 4.2 LM检验的具体步骤

1. 估计受约束模型 → $G_i, \hat{u}_i = y_i - G(\mathbf{x}_i\hat{\boldsymbol{\beta}})$
2. OLS回归: $\frac{\hat{u}_i}{\sqrt{G_i(1-G_i)}}$ 对 $\frac{g_i \mathbf{x}_i}{\sqrt{G_i(1-G_i)}}$ 和 $\frac{g_i \mathbf{z}_i}{\sqrt{G_i(1-G_i)}}$
3. $LM = N \cdot R_u^2 \xrightarrow{d} \chi^2_Q$

LM检验在 $\mathbf{z}$ 维数大时特别方便（无需估计无约束模型）。

### 4.3 非线性假设

$H_0: \mathbf{c}(\boldsymbol{\beta}) = \mathbf{0}$ ($Q \times 1$):

$$W = \mathbf{c}(\hat{\boldsymbol{\beta}})' [\nabla_{\boldsymbol{\beta}}\mathbf{c}(\hat{\boldsymbol{\beta}}) \cdot \hat{\mathbf{V}} \cdot \nabla_{\boldsymbol{\beta}}\mathbf{c}(\hat{\boldsymbol{\beta}})']^{-1} \mathbf{c}(\hat{\boldsymbol{\beta}})$$

---

## 5. 偏效应：APE vs PEA

### 5.1 均值处偏效应 (PEA)

$$\widehat{\text{PEA}}_j = g(\bar{\mathbf{x}}\hat{\boldsymbol{\beta}}) \cdot \hat{\beta}_j$$

问题: $\bar{\mathbf{x}}$ 对离散变量无意义。

### 5.2 平均偏效应 (APE)

**连续变量**:

$$\widehat{\text{APE}}_j = \hat{\beta}_j \cdot \frac{1}{N} \sum_{i=1}^{N} g(\mathbf{x}_i\hat{\boldsymbol{\beta}})$$

**离散变量**:

$$\widehat{\text{APE}}_K = \frac{1}{N} \sum_{i=1}^{N} \left[ G(\hat{\beta}_1 + \hat{\beta}_2 x_{i2} + \cdots + \hat{\beta}_K) - G(\hat{\beta}_1 + \hat{\beta}_2 x_{i2} + \cdots + 0) \right]$$

### 5.3 Probit与Logit系数的缩放关系

| 模型 | $g(0)$ | 关系 |
|------|--------|------|
| LPM | 1 | 基准 |
| Probit | $\phi(0) \approx 0.4$ | $\hat{\beta}_{Logit} \approx 1.6 \hat{\beta}_{Probit}$ |
| Logit | $\Lambda(0)[1-\Lambda(0)] = 0.25$ | $\hat{\beta}_{Logit} \approx 4 \hat{\beta}_{LPM}$ |

→ 不同模型的系数量级不同，但**APE可比**。

### 5.4 Stoker (1986) 的结果

若 $\mathbf{x}$ 服从多元正态（凸无界支撑），$y$ 对 $(1,\mathbf{x})$ 线性投影的斜率系数可以**一致估计APE**——为LPM的广泛适用性提供了理论支撑。

---

## 6. 设定问题

### 6.1 遗漏异质性 (Neglected Heterogeneity)

$$P(y=1|\mathbf{x}, c) = \Phi(\mathbf{x}\boldsymbol{\beta} + \gamma c), \quad c | \mathbf{x} \sim N(0, \tau^2)$$

→ $P(y=1|\mathbf{x}) = \Phi(\mathbf{x}\boldsymbol{\beta} / \sigma)$，$\sigma^2 = \gamma^2\tau^2 + 1$

**衰减偏误**: $|\beta_j / \sigma| < |\beta_j|$——系数被缩放。

**但APE仍可一致估计**: $E[\beta_j \phi(\mathbf{x}\boldsymbol{\beta} + \gamma c)] = (\beta_j/\sigma) \phi(\mathbf{x}\boldsymbol{\beta}/\sigma)$

若 $c$ 与 $\mathbf{x}$ 相关 ($c|\mathbf{x} \sim N(\mathbf{x}\boldsymbol{\delta}, \eta^2)$) → Probit有偏且不一致（偏误+异质性）。

### 6.2 连续内生解释变量 (Rivers & Vuong 1988 控制函数)

$$y_1 = \mathbf{1}[\mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + u_1 > 0] \tag{15.39}$$
$$y_2 = \mathbf{z}\boldsymbol{\delta}_2 + v_2 \tag{15.40}$$

$(u_1, v_2)$ 独立于 $\mathbf{z}$，二元正态，$\text{Corr}(u_1, v_2) = \rho_1$。

$u_1 = \theta_1 v_2 + e_1$, $\theta_1 = \text{Cov}(v_2, u_1)/\text{Var}(v_2)$

$$y_1^* = \mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + \theta_1 v_2 + e_1$$

**Procedure 15.1**:
1. $y_2$ 对 $\mathbf{z}$ OLS → $\hat{v}_2$
2. $y_1$ 对 $\mathbf{z}_1, y_2, \hat{v}_2$ Probit → 得缩放系数

$\hat{v}_2$ 的t统计量检验 $H_0: \theta_1 = 0$（$y_2$ 外生）。

### 6.3 二元内生解释变量 (Bivariate Probit)

$$y_1 = \mathbf{1}[\mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1 y_2 + u_1 > 0] \tag{15.51}$$
$$y_2 = \mathbf{1}[\mathbf{z}\boldsymbol{\delta}_2 + v_2 > 0] \tag{15.52}$$

$(u_1, v_2)$ 独立于 $\mathbf{z}$，二元标准正态，$\rho_1 = \text{Corr}(u_1, v_2)$。

通过联合分布 $f(y_1, y_2 | \mathbf{z}) = f(y_1 | y_2, \mathbf{z}) \cdot f(y_2 | \mathbf{z})$ 进行MLE（**二元Probit模型**）。

**平均处理效应**: $\Phi(\mathbf{z}_1\boldsymbol{\delta}_1 + \alpha_1) - \Phi(\mathbf{z}_1\boldsymbol{\delta}_1)$

### 6.4 异方差Probit

若 $D(e|\mathbf{x}) = N(0, h(\mathbf{x}))$:

$$P(y=1|\mathbf{x}) = \Phi\left(\frac{\mathbf{x}\boldsymbol{\beta}}{\sqrt{h(\mathbf{x})}}\right)$$

→ 标准Probit不一致。

**APE估计**（对连续变量）:

$$\widehat{\text{APE}}_j = \frac{1}{N} \sum_{i=1}^{N} \phi\left(\frac{\mathbf{x}_i\hat{\boldsymbol{\beta}}}{\exp(\mathbf{x}_{i1}\hat{\boldsymbol{\delta}})}\right) \cdot \frac{\hat{\beta}_j}{\exp(\mathbf{x}_{i1}\hat{\boldsymbol{\delta}})}$$

---

## 7. 拟合优度

- **正确预测百分比**: $\tilde{y}_i = \mathbf{1}[G(\mathbf{x}_i\hat{\boldsymbol{\beta}}) \geq 0.5]$，$\frac{1}{N}\sum \mathbf{1}[y_i = \tilde{y}_i]$
- **McFadden's Pseudo-$R^2$**: $1 - \mathcal{L}_{ur} / \mathcal{L}_0$（$\mathcal{L}_0$ 为仅含截距的对数似然）

---

## 8. 面板二元响应模型

### 8.1 混合Probit/Logit

$$P(y_{it}=1 | \mathbf{x}_{it}) = G(\mathbf{x}_{it}\boldsymbol{\beta}) \tag{15.62}$$

用混合MLE估计。动态完备性: $P(y_{it}=1|\mathbf{x}_{it}, y_{i,t-1}, \mathbf{x}_{i,t-1},...) = P(y_{it}=1|\mathbf{x}_{it})$ (15.63)

检验: 加 $\hat{u}_{i,t-1}$ 到Probit中检验显著性。

### 8.2 随机效应Probit

$$P(y_{it}=1 | \mathbf{x}_{it}, c_i) = \Phi(\mathbf{x}_{it}\boldsymbol{\beta} + c_i) \tag{15.66}$$

$c_i | \mathbf{x}_i \sim N(0, \sigma_c^2)$ (15.69)

$\rho = \frac{\sigma_c^2}{\sigma_c^2 + 1}$: 复合误差 $c_i + e_{it}$ 跨时期的相关系数。

CMLE通过对 $c_i$ 积分求得。

### 8.3 Chamberlain相关随机效应Probit

$$c_i | \mathbf{x}_i \sim N(\psi + \bar{\mathbf{x}}_i\boldsymbol{\xi}, \sigma_a^2) \tag{15.73}$$

$$y_{it}^* = \psi + \mathbf{x}_{it}\boldsymbol{\beta} + \bar{\mathbf{x}}_i\boldsymbol{\xi} + a_i + e_{it}$$

$H_0: \boldsymbol{\xi} = \mathbf{0}$ 检验RE vs CRE。

**APE**: $\frac{1}{N}\sum \Phi(\hat{\psi}_a + \tilde{\mathbf{x}}\hat{\boldsymbol{\beta}}_a + \bar{\mathbf{x}}_i\hat{\boldsymbol{\xi}}_a)$ (15.76)

### 8.4 固定效应Logit (条件Logit)

**优势**: 无需 $D(c_i|\mathbf{x}_i)$ 假设即可得 $\sqrt{N}$-一致的 $\hat{\boldsymbol{\beta}}$。

**T=2**: 仅用 $n_i = y_{i1}+y_{i2} = 1$ 的观测。对数似然:

$$\ell_i(\boldsymbol{\beta}) = \mathbf{1}[n_i=1] \cdot \left(w_i \log\Lambda(\Delta\mathbf{x}_i\boldsymbol{\beta}) + (1-w_i)\log[1-\Lambda(\Delta\mathbf{x}_i\boldsymbol{\beta})]\right) \tag{15.81}$$

$w_i = 1$ 若 $(y_{i1}=0, y_{i2}=1)$，否则 $0$。

**局限**: 不能估计偏效应（需知道 $c_i$）；需条件独立假设 (15.67)。

### 8.5 动态非观测效应模型

$$P(y_{it}=1 | y_{i,t-1},...,y_{i0}, \mathbf{z}_i, c_i) = G(\mathbf{z}_{it}\boldsymbol{\delta} + \rho y_{i,t-1} + c_i)$$

**状态依赖**: $H_0: \rho = 0$（控制 $c_i$ 后）。

**初始条件问题**: $y_{i0}$ 不是外生给定的——与 $c_i$ 可能相关。

**Heckman (1981) 方法**: 假设 $a_i | (y_{i0}, \mathbf{z}_i) \sim N(\psi + \xi_0 y_{i0} + \mathbf{z}_i\boldsymbol{\xi}, \sigma_a^2)$，然后用标准RE Probit软件。

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $P(y=1\|\mathbf{x}) = G(\mathbf{x}\boldsymbol{\beta})$ | 指数模型 (15.8) |
| $\ell_i = y_i\log G_i + (1-y_i)\log(1-G_i)$ | 对数似然 |
| $\partial P/\partial x_j = g(\mathbf{x}\boldsymbol{\beta})\beta_j$ | 连续变量偏效应 |
| $\widehat{\text{APE}}_j = \hat{\beta}_j \cdot \frac{1}{N}\sum g(\mathbf{x}_i\hat{\boldsymbol{\beta}})$ | 连续变量APE |
| $\widehat{\text{Avar}} = [\sum \mathbf{A}(\mathbf{x}_i,\hat{\boldsymbol{\beta}})]^{-1}$ | MLE渐近方差 |
| $LR = 2(\mathcal{L}_{ur} - \mathcal{L}_r)$ | 似然比检验 |
| $y_1^* = \mathbf{z}_1\delta_1 + \alpha_1 y_2 + \theta_1 \hat{v}_2 + e_1$ | Rivers-Vuong CF |
| $\text{ATE} = \Phi(\mathbf{z}_1\delta_1 + \alpha_1) - \Phi(\mathbf{z}_1\delta_1)$ | 二元内生变量ATE |
| $\ell_i = \mathbf{1}[n_i=1] \cdot (w_i\log\Lambda(\Delta\mathbf{x}_i\beta) + ...)$ | 条件Logit (15.81) |
