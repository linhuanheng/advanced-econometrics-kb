---
title: 处理效应评估
type: framework
domain: econometrics
source: Chapter VIII — Treatment Evaluation
created: 2026-05-19
updated: 2026-05-24
---

# Ch8: 处理效应评估 (Treatment Evaluation)

## 核心问题

评估二元处理（政策、项目）的因果效应。核心挑战：对每个个体只能观测到一个潜在结果——**反事实缺失问题**。

三种方法体系：(I) 可忽略性假设下的方法，(II) IV方法，(III) 断点回归设计。

---

## 1. 反事实框架 (Rubin Causal Model)

### 1.1 潜在结果

- $y_1$: 接受处理时的结果
- $y_0$: 未接受处理时的结果
- $w \in \{0, 1\}$: 处理指示变量

**观测结果**:

$$y = (1-w)y_0 + w y_1 = y_0 + w(y_1 - y_0) \tag{4}$$

### 1.2 处理效应定义

| 效应 | 定义 | 含义 |
|------|------|------|
| **ATE** | $\tau_{ate} = E[y_1 - y_0]$ | 全人群平均处理效应 |
| **ATT** | $\tau_{att} = E[y_1 - y_0 \mid w=1]$ | 处理组平均处理效应 |
| **LATE** | $\tau_{late} = E[y_1 - y_0 \mid w_1=1, w_0=0]$ | 依从者局部平均处理效应 |

**ATT与ATE的关系**:

令 $y_0 = \mu_0 + v_0$, $y_1 = \mu_1 + v_1$ ($\mu_g = E[y_g]$):

$$\tau_{att} = \tau_{ate} + E[v_1 - v_0 \mid w=1] \tag{7}$$

第二项是**个体特定处理收益**——若选择处理基于预期收益，则 $\tau_{att} > \tau_{ate}$。

### 1.3 三种观测机制

**Case I — 随机化处理** ($w \perp\!\!\!\perp (y_0, y_1)$):

$$\tau_{ate} = \tau_{att} = E[y|w=1] - E[y|w=0] \tag{5}$$

均值差估计量无偏、一致、渐近正态。

**Case II — 仅 $w \perp\!\!\!\perp y_0$**:

$$E[y|w=1] - E[y|w=0] = \tau_{att} + \underbrace{E[y_0|w=1] - E[y_0|w=0]}_{\text{选择偏误}} \tag{6}$$

若 $w$ 与 $y_0$ 独立 → 选择偏误为零 → 识别 $\tau_{att}$。

**Case III — 一般情况**: 需要通过控制 $\mathbf{x}$ 或使用IV来消除选择偏误。

---

## 2. 可忽略性假设下的方法

### 2.1 假设

#### 2.1.1 为什么需要可忽略性假设——因果推断的根本难题

处理效应评估的根本挑战在于**反事实缺失**：对每个个体，我们只能观测到一个潜在结果——要么 $y_1$（如果 $w=1$），要么 $y_0$（如果 $w=0$）。另一个潜在结果永远无法观测。因此，个体层面的处理效应 $y_1 - y_0$ 是**不可识别的**。

在总体层面，ATE = $E[y_1] - E[y_0]$ 涉及两个潜在结果的**边际分布**。问题在于：$E[y_1]$ 不能简单地用处理组样本均值来估计——因为处理组可能不是从总体中随机抽取的，其 $y_1$ 的分布可能与总体的 $y_1$ 分布不同。

**随机化实验的基准**：如果处理 $w$ 是通过抛硬币随机分配的（$w \perp\!\!\!\perp (y_0, y_1)$），那么：

$$E[y_1] = E[y_1 | w=1] = E[y | w=1]$$
$$E[y_0] = E[y_0 | w=0] = E[y | w=0]$$

此时简单的均值差 $E[y|w=1] - E[y|w=0]$ 就是 ATE。这是因果推断的**黄金标准**。

但在观测数据中，处理分配不是随机的——个体是否接受处理通常与其潜在结果相关（自选择）。**可忽略性假设是将观测研究"拉近"随机化实验的关键桥梁**。

#### 2.1.2 假设 ATE.1（全条件独立性/强可忽略性）

$$(y_0, y_1) \perp\!\!\!\perp w \mid \mathbf{x}$$

**数学含义**：给定协变量向量 $\mathbf{x}$，处理分配 $w$ 与潜在结果对 $(y_0, y_1)$ 在统计上独立。

**经济学含义（"基于可观测变量的选择"）**：一旦我们控制了 $\mathbf{x}$ 中的变量，处理组和对照组之间就不再有系统性的差异——两组在潜在结果分布上变得可比。任何影响处理选择的因素，只要同时也影响潜在结果，都已被 $\mathbf{x}$ 所包含。

**为什么称为"可忽略性"？** 因为在给定 $\mathbf{x}$ 后，处理分配机制可以被"忽略"——我们不需要显式地建模个体如何选择进入处理组。条件于 $\mathbf{x}$，处理组和对照组就像是从同一个分布中随机抽取的两个子样本。

**假设 ATE.1'（均值可忽略性/弱可忽略性）**：

$$E[y_0 | \mathbf{x}, w] = E[y_0 | \mathbf{x}], \quad E[y_1 | \mathbf{x}, w] = E[y_1 | \mathbf{x}]$$

这是一个更弱的版本——仅要求条件均值独立，不要求整个分布独立。ATE.1' 足以识别 ATE 和 ATT，因为 ATE 仅涉及 $y_0$ 和 $y_1$ 的**一阶矩**。但如果需要估计分位数处理效应或整个反事实分布，则需要完整的 ATE.1。

#### 2.1.3 假设 ATE.2（重叠/共同支撑条件）

$$0 < P(w=1 | \mathbf{x}) < 1, \quad \forall \mathbf{x}$$

**含义**：对于协变量空间中的每一个 $\mathbf{x}$ 值，都同时存在处理组和对照组的个体。即每个"类型"的人都有被处理和不被处理的可能。

**为什么不可或缺？** 如果对某些 $\mathbf{x}$，$P(w=1 | \mathbf{x}) = 0$（该类型从未被处理），那么我们永远无法观测到该类型的 $y_1$ → $\mu_1(\mathbf{x})$ 不可识别。同理，如果 $P(w=1 | \mathbf{x}) = 1$（该类型总是被处理），$\mu_0(\mathbf{x})$ 也不可识别。

在倾向得分加权中，重叠失效直接导致分母 $p(\mathbf{x})$ 或 $1-p(\mathbf{x})$ 为零 → IPW 估计量无定义。在匹配中，重叠失效意味着无法为某些处理个体找到可比的对照个体。

#### 2.1.4 关键实践问题：哪些变量应纳入 $\mathbf{x}$？

**应该纳入的变量**：
- **处理前测定的变量**（pretreatment variables）：年龄、教育（培训前）、既往收入史、基准测试成绩等。这些变量在处理前已经确定，不会被处理本身改变。
- **同时影响处理选择和潜在结果的变量**（confounders）：例如，在培训项目中，既往收入低的人更可能参加培训，且既往收入也预测未来收入。
- **过去的结果变量**（lagged outcomes）：$y$ 在处理前的取值通常是极强的混杂因素。

**绝不应纳入的变量**：
- **受处理影响的变量**（post-treatment variables）：如培训后获得的教育、培训后的职业资格。这些变量是处理的结果，纳入它们会"控制掉"处理效应的一部分——导致**过度控制偏误**（over-controlling bias）。

**例子**：$w$ 为职业培训指标，$y$ 为两年后的收入。不应在 $\mathbf{x}$ 中包含"培训后六个月内获得的证书"——因为获得证书本身就是培训的效果之一，控制它会低估培训对收入的总效应。

#### 2.1.5 可忽略性假设的可检验性

**核心结论：可忽略性假设是不可检验的。**

原因：我们只能观测到 $(y, w, \mathbf{x})$，即每个个体的**一个**潜在结果。我们永远无法同时观测 $y_0$ 和 $y_1$ 对同一个体。因此 $E[y_0 | \mathbf{x}, w=1]$（处理组的反事实均值）是无法从数据中直接计算的——它本身就需要可忽略性假设来识别。

**可检验的内容**：
- **重叠条件**（ATE.2）：可直接检查 $\hat{p}(\mathbf{x})$ 是否在 $(0,1)$ 内
- **协变量平衡**：在匹配或加权后，检查协变量在处理组和（加权后的）对照组之间是否平衡。若严重不平衡 → 可忽略性假设可能有问题（或重叠条件不足）
- **伪结果检验**（falsification test）：如果有处理前的 $y$ 值，可检验处理是否"预测"处理前的 $y$ 变化——如果预测，则暗示存在未观测混杂

**实践中**：可忽略性的可信度主要依赖**研究设计**（是否收集了足够丰富的处理前协变量）和**领域知识**（是否知道还有哪些重要的混杂因素未被观测），而非统计检验。

#### 2.1.6 可忽略性下 ATE 和 ATT 的可识别性总结

| 参数 | 可忽略性版本 | 重叠条件 | 识别公式 |
|------|------------|---------|---------|
| $\tau_{ate}$ | ATE.1' | ATE.2 ($0 < p(\mathbf{x}) < 1$) | $E[m_1(\mathbf{x}) - m_0(\mathbf{x})]$ |
| $\tau_{att}$ | ATT.1' |（仅需 $E[y_0\|\mathbf{x},w]=E[y_0\|\mathbf{x}]$） | ATT.2 ($p(\mathbf{x}) < 1$) | $E[m_1(\mathbf{x}) - m_0(\mathbf{x}) \mid w=1]$ |

注意 ATT 仅需要关于 $y_0$ 的均值可忽略性（不需要关于 $y_1$ 的假设）和部分重叠（仅要求对照组中存在与处理组相似的个体，不要求处理组中存在与所有对照个体相似的个体）。因此 ATT 的识别条件**弱于** ATE。

### 2.2 条件均值的识别

在均值可忽略性（ATE.1'）下：

$$E[y | \mathbf{x}, w] = \mu_0(\mathbf{x}) + w [\mu_1(\mathbf{x}) - \mu_0(\mathbf{x})] \tag{8}$$

其中 $\mu_g(\mathbf{x}) = E[y_g | \mathbf{x}]$。

**推导**：由 $y = (1-w)y_0 + w y_1 = y_0 + w(y_1 - y_0)$，取 $(\mathbf{x}, w)$ 的条件期望：

$$\begin{aligned}
E[y | \mathbf{x}, w] &= E[y_0 | \mathbf{x}, w] + w \cdot E[y_1 - y_0 | \mathbf{x}, w] \\
&= E[y_0 | \mathbf{x}] + w \cdot (E[y_1 | \mathbf{x}] - E[y_0 | \mathbf{x}]) \quad (\text{由 ATE.1'}) \\
&= \mu_0(\mathbf{x}) + w \cdot (\mu_1(\mathbf{x}) - \mu_0(\mathbf{x}))
\end{aligned}$$

令 $m_0(\mathbf{x}) = E[y|\mathbf{x}, w=0] = \mu_0(\mathbf{x})$, $m_1(\mathbf{x}) = E[y|\mathbf{x}, w=1] = \mu_1(\mathbf{x})$。注意 $m_0$ 和 $m_1$ 都是**可观测条件期望**——$m_0$ 从对照组中估计，$m_1$ 从处理组中估计。

**ATE识别**:

$$\tau_{ate}(\mathbf{x}) = m_1(\mathbf{x}) - m_0(\mathbf{x}) \tag{9}$$

$$\tau_{ate} = E[m_1(\mathbf{x}) - m_0(\mathbf{x})] \tag{10}$$

(10) 的期望是对 $\mathbf{x}$ 的总体分布取的，即 $\tau_{ate} = \int [m_1(\mathbf{x}) - m_0(\mathbf{x})] \, dF_{\mathbf{x}}(\mathbf{x})$。

**ATT识别**:

$$\tau_{att} = E[m_1(\mathbf{x}) - m_0(\mathbf{x}) \mid w = 1] \tag{12}$$

(12) 的期望是对**处理组**中 $\mathbf{x}$ 的条件分布取的，即仅对 $w=1$ 的子总体积分。

### 2.3 逆概率加权 (Inverse Propensity Score Weighting, IPW)

**倾向得分**: $p(\mathbf{x}) = P(w=1 | \mathbf{x})$

#### 2.3.1 基本矩条件的推导

IPW的核心是将观测结果按处理概率的倒数加权，从而在期望意义上恢复潜在的 $\mu_1(\mathbf{x})$ 和 $\mu_0(\mathbf{x})$。

**推导 (13)——处理组的逆概率加权**:

从观测结果 $y = y_0 + w(y_1 - y_0)$ 可知：$w y = w y_1$（只有处理组观测到 $y_1$）。利用可忽略性（$w \perp\!\!\!\perp (y_0, y_1) \mid \mathbf{x}$），计算给定 $\mathbf{x}$ 下的条件期望：

$$\begin{aligned}
E\left[\frac{w y}{p(\mathbf{x})} \;\middle|\; \mathbf{x}\right]
&= E\left[\frac{w y_1}{p(\mathbf{x})} \;\middle|\; \mathbf{x}\right] \\
&= \frac{1}{p(\mathbf{x})} \cdot E[w y_1 \mid \mathbf{x}] \\
&= \frac{1}{p(\mathbf{x})} \cdot \underbrace{E[w \mid \mathbf{x}]}_{p(\mathbf{x})} \cdot \underbrace{E[y_1 \mid \mathbf{x}]}_{\mu_1(\mathbf{x})} \quad (\text{由 } w \perp\!\!\!\perp y_1 \mid \mathbf{x}) \\
&= \mu_1(\mathbf{x})
\end{aligned}$$

这里的关键步骤在于：给定 $\mathbf{x}$ 后 $w$ 与 $y_1$ 独立，因此 $E[w y_1 | \mathbf{x}] = E[w | \mathbf{x}] \cdot E[y_1 | \mathbf{x}] = p(\mathbf{x}) \cdot \mu_1(\mathbf{x})$。

**推导 (14)——对照组的逆概率加权**:

同理，$(1-w)y = (1-w)y_0$，且 $E[1-w \mid \mathbf{x}] = 1-p(\mathbf{x})$：

$$\begin{aligned}
E\left[\frac{(1-w)y}{1-p(\mathbf{x})} \;\middle|\; \mathbf{x}\right]
&= \frac{1}{1-p(\mathbf{x})} \cdot E[(1-w)y_0 \mid \mathbf{x}] \\
&= \frac{1}{1-p(\mathbf{x})} \cdot \underbrace{E[1-w \mid \mathbf{x}]}_{1-p(\mathbf{x})} \cdot \underbrace{E[y_0 \mid \mathbf{x}]}_{\mu_0(\mathbf{x})} \\
&= \mu_0(\mathbf{x})
\end{aligned}$$

**直觉**: 逆概率加权的本质是**权重校正**——处理组中的个体以权重 $1/p(\mathbf{x})$ 重新加权，使得处理组在加权后代表整个人群（包括未被处理的人）。$p(\mathbf{x})$ 小的处理个体（不太可能被处理的类型）获得更大的权重，因为他们在处理组中稀有但在人群中确实存在。

#### 2.3.2 由 (13)(14) 推导 ATE 公式 (16)

**第一步：条件ATE的表达式 (15)**。由 (13) 和 (14) 相减：

$$\tau_{ate}(\mathbf{x}) = \mu_1(\mathbf{x}) - \mu_0(\mathbf{x}) = E\left[\frac{w y}{p(\mathbf{x})} - \frac{(1-w)y}{1-p(\mathbf{x})} \;\middle|\; \mathbf{x}\right]$$

将两项通分，注意到 $w^2 = w$（$w$ 是0-1变量），$w(1-w) = 0$：

$$\begin{aligned}
\frac{w y}{p(\mathbf{x})} - \frac{(1-w)y}{1-p(\mathbf{x})}
&= \frac{w(1-p(\mathbf{x}))y - (1-w)p(\mathbf{x})y}{p(\mathbf{x})(1-p(\mathbf{x}))} \\
&= \frac{[w - w p(\mathbf{x}) - p(\mathbf{x}) + w p(\mathbf{x})] y}{p(\mathbf{x})(1-p(\mathbf{x}))} \\
&= \frac{[w - p(\mathbf{x})] y}{p(\mathbf{x})(1-p(\mathbf{x}))}
\end{aligned}$$

因此：

$$\tau_{ate}(\mathbf{x}) = E\left[\frac{[w - p(\mathbf{x})] y}{p(\mathbf{x})[1-p(\mathbf{x})]} \;\middle|\; \mathbf{x}\right] \tag{15}$$

**第二步：无条件ATE (16)**。对 $\mathbf{x}$ 取无条件期望（迭代期望律）：

$$\tau_{ate} = E[\tau_{ate}(\mathbf{x})] = E\left[\frac{[w - p(\mathbf{x})] y}{p(\mathbf{x})[1-p(\mathbf{x})]}\right] \tag{16}$$

**解读**: $k_i = \frac{[w_i - p(\mathbf{x}_i)] y_i}{p(\mathbf{x}_i)[1-p(\mathbf{x}_i)]}$ 是一个"伪观测"——对处理组观测 $(w_i=1)$，$k_i = y_i / p(\mathbf{x}_i)$；对对照组观测 $(w_i=0)$，$k_i = -y_i / (1-p(\mathbf{x}_i))$。$k_i$ 的无条件期望恰好等于 $\tau_{ate}$。

#### 2.3.3 ATT的IPW推导 (22)

从 $y = y_0 + w(y_1 - y_0)$ 出发，考虑处理个体的权重。

**第一步**：分析 $\frac{[w - p(\mathbf{x})] y}{1-p(\mathbf{x})}$ 的分子。

利用 $y = y_0 + w(y_1 - y_0)$，有 $[w - p(\mathbf{x})] y = [w - p(\mathbf{x})] y_0 + w[w - p(\mathbf{x})] (y_1 - y_0)$。由于 $w^2 = w$：

$$w[w - p(\mathbf{x})] = w - w p(\mathbf{x}) = w(1 - p(\mathbf{x}))$$

因此：

$$\begin{aligned}
\frac{[w - p(\mathbf{x})] y}{1-p(\mathbf{x})}
&= \frac{[w - p(\mathbf{x})] y_0}{1-p(\mathbf{x})} + w(y_1 - y_0) \tag{18}
\end{aligned}$$

**第二步**：证明第一项的条件期望为零。

$$\begin{aligned}
E\left[\frac{[w - p(\mathbf{x})] y_0}{1-p(\mathbf{x})} \;\middle|\; \mathbf{x}\right]
&= E\left[\frac{[w - p(\mathbf{x})] \cdot E[y_0 \mid \mathbf{x}]}{1-p(\mathbf{x})} \;\middle|\; \mathbf{x}\right] \quad (\text{可忽略性}) \\
&= \frac{\mu_0(\mathbf{x})}{1-p(\mathbf{x})} \cdot \underbrace{E[w - p(\mathbf{x}) \mid \mathbf{x}]}_{= p(\mathbf{x}) - p(\mathbf{x}) = 0} = 0
\end{aligned}$$

**第三步**：因此，对 (18) 取条件期望后，第一项消失，仅剩第二项：

$$E\left[\frac{[w - p(\mathbf{x})] y}{1-p(\mathbf{x})} \;\middle|\; \mathbf{x}\right] = E[w(y_1 - y_0) \mid \mathbf{x}] \tag{20}$$

**第四步**：取无条件期望。

$$\begin{aligned}
E\left[\frac{[w - p(\mathbf{x})] y}{1-p(\mathbf{x})}\right]
&= E[w(y_1 - y_0)] \\
&= \underbrace{P(w=0) \cdot E[w(y_1 - y_0) \mid w=0]}_{= 0} \;+\; P(w=1) \cdot E[y_1 - y_0 \mid w=1] \\
&= P(w=1) \cdot \tau_{att}
\end{aligned}$$

设 $\Omega = P(w=1)$，则：

$$\tau_{att} = \frac{1}{\Omega} \cdot E\left[\frac{[w - p(\mathbf{x})] y}{1-p(\mathbf{x})}\right] \tag{22}$$

**IPW估计量**:

$$\hat{\tau}_{ate, ipw} = \frac{1}{N} \sum_{i=1}^{N} \left[ \frac{w_i y_i}{\hat{p}(\mathbf{x}_i)} - \frac{(1-w_i)y_i}{1-\hat{p}(\mathbf{x}_i)} \right]$$

$$\hat{\tau}_{att, ipw} = \frac{1}{N} \sum_{i=1}^{N} \frac{[w_i - \hat{\Omega}^{-1} \hat{p}(\mathbf{x}_i)] y_i}{1 - \hat{p}(\mathbf{x}_i)}, \quad \hat{\Omega} = N_1/N$$

### 2.4 回归调整 (Regression Adjustment)

分别估计 $m_0(\mathbf{x})$ (用 $w=0$ 观测) 和 $m_1(\mathbf{x})$ (用 $w=1$ 观测)：

$$\hat{\tau}_{ate, reg} = \frac{1}{N} \sum_{i=1}^{N} [\hat{m}_1(\mathbf{x}_i) - \hat{m}_0(\mathbf{x}_i)] \tag{23}$$

$$\hat{\tau}_{att, reg} = \frac{1}{N_1} \sum_{i=1}^{N} w_i [\hat{m}_1(\mathbf{x}_i) - \hat{m}_0(\mathbf{x}_i)] \tag{24}$$

**线性模型特例**: $m_0(\mathbf{x}) = \alpha_0 + \mathbf{x}\boldsymbol{\beta}_0$, $m_1(\mathbf{x}) = \alpha_1 + \mathbf{x}\boldsymbol{\beta}_1$

$$\hat{\tau}_{ate, reg} = (\hat{\alpha}_1 - \hat{\alpha}_0) + \bar{\mathbf{x}}(\hat{\boldsymbol{\beta}}_1 - \hat{\boldsymbol{\beta}}_0)$$

$$\hat{\tau}_{att, reg} = (\hat{\alpha}_1 - \hat{\alpha}_0) + \bar{\mathbf{x}}_1(\hat{\boldsymbol{\beta}}_1 - \hat{\boldsymbol{\beta}}_0)$$

$\bar{\mathbf{x}}_1$ 为处理组子样本的均值。

### 2.5 倾向得分方法 (Rosenbaum & Rubin, 1983)

倾向得分方法的突破性贡献在于：将控制**高维向量 $\mathbf{x}$** 的问题简化为控制**标量 $p(\mathbf{x})$** 的问题。这一结果由以下命题保证。

#### 命题 2.1（倾向得分的平衡性质）

在可忽略性假设 ATE.1（$w \perp\!\!\!\perp (y_0, y_1) \mid \mathbf{x}$）下，处理分配 $w$ 与潜在结果 $(y_0, y_1)$ 条件于**倾向得分 $p(\mathbf{x})$** 也是独立的：

$$w \perp\!\!\!\perp (y_0, y_1) \mid p(\mathbf{x})$$

**证明**：

**Step 1**：证明 $E[w \mid y_0, y_1, p(\mathbf{x})] = p(\mathbf{x})$。

由迭代期望律：

$$E[w \mid y_0, y_1, p(\mathbf{x})] = E\big[ E[w \mid y_0, y_1, \mathbf{x}] \;|\; y_0, y_1, p(\mathbf{x}) \big]$$

内层期望 $E[w \mid y_0, y_1, \mathbf{x}]$：由于 $w \perp\!\!\!\perp (y_0, y_1) \mid \mathbf{x}$（ATE.1），$w$ 在给定 $\mathbf{x}$ 后与 $(y_0, y_1)$ 独立，因此：

$$E[w \mid y_0, y_1, \mathbf{x}] = E[w \mid \mathbf{x}] = p(\mathbf{x})$$

这里 $E[w \mid \mathbf{x}] = P(w=1 \mid \mathbf{x}) = p(\mathbf{x})$ 就是倾向得分。

将这一结果代入外层期望：

$$E[w \mid y_0, y_1, p(\mathbf{x})] = E[p(\mathbf{x}) \mid y_0, y_1, p(\mathbf{x})] = p(\mathbf{x})$$

最后等号成立是因为 $p(\mathbf{x})$ 在给定 $p(\mathbf{x})$ 的条件下是**退化的**（已知自己的值）。换句话说，$p(\mathbf{x})$ 是 $\mathbf{x}$ 的函数，其值已经完全确定了 $p(\mathbf{x})$ 本身。

**Step 2**：由Step 1得出独立性的结论。

因为 $w$ 是二值变量，$E[w \mid y_0, y_1, p(\mathbf{x})] = p(\mathbf{x})$ 意味着：

$$P(w=1 \mid y_0, y_1, p(\mathbf{x})) = p(\mathbf{x}) = P(w=1 \mid p(\mathbf{x}))$$

$w=1$ 的条件概率不依赖于 $(y_0, y_1)$ 的取值 → $w$ 与 $(y_0, y_1)$ 条件于 $p(\mathbf{x})$ 独立。$\square$

**直观含义**：一旦你知道个体的倾向得分 $p(\mathbf{x})$，$\mathbf{x}$ 中的额外信息就不再影响处理分配的概率——倾向得分已经"汇总"了 $\mathbf{x}$ 关于处理选择的所有信息。

#### 由命题 2.1 推导 ATE

**Step 1**：用 $p(\mathbf{x})$ 替代 $\mathbf{x}$ 作为条件变量。

由命题 2.1 的独立性结果，结合观测结果分解 $y = (1-w)y_0 + w y_1$：

$$E[y \mid w, p(\mathbf{x})] = (1-w)E[y_0 \mid w, p(\mathbf{x})] + w E[y_1 \mid w, p(\mathbf{x})]$$

利用 $w \perp\!\!\!\perp (y_0, y_1) \mid p(\mathbf{x})$ → $E[y_0 \mid w, p(\mathbf{x})] = E[y_0 \mid p(\mathbf{x})]$ 且 $E[y_1 \mid w, p(\mathbf{x})] = E[y_1 \mid p(\mathbf{x})]$。因此：

$$E[y \mid w, p(\mathbf{x})] = (1-w)E[y_0 \mid p(\mathbf{x})] + w E[y_1 \mid p(\mathbf{x})] \tag{25}$$

**Step 2**：识别 $\mu_0(p)$ 和 $\mu_1(p)$。

定义 $r_0(p) = E[y \mid w=0, p(\mathbf{x})=p] = E[y_0 \mid p(\mathbf{x})=p]$
定义 $r_1(p) = E[y \mid w=1, p(\mathbf{x})=p] = E[y_1 \mid p(\mathbf{x})=p]$

$r_0(p)$ 从对照组中识别（回归 $y$ 对 $p$ 仅用 $w=0$ 的观测）
$r_1(p)$ 从处理组中识别（回归 $y$ 对 $p$ 仅用 $w=1$ 的观测）

**Step 3**：由 $r_0(p)$ 和 $r_1(p)$ 构造 ATE。

在倾向得分值 $p$ 处的条件ATE为：

$$\tau_{ate}(p) = E[y_1 - y_0 \mid p(\mathbf{x}) = p] = r_1(p) - r_0(p)$$

无条件ATE由对 $p(\mathbf{x})$ 的分布取期望得出：

$$\tau_{ate} = E[\tau_{ate}(p(\mathbf{x}))] = E[r_1(p(\mathbf{x})) - r_0(p(\mathbf{x}))]$$

用更直接的形式：

$$\tau_{ate} = E\big[E[y_1 \mid p(\mathbf{x})] - E[y_0 \mid p(\mathbf{x})]\big]$$

其中外层期望是对 $p(\mathbf{x})$ 的总体分布取的。

**Step 4**：倾向得分回归估计。

若假设 $r_0(p)$ 和 $r_1(p)$ 是 $p$ 的线性函数（$r_g(p) = \nu_g + \pi_g p$, $g=0,1$），则：

$$\tau_{ate}(p) = (\nu_1 - \nu_0) + (\pi_1 - \pi_0)p$$

$$\tau_{ate} = (\nu_1 - \nu_0) + (\pi_1 - \pi_0) \cdot E[p(\mathbf{x})]$$

等价地，可用OLS估计 $y_i$ 对 $1, w_i, \hat{p}(\mathbf{x}_i)$ 的回归——$w_i$ 的系数一致估计 $\tau_{ate}$（在线性设定下）。

**关键识别条件——重叠假设 (Overlap)**：对每个 $p(\mathbf{x}) \in (0,1)$，必须同时存在处理和对照观测。若某个倾向得分值附近只有处理组或只有对照组，则无法在该处估计 $r_1(p) - r_0(p)$——倾向得分方法失效。这就是为什么 ATE.2（$0 < p(\mathbf{x}) < 1$）是不可或缺的。

#### 方法对比：IPW vs 倾向得分回归 vs 匹配

| 方法 | 使用 $p(\mathbf{x})$ 的方式 | 核心公式 |
|------|-------------------------|---------|
| IPW | 作为**权重**的倒数 | $\tau_{ate} = E\left[\frac{(w-p)y}{p(1-p)}\right]$ |
| 倾向得分回归 | 作为**回归元** | $\tau_{ate} = E[r_1(p) - r_0(p)]$ |
| 倾向得分匹配 | 作为**匹配距离** | 对每个处理个体找 $p$ 最接近的对照个体 |

三种方法的共同基础都是 Proposition 2.1——倾向得分足以消除选择偏误。

### 2.6 匹配方法 (Matching)

对每个 $i$，补全反事实值 $\hat{y}_{i0}$ 和 $\hat{y}_{i1}$:

$$\hat{\tau}_{ate, match} = \frac{1}{N} \sum_{i=1}^{N} [\hat{y}_{i1} - \hat{y}_{i0}] \tag{26}$$

$$\hat{\tau}_{att, match} = \frac{1}{N_1} \sum_{i=1}^{N} w_i [\hat{y}_{i1} - \hat{y}_{i0}] \tag{27}$$

**对处理观测 $i$**: $\hat{y}_{i1} = y_i$, $\hat{y}_{i0} = y_{h(i)}$ ($h(i)$ 为最接近的对照观测)

**马氏距离**: $d(\mathbf{x}_h, \mathbf{x}_i) = \sqrt{(\mathbf{x}_h - \mathbf{x}_i)' \hat{\boldsymbol{\Sigma}}_x^{-1} (\mathbf{x}_h - \mathbf{x}_i)}$

基于命题2.1，可在**倾向得分**上匹配（标量 → 简单得多）。

---

## 3. 工具变量方法

### 3.1 IV估计ATE

当可忽略性不成立时使用。假设 $v_1 = v_0$（处理效应的随机部分相同）:

$$y = \delta + \tau_{ate} w + \mathbf{x}\boldsymbol{\beta}_0 + u_0 \tag{28}$$

$u_0$ 零均值，与 $(\mathbf{x}, \mathbf{z})$ 不相关；$w$ 与 $u_0$ 一般相关 → 2SLS。

**假设 ATEIV.1**:
- (a) $v_1 = v_0$
- (b) $L(v_0 | \mathbf{x}, \mathbf{z}) = L(v_0 | \mathbf{x})$: $\mathbf{z}$ 与不可观测异质性无关
- (c) $L(w | \mathbf{x}, \mathbf{z}) \neq L(w | \mathbf{x})$: $\mathbf{z}$ 有预测力

**Procedure 21.1 (两步法)**:
1. $P(w=1|\mathbf{x}, \mathbf{z}) = G(\mathbf{x}, \mathbf{z}; \boldsymbol{\gamma})$ 通过MLE → $\hat{G}_i$
2. (28) 的IV估计，IV为 $(1, \hat{G}_i, \mathbf{x}_i)$

第一阶段估计可忽略不计影响推断 → IV渐近有效（在使用 $(\mathbf{x}_i, \mathbf{z}_i)$ 函数的IV类中）。

### 3.2 放缩 $v_1 = v_0$ — 含交互项的模型

当 $v_1 \neq v_0$:

$$y = \mu_0 + \tau w + \mathbf{x}\boldsymbol{\beta}_0 + w(\mathbf{x} - \boldsymbol{\xi})\boldsymbol{\delta} + e_0 + w(e_1 - e_0) \tag{31}$$

**Procedure 21.2**:
1. 第一阶段得 $\hat{G}_i$
2. IV估计，IV为 $(1, \hat{G}_i, \mathbf{x}_i, \hat{G}_i(\mathbf{x}_i - \bar{\mathbf{x}}))$

若进一步放松 $e_1 = e_0$，仅需 $E[w(e_1 - e_0) | \mathbf{x}, \mathbf{z}] = E[w(e_1 - e_0)]$，Procedure 21.2仍一致。

### 3.3 LATE (Imbens & Angrist, 1994)

设 $\mathbf{z}$ 为二元IV。$\mathbf{z}$ 的每一取值对应反事实处理状态 $w_0, w_1$。

**四种人群**:
- $w_1=1, w_0=0$: **依从者 (Compliers)**
- $w_1=1, w_0=1$: 总是接受者 (Always-takers)
- $w_1=0, w_0=0$: 从不接受者 (Never-takers)
- $w_1=0, w_0=1$: **违抗者 (Defiers)** — 单调性排除此类

**观测处理**: $w = w_0 + z(w_1 - w_0)$

**单调性假设**: $w_1 \geq w_0$（无违抗者）

在 $\mathbf{z} \perp\!\!\!\perp (y_0, y_1, w_0, w_1)$ 和单调性下:

$$\tau_{late} = \frac{E[y | z=1] - E[y | z=0]}{E[w | z=1] - E[w | z=0]} \tag{37}$$

**Wald估计量**:

$$\hat{\tau}_{late} = \frac{\bar{y}_1 - \bar{y}_0}{\bar{w}_1 - \bar{w}_0} \tag{38}$$

LATE仅适用于**依从者**——那些受IV诱导而改变处理状态的人。

---

## 4. 断点回归设计 (RDD)

### 4.1 尖锐断点回归 (Sharp RD)

**处理规则**: $w_i = \mathbf{1}[x_i \geq c]$（$x_i$ 为驱动变量）

**关键假设**: $\mu_g(x) = E[y_g | x]$ ($g=0,1$) 在 $x=c$ 处**连续**。

**重叠完全失败**: $P(w=1|x<c)=0$, $P(w=1|x \geq c)=1$

→ RD关注断点处的处理效应:

$$\tau_c = E[y_1 - y_0 | x = c] = \mu_1(c) - \mu_0(c)$$

从 $m(x) = E[y|x] = \mathbf{1}[x<c]\mu_0(x) + \mathbf{1}[x \geq c]\mu_1(x)$:

$$m_-(c) \equiv \lim_{x \uparrow c} m(x) = \mu_0(c)$$
$$m_+(c) \equiv \lim_{x \downarrow c} m(x) = \mu_1(c)$$
$$\tau_c = m_+(c) - m_-(c)$$

**局部线性回归估计**:

$$y_i = \alpha_{0c} + \tau_c w_i + \theta_0 (x_i - c) + \delta w_i (x_i - c) + \varepsilon_i$$

仅使用 $c-h < x_i < c+h$ 的数据（$h$ 为带宽）。

### 4.2 模糊断点回归 (Fuzzy RD)

倾向得分 $F(x) = P(w=1|x)$ 在 $x=c$ 处有跳跃但不为 $0 \to 1$:

$$0 < F_+(c) - F_-(c) < 1$$

**假设**: $(y_1 - y_0) \perp\!\!\!\perp w \mid x$（个体特定收益与处理无关条件于 $x$）

从 $E[y|x] = \mu_0(x) + F(x) \cdot \tau(x)$ (39):

$$\tau_c = \frac{m_+(c) - m_-(c)}{F_+(c) - F_-(c)}$$

**Sharp RD是Fuzzy RD的特例** ($F_+(c)-F_-(c)=1$)。

**局部IV估计**（与上等价）:

$$y_i = \alpha_{0c} + \tau_c w_i + \theta_0 (x_i - c) + \delta \mathbf{1}[x_i \geq c](x_i - c) + \varepsilon_i$$

其中 $\mathbf{1}[x_i \geq c]$ 为 $w_i$ 的IV，仅用 $c-h < x_i < c+h$ 数据。

### 4.3 带宽选择

$h$ 的选择是RD实践中的核心权衡: $h$ 小 → 偏误小但方差大；$h$ 大 → 偏误大但方差小。常用方法: 交叉验证 (CV)、Imbens-Kalyanaraman最优带宽。

---

## 5. 方法对比总结

| 方法 | 识别来源 | 核心假设 | 估计的目标参数 |
|------|---------|---------|--------------|
| 回归调整/IPW/匹配 | $\mathbf{x}$ 中的可观测变量 | 可忽略性 + 重叠 | ATE / ATT |
| IV-ATE | $\mathbf{z}$ 的外生变异 | $v_1 = v_0$ + $\mathbf{z}$ 排除约束 | ATE |
| LATE (Wald) | $\mathbf{z}$ 的外生变异 | 单调性 + $\mathbf{z}$ 排除约束 | LATE (依从者) |
| Sharp RD | 断点处的局部随机化 | $\mu_g(x)$ 在断点处连续 | $\tau_c$ (断点处) |
| Fuzzy RD | 断点处的局部随机化 + IV | $\mu_0(x)$ 连续 + 个体收益独立 | $\tau_c$ (断点处) |

---

## 关键公式速查

| 公式 | 含义 |
|------|------|
| $y = y_0 + w(y_1 - y_0)$ | 观测结果分解 (4) |
| $\tau_{ate} = E[y_1 - y_0]$ | ATE定义 |
| $\tau_{att} = E[y_1 - y_0 \mid w=1]$ | ATT定义 |
| $(y_0, y_1) \perp\!\!\!\perp w \mid \mathbf{x}$ | 可忽略性 (ATE.1) |
| $p(\mathbf{x}) = P(w=1 \mid \mathbf{x})$ | 倾向得分 |
| $\tau_{ate} = E[(w-p(\mathbf{x}))y / (p(\mathbf{x})(1-p(\mathbf{x})))]$ | IPW识别 (16) |
| $\tau_{att} = E[(w-p(\mathbf{x}))y / (\Omega\cdot(1-p(\mathbf{x})))]$ | ATT的IPW (22) |
| $\hat{\tau}_{ate,reg} = \frac{1}{N}\sum[\hat{m}_1(\mathbf{x}_i)-\hat{m}_0(\mathbf{x}_i)]$ | 回归调整 (23) |
| $\tau_{late} = \frac{E[y\mid z=1]-E[y\mid z=0]}{E[w\mid z=1]-E[w\mid z=0]}$ | LATE Wald (37) |
| $\tau_c = m_+(c) - m_-(c)$ | Sharp RD |
| $\tau_c = \frac{m_+(c)-m_-(c)}{F_+(c)-F_-(c)}$ | Fuzzy RD |
