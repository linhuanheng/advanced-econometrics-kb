---
title: 高级计量经济学知识库 (Graduate Econometrics III)
type: framework
domain: econometrics
source: IESR, Jinan University — Graduate Econometrics III Lecture Notes (Spring 2025)
textbooks:
  - "Wooldridge, J. M., Introductory Econometrics: A Modern Approach"
  - "Wooldridge, J., Econometric Analysis of Cross-Section and Panel Data, MIT Press, 2010"
  - "A. Colin Cameron and Pravin K. Trivedi, Microeconometrics: Methods and Applications"
created: 2026-05-19
---

# 高级计量经济学知识库

## 方法论总览

本知识库涵盖高级计量经济学（Graduate Econometrics III）的核心内容，按方法论逻辑链组织：

```
OLS回归基础 (heteroskedasticity)
    ↓
设定与数据问题 (specification, measurement error, proxy variables)
    ↓
内生性处理 (IV, 2SLS, Control Function)
    ↓
单方程扩展 (generated regressors, CRC, Diff-in-Diff)
    ↓
面板数据模型 (RE, FE, FD, Hausman-Taylor)
    ↓
高级面板专题 (sequential exogeneity, dynamic panels)
    ↓
非线性模型 (binary response: Logit/Probit)
    ↓
处理效应评估 (ATE, propensity score, IV, RDD)
```

## 知识领域索引

| 编号 | 领域 | 核心问题 |
|------|------|---------|
| 00 | [OLS回归与异方差](00-ols-heteroskedasticity/index.md) | 异方差下的稳健推断与加权最小二乘 |
| 01 | [设定与数据问题](01-specification-data-problems/index.md) | 模型误设、代理变量、测量误差 |
| 02 | [工具变量估计](02-instrumental-variables/index.md) | 内生性的IV/2SLS解决方案 |
| 03 | [单方程扩展专题](03-additional-single-equation/index.md) | 控制函数、随机系数、DID |
| 04 | [面板数据模型](04-panel-data-models/index.md) | 非观测效应、RE/FE/FD估计 |
| 05 | [高级面板专题](05-advanced-panel-data/index.md) | 序贯外生、动态面板、Hausman-Taylor |
| 06 | [二元响应模型](06-binary-response-models/index.md) | LPM、Logit、Probit |
| 07 | [处理效应评估](07-treatment-evaluation/index.md) | ATE、倾向得分、IV、断点回归 |

## 核心方法论决策树

### 1. 选择估计方法

```
是否存在内生性？
├── 否 → OLS (检查异方差 → robust SE / WLS)
└── 是 → 是否有有效IV？
    ├── 是 → 2SLS / Control Function
    └── 否 → 是否有面板数据？
        ├── 是 → FE/FD (消除时间不变非观测效应)
        └── 否 → 代理变量 / 滞后因变量
```

### 2. 面板数据方法选择

```
cᵢ 是否与 xᵢₜ 相关？
├── 否 → Random Effects (GLS, 更有效)
└── 是 → Fixed Effects / First Differencing
    └── 严格外生是否成立？
        ├── 是 → FE (T>2时偏好) / FD
        └── 否 → 序贯外生? → IV-FD / GMM
```

### 3. 处理效应评估

```
处理分配机制？
├── 基于可观测变量 (ignorability)
│   ├── 回归调整
│   ├── 倾向得分匹配
│   └── 逆概率加权
├── 基于不可观测变量 (selection on unobservables)
│   ├── IV → ATE / LATE
│   └── Control Function
└── 基于断点 (discontinuity)
    ├── Sharp RD
    └── Fuzzy RD (IV within RD)
```

## 交叉引用索引

- [按估计方法索引](cross-ref/by-estimation-method.md)
- [按假设检验索引](cross-ref/by-hypothesis-test.md)
- [按数据结构索引](cross-ref/by-data-structure.md)
