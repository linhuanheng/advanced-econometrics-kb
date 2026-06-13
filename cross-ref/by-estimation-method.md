---
title: 按估计方法索引
type: application
domain: econometrics
created: 2026-05-19
---

# 估计方法索引

## OLS 及相关

| 方法 | 使用场景 | 详见 |
|------|---------|------|
| OLS | 基准线性回归 | [异方差](../00-ols-heteroskedasticity/index.md) |
| 异方差稳健OLS | 异方差但无内生性 | [异方差稳健推断](../00-ols-heteroskedasticity/index.md#2-异方差稳健推断-heteroskedasticity-robust-inference) |
| WLS | 已知异方差形式 | [WLS/GLS](../00-ols-heteroskedasticity/index.md#4-加权最小二乘-wlsgls) |
| FGLS | 未知异方差形式 | [FGLS](../00-ols-heteroskedasticity/index.md#可行gls-fgls) |
| LAD | 离群值稳健 | [LAD](../01-specification-data-problems/index.md#离群值-outliers) |

## IV 及 2SLS

| 方法 | 使用场景 | 详见 |
|------|---------|------|
| IV | 单一内生变量+恰好识别 | [IV估计](../02-instrumental-variables/index.md#1-iv的基本原理) |
| 2SLS | 多IV/过度识别 | [2SLS](../02-instrumental-variables/index.md#2-两阶段最小二乘-2sls) |
| Control Function | 线性/非线性内生性 | [CF方法](../03-additional-single-equation/index.md#2-控制函数方法-control-function-approach) |
| 2SLS with Generated IV | IV需第一阶估计 | [Generated Instruments](../03-additional-single-equation/index.md#2sls-with-generated-instruments) |

## 面板数据

| 方法 | 使用场景 | 详见 |
|------|---------|------|
| Pooled OLS | $c_i$ 与 $X$ 无关 | [面板数据](../04-panel-data-models/index.md) |
| Random Effects (RE) | $c_i$ 与 $X$ 无关，效率优先 | [RE](../04-panel-data-models/index.md#随机效应-random-effects-re) |
| Fixed Effects (FE) | $c_i$ 与 $X$ 相关 | [FE](../04-panel-data-models/index.md#固定效应-fixed-effects-fe) |
| First Differencing (FD) | 序贯外生/弱工具变量 | [FD](../04-panel-data-models/index.md#一阶差分-first-differencing-fd) |
| IV-FD (Anderson-Hsiao) | 动态面板 | [动态面板](../05-advanced-panel-data/index.md#2-动态面板模型-dynamic-panel-models) |
| Chamberlain System OLS | FE框架下的投影方法 | [Chamberlain](../05-advanced-panel-data/index.md#4-chamberlain方法) |
| Hausman-Taylor | 部分外生时间不变变量 | [HT模型](../05-advanced-panel-data/index.md#5-hausman-taylor模型) |

## 政策评估

| 方法 | 使用场景 | 详见 |
|------|---------|------|
| DID | 自然实验 (两期两群) | [DID](../03-additional-single-equation/index.md#5-双重差分-difference-in-differences-did) |
| 回归调整 | 可忽略性假设 | [回归调整](../07-treatment-evaluation/index.md#回归调整-regression-adjustment) |
| 倾向得分匹配 | 可忽略性假设 | [PSM](../07-treatment-evaluation/index.md#倾向得分匹配-propensity-score-matching) |
| IPW | 可忽略性假设 | [IPW](../07-treatment-evaluation/index.md#逆概率加权-inverse-propensity-score-weighting-ipw) |
| IV-ATE/LATE | 未观测混杂 | [IV方法](../07-treatment-evaluation/index.md#3-工具变量方法) |
| Sharp RD | 断点处的因果效应 | [SRD](../07-treatment-evaluation/index.md#尖锐断点回归-sharp-rd) |
| Fuzzy RD | 断点+IV | [FRD](../07-treatment-evaluation/index.md#模糊断点回归-fuzzy-rd) |
