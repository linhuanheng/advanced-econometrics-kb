---
title: 按数据结构索引
type: application
domain: econometrics
created: 2026-05-19
---

# 数据结构与方法匹配

## 横截面数据 (Cross-Sectional Data)

| 问题 | 方法 |
|------|------|
| 异方差 | 稳健OLS / WLS / FGLS |
| 遗漏变量(无IV) | 代理变量 / 滞后因变量 |
| 内生性(有IV) | 2SLS / Control Function |
| 函数形式不确定 | RESET / 非嵌套检验 |
| 测量误差 | IV (CEV下OLS不一致) |
| 离群值 | LAD |
| 二元因变量 | LPM / Logit / Probit |

## 面板数据 (Panel Data)

| 问题 | 方法 |
|------|------|
| $c_i$ 与 $X$ 无关 | Pooled OLS / RE |
| $c_i$ 与 $X$ 相关 (严格外生) | FE / FD |
| $c_i$ 与 $X$ 相关 (序贯外生) | FD + IV / GMM |
| 动态面板 (AR) | Anderson-Hsiao / Arellano-Bond |
| 随机趋势 | 二次差分 / Mᵢ变换 |
| 个体特定斜率 | $M_i$ 组内变换 |
| 部分外生时间不变变量 | Hausman-Taylor |
| 聚类样本 | 聚类内去均值 |

## 混合横截面 (Pooled Cross Sections)

| 问题 | 方法 |
|------|------|
| 政策评估 (两期两组) | DID |
| 时间趋势控制 | 年份虚拟变量 |
| 偏效应随时间变化 | 时间虚拟变量 × 自变量 |

## 处理效应数据

| 机制 | 方法 |
|------|------|
| 可忽略性 (基于可观测) | 回归调整 / IPW / PSM |
| 未观测混杂 (有IV) | IV → ATE / LATE |
| 断点 (处理概率跳跃) | Sharp RD / Fuzzy RD |
