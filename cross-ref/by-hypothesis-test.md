---
title: 按假设检验索引
type: application
domain: econometrics
created: 2026-05-19
---

# 假设检验索引

## 异方差检验

| 检验 | 原假设 | 统计量 | 详见 |
|------|--------|--------|------|
| Breusch-Pagan | 同方差 | $LM = n R^2_{\hat{u}^2} \sim \chi^2_k$ | [Ch1 BP检验](../00-ols-heteroskedasticity/index.md#breusch-pagan-bp-test) |
| White检验 | 同方差 | $LM = n R^2_{\hat{u}^2} \sim \chi^2_Q$ | [Ch1 White检验](../00-ols-heteroskedasticity/index.md#white-test) |
| White特例 ($\hat{y}, \hat{y}^2$) | 同方差 | $F$ 或 $LM$ (2约束) | [Ch1 特例](../00-ols-heteroskedasticity/index.md#white检验的特例) |

## 函数形式检验

| 检验 | 原假设 | 详见 |
|------|--------|------|
| RESET | 原模型无遗漏非线性项 | [Ch2 RESET](../01-specification-data-problems/index.md#reset检验-ramsey-1969) |
| Davidson-MacKinnon | 模型A正确 (vs 模型B) | [Ch2 非嵌套检验](../01-specification-data-problems/index.md#非嵌套模型检验) |

## 内生性/过度识别检验

| 检验 | 原假设 | 统计量 | 详见 |
|------|--------|--------|------|
| Durbin-Wu-Hausman | 所有变量外生 | $\chi^2$ 二次型 | [Ch4 DWH](../03-additional-single-equation/index.md#durbin-wu-hausman-dwh-内生性检验) |
| 回归基础Hausman (CF) | $\rho_1 = 0$ (外生) | t检验 | [Ch4 CF检验](../03-additional-single-equation/index.md#基于回归的内生性检验-cf方法) |
| Sargan过度识别 | 所有IV有效 | $N R^2_u \sim \chi^2_{L_2-G_1}$ | [Ch4 Sargan](../03-additional-single-equation/index.md#过度识别约束检验-sargan-1958) |
| 第一阶段F检验 | 弱IV ($\mu_j = 0$) | F统计量 | [Ch3 弱IV](../02-instrumental-variables/index.md#5-弱工具变量问题) |

## 面板数据检验

| 检验 | 原假设 | 详见 |
|------|--------|------|
| Hausman (FE vs RE) | $\text{Cov}(X, c_i) = 0$ | [Ch5 Hausman](../04-panel-data-models/index.md#hausman检验) |
| 内生性 (CRC模型) | $\Omega_1 = 0$ (外生) | [Ch6 CRC](../05-advanced-panel-data/index.md) |

## 联合假设

| 检验 | 场景 | 详见 |
|------|------|------|
| 2SLS F检验 | 2SLS中多约束 | [Ch3 2SLS检验](../02-instrumental-variables/index.md#4-假设检验) |
| 异方差稳健Wald | 异方差下多约束 | [Ch3 稳健推断](../02-instrumental-variables/index.md#异方差稳健推断) |
