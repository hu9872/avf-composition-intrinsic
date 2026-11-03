# 投稿包验收清单（iScience · Manuscript_iScience_v1）

**版本**：v1.1（2026-11-03 复核）
**稿件**：`manuscript/Manuscript_iScience_v1.docx`（源：同名 `.md`）
**图例**：✅ 已核验通过；⚠️ 投稿当日需现场处理；❌ 未做（非投稿阻断项）

本版所有数字均由脚本从成品文件重新测量，不再沿用旧记录（v1.0 中"Summary 150 词""Data S1 14 sheet"等均为估值，已更正）。

---

## A. 前置材料

| # | 项目 | 刊方标准 | 实测 | 状态 |
|---|---|---|---|---|
| A1 | 标题 | 简洁、含对象 | 113 字符 | ✅ |
| A2 | Summary | 单段 ≤150 词 | **147 词**，单段 | ✅ |
| A3 | Highlights | 3–4 条，每条 ≤85 字符 | 4 条：68 / 71 / 77 / 72 字符 | ✅ |
| A4 | eTOC blurb | ≤50 词、第三人称 | **41 词**，以 `Tang et al.` 开头 | ✅ |
| A5 | Graphical abstract | 1200×1200 px、300 dpi | 1200×1200，300 dpi | ✅ |
| A6 | 作者/单位/通讯/Lead contact | 全部落实 | Yuewu Tang、Daijun Hu（\* 通讯 + † lead contact，283624940@qq.com）；重庆大学附属三峡医院肾内科 404000 | ✅ |
| A7 | 关键词 | 6 个 | 6 个 | ✅ |
| A8 | 期刊元信息残留 | 无 | 全文无旧目标期刊字样 | ✅ |
| A9 | 致谢与基金 | 无占位符 | 致谢三个公共数据集；标准"无特定基金"句 | ✅ |
| A10 | CRediT | 与基金声明自洽 | 已删 `funding acquisition`（因无基金） | ✅ |
| A11 | 封面信 | 背景 / 进展 / 意义 | 三段式齐备；`[Date]` 待投稿当日填写 | ⚠️ A11 |

## B. 正文结构

| # | 项目 | 刊方标准 | 实测 | 状态 |
|---|---|---|---|---|
| B1 | 章节顺序 | Title → Authors → Summary → Introduction → Results → Limitations → Discussion → Resource availability → Acknowledgments → Author contributions → Declaration of interests → AI 声明 → Supplemental info → Figure legends → Tables → STAR Methods → RRID → References | 完全一致 | ✅ |
| B2 | Limitations of study | 独立章节 | 独立成节 | ✅ |
| B3 | 正文篇幅 | 无硬性上限 | **约 4,490 词**（Introduction→Discussion，含 Limitations） | ✅ |
| B4 | 过强表述 | 无 overclaim | 阴性结果统一写"排除大效应"；FN-integrin 标为假设生成 | ✅ |
| B5 | 行号 / 页码 | 连续行号、页码可引 | 连续行号 + 页脚居中页码 | ✅ |
| B6 | 参考文献 | 上标、连续、含 DOI | 37 条，上标连续，DOI 齐备 | ✅ |

## C. STAR Methods 与资源

| # | 项目 | 刊方标准 | 实测 | 状态 |
|---|---|---|---|---|
| C1 | 位置 | 图注与表格之后、文献之前 | 一致 | ✅ |
| C2 | 结构 | 四标题齐全 | Resource availability / Experimental model / Method details / Quantification and statistical analysis | ✅ |
| C3 | Key Resources Table | 每项 source + identifier | 含 GEO 号、软件版本、RRID | ✅ |
| C4 | Data and code availability | GEO + 仓库 + DOI | 3 个 GEO + GitHub + concept DOI `10.5281/zenodo.23887195` | ✅ |
| C5 | 统计节集中化 | 检验/n/误差/校正/版本 | 已集中，含把握度与 TOST | ✅ |
| C6 | SGBA | 报告或说明不可得 | 已报各组性别构成并入局限 | ✅ |
| C7 | 预注册 | 有则报号，无则说明 | 明确"未正式预注册" | ✅ |

## D. 声明与合规

| # | 项目 | 刊方标准 | 实测 | 状态 |
|---|---|---|---|---|
| D1 | 伦理 | 委员会或豁免依据 | 公开去标识数据、无需审批，已写明 | ✅ |
| D2 | Declaration of interests | 独立节 + **官方表 PDF** | 章节已写；官方表需在 EM 内下载、签字后随稿上传 | ⚠️ D2 |
| D3 | 生成式 AI 声明 | 官方章节名 | 使用官方标题，说明仅辅助代码与语言 | ✅ |
| D4 | 二次分析声明 | 封面信 + 正文 | 已明确 secondary analysis、非重复发表 | ✅ |
| D5 | 查重 | 全文 <15%，单源 <3% | 内部 8-gram 比对：对 ref 8 / ref 12 均为 **0.00%**；自比 48.8%（文本回收，非他人文本）。机构 iThenticate 未跑 | ⚠️ D5 |
| D6 | 干净环境复现 | 从 raw 到图零报错 | 未整跑（约 6–8 h）；单步与断言脚本均已通过 | ❌ D6 |

## E. 图表与补充材料

| # | 项目 | 刊方标准 | 实测 | 状态 |
|---|---|---|---|---|
| E1 | 主图 | 单独文件、单页、≥300 dpi | 8 张，PNG 300 dpi + 矢量 PDF，合计 4.6 MB | ✅ |
| E2 | 图注 | 检验/n/误差/符号定义 | 8 条图注齐备 | ✅ |
| E3 | 主表 | 可编辑、无合并底纹 | 3 张（Table grid） | ✅ |
| E4 | 补充材料 | 单一 Document S1 + Data S1 | `Document_S1.docx`（Methods S1–S13、Figures S1–S3、Tables S1–S12）+ `Data_S1.xlsx`（**13 个工作表**：S1–S5、S7–S14） | ✅ |
| E5 | 补充项正文引用 | 每项至少引一次 | S1–S12、FigS1–S3、Document/Data S1 均已引用 | ✅ |

> E4 说明：工作表编号跳过 S6。成因是生成脚本先写了一版重复表再删除；**内容无缺项**，仅为编号不连续。为避免归档快照与投稿文件不一致，本次不改，如需连号可改 `44_make_supp_tables.py` 重排（会改动已归档的 Data S1）。

## F. 科学诚信

| # | 项目 | 标准 | 实测 | 状态 |
|---|---|---|---|---|
| F1 | 数字溯源 | 全部断言通过 | `90_verify_numbers.py` **48/48 PASS**（144 个表格单元格） | ✅ |
| F2 | 不可复现数字清理 | 无 | 已更正旧稿 8 处与结果不一致的数值 | ✅ |
| F3 | 敏感性分析 | 阈值/深度/门控/下限 | 把握度+TOST、来源规则 4 档、下限 3 档、门控 ±20%、等深度 | ✅ |
| F4 | 外部验证 | 独立队列 | GSE220796 反卷积不可辨识 → 改签名评分；结果不一致已如实报告 | ✅ |

## G. 投稿系统（Editorial Manager）

| # | 项目 | 现状 | 状态 |
|---|---|---|---|
| G1 | 投稿入口 | https://www.editorialmanager.com/iscience/ | ✅ 已确认 |
| G2 | 建议 / 回避审稿人 | 3–4 位建议 + ≤2 位回避 | ⚠️ 待定名单 |
| G3 | 备选联系人（Alternate contact） | 需姓名 + 邮箱 | ⚠️ 缺 Tang 的邮箱 |
| G4 | APC | USD 3,240（不含税），接受后缴 | ⚠️ 待确认预算/机构协议 |
| G5 | 官方 Declaration of Interests 表 | 需下载签字 | ⚠️ 未做 |

---

## 上传前的三条硬门槛

1. `python code/90_verify_numbers.py` → 48/48 PASS（当前已通过）。
2. 全部成品文件无 `[...]` 占位符（唯一例外：封面信 `[Date]`，投稿当日填）。
3. 上传前用 EM 生成的 PDF 校样逐页核对图、表、上标与页码。
