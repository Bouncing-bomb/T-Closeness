# T-Closeness 的实现与分析

## 一、实验基本信息

- **实验题目**：T-Closeness 的实现与分析
- **实验内容**：构建满足 T-Closeness 的 Mondrian 算法，并结合信息损失（GCP）与不同 $t$ 值进行实验分析
- **实验目标**：
  1. 完成 Earth Mover's Distance（EMD）的计算；
  2. 在基础 Mondrian 实现的基础上，构建满足 T-Closeness 的 Mondrian 算法；
  3. 结合信息损失（GCP）与不同 $t$ 值进行实验分析。
- **评价标准**：补充代码 50% + 实验报告 50%
- **提交内容**：实验报告 PDF + 代码压缩包
- **截止时间**：11.10 23:59（两周后）

---

## 二、前置知识

### 1. 隐私攻击模型

- **同质化攻击（Homogeneity Attack）**：如果一个等价类中的所有记录在敏感属性上的取值非常相似，攻击者就可以通过该等价类推测出其中所有记录的敏感信息。k-anonymity 没有考虑敏感属性的多样性，因此可能出现这种情况。
- **背景攻击（Background Knowledge Attack）**：攻击者拥有额外的背景信息，可以将其与其他信息源结合，从而推测出个体的敏感信息。k-anonymity 只关心记录能否与其他记录“混淆”，没有考虑敏感属性的取值分布，因此容易受到背景攻击。
- **相似性攻击（Similarity Attack）**：l-diversity 要求等价类中敏感属性至少有 $l$ 个“不同”的取值，但不关心这些取值在语义上是否相似。例如“胃溃疡、胃炎、胃癌”虽然满足多样性要求，但语义相近，攻击者仍可推断出大致病情。因此 l-diversity 无法抵御相似性攻击。
- **偏斜攻击（Skewness Attack）**：l-diversity 不比较等价类内分布与全局分布的偏离程度。当等价类内的分布高度偏斜时，攻击者仍能以很高的概率猜中敏感值。

### 2. T-Closeness 简介

- **定义**：每个等价类（equivalence class）内敏感属性的分布应当“接近”其在整个原始数据库中的分布。设某等价类内敏感属性的分布为 $P$，全表（全局）分布为 $Q$，则要求

$$
D(P,Q)\le t,
$$

其中 $D$ 为下面定义的 EMD 距离。

- **距离度量：Earth Mover's Distance（EMD）**：设敏感属性共有 $m$ 个不同取值，$p_i$、$q_i$ 分别为分布 $P$、$Q$ 在第 $i$ 个取值上的概率，$d_{ij}$ 为第 $i$ 个取值与第 $j$ 个取值之间的 ground distance，$f_{ij}$ 为把“土方”从 $i$ 搬运到 $j$ 的流量，则

$$
D[P,Q]=\min_{F}\ \mathrm{WORK}(P,Q,F)=\min_{F}\sum_{i=1}^{m}\sum_{j=1}^{m}d_{ij}f_{ij},
$$

满足如下约束：

$$
\begin{aligned}
& f_{ij}\ge 0,\quad 1\le i\le m,\ 1\le j\le m,\\
& p_i-\sum_{j=1}^{m}f_{ij}+\sum_{j=1}^{m}f_{ji}=q_i,\quad 1\le i\le m,\\
& \sum_{i=1}^{m}\sum_{j=1}^{m}f_{ij}=\sum_{i=1}^{m}p_i=\sum_{i=1}^{m}q_i=1.
\end{aligned}
$$

### 3. Mondrian 划分技术

- **划分（Partitioning）**：
  - **单属性（数值型）**：以中位数划分。例如 $k=2$、数据集为 $[1,2,3,3,4,5]$ 时，严格划分得到 $[1,2,3,3]$ 与 $[4,5]$；松弛划分得到 $[1,2,3]$ 与 $[3,4,5]$（位于中位数位置上的并列值被平均分配到左右两侧）。
  - **多属性（数值型）**：每个 partition 默认选择当前归一化跨度最大的属性进行划分，直到无法继续划分为止。

- **合法切分条件**：在位置 $x$ 处垂直于坐标轴 $X$ 切分，当且仅当切分后左右两个子集的记录数均不少于 $k$，且两个子集各自的敏感属性分布与全集分布的 EMD 距离均不超过 $t$。

- **算法（严格多维划分）**：

```text
Anonymize(partition)
  if (no allowable multidimensional cut for partition)
    return φ : partition → summary
  else
    dim ← choose_dimension()
    fs ← frequency_set(partition, dim)
    splitVal ← find_median(fs)
    lhs ← {t ∈ partition : t.dim ≤ splitVal}
    rhs ← {t ∈ partition : t.dim > splitVal}
    return Anonymize(rhs) ∪ Anonymize(lhs)
```

其中 “no allowable multidimensional cut” 指该 partition 的大小不足 $2k$，或不存在同时满足上述 $k$ 约束与 $t$-closeness 约束的切分。

- **算法（宽松/松弛多维划分）**：落在中位数位置上的数据点被平均分配到左右子集。

### 4. 评价指标

- **NCP（归一化确定性惩罚，Normalized Certainty Penalty）**：

$$
NCP_{A_i}(t)=\frac{z_i-y_i}{|A_i|},
$$

其中 $t$ 为一条记录，$[y_i,z_i]$ 为记录 $t$ 在属性 $A_i$ 上被泛化后所处的区间，$|A_i|$ 为该属性在整个数据表上的取值跨度（即 $\max-\min$）。

- **GCP（全局确定性惩罚，Global Certainty Penalty）**：

$$
GCP(T)=\frac{\sum_{t\in T}\sum_{i=1}^{n}NCP_{A_i}(t)}{|T|\cdot n},
$$

其中 $n$ 为准标识符属性的个数，$|T|$ 为记录总数。$GCP=0$ 表示没有发生泛化，$GCP=1$ 表示所有准标识符都被泛化到最顶层（全表取值范围）。

### 5. 数据集

- **Adult Dataset**：共 14 个属性，其中 6 个数值属性、8 个类别属性，另有 1 个类别标签 `income`。
- 下载链接：<https://archive.ics.uci.edu/ml/datasets/adult>

---

## 三、实验要求

### 1. 补充代码（50%）

- 正确实现全集分布计算（10 分）
- 正确实现 EMD 距离计算（10 分）
- 正确实现满足 T-Closeness 的 Mondrian 算法（20 分）
- 代码清晰、有注释（10 分）

### 2. 实验报告（50%）

- 问题描述（5 分）、解决思路（5 分）、算法具体流程图（10 分）
- 分析 $k$、`qi`、$t$、`datasize` 对 GCP、max EMD、RunTime、等价类数量的影响（15 分）
- 实验问题：T-Closeness 条件对 Mondrian 算法的隐私保护效果有什么影响？k-anonymity 与 t-closeness 这两种隐私规则在影响实验效果的表征上有什么异同点？（10 分）
- 总结满足 T-Closeness 的 Mondrian 算法的特点（5 分）

---

## 四、实验说明

### 1. 全局分布计算（Part 1）

- 对整个数据集中的敏感属性进行分布统计，并储存在全局变量 `GLOBAL_SA_DIST` 中。
- 例如，敏感属性集为 `[A, A, B, C, C, C]`，则 `GLOBAL_SA_DIST = {'A': 2/6, 'B': 1/6, 'C': 3/6}`。
- 按照说明补充 `init(data, k, QI_num=-1)` 中 `GLOBAL_SA_DIST` 的计算。

### 2. EMD 距离计算（Part 2）

- 对某个划分的分布与全局分布计算 EMD 距离：

$$
EMD(pv,qv)=\frac{1}{2}\sum_{i}\left|pv_i-qv_i\right|.
$$

- 说明：该式成立的前提是敏感属性为类别型，且任意两个不同取值之间的 ground distance 均取 1（此时 EMD 等于总变差距离）。若敏感属性为有序或数值型，则应改用累积分布形式：

$$
EMD(pv,qv)=\sum_{i}\left|F_{pv}(i)-F_{qv}(i)\right|,
$$

其中 $F_{pv}$、$F_{qv}$ 分别为 $pv$、$qv$ 的累积分布函数。

- 按照说明补充 `compute_emd(partition)`。

### 3. 满足 T-Closeness 的 Mondrian（Part 3 & 4）

- 只需完成面向数值属性的基础型 Mondrian。
- 需要在 k-anonymity 检查之后增加 t-closeness 检查。
- 可使用预设函数 `satisfy_t_closeness(partition, t)`。
- 按照说明补充 `anonymize_strict()` 与 `anonymize_relaxed()`。若某个划分不满足 t-closeness，则拒绝该划分，将该维度标记为不可用，并继续尝试其他维度。

### 4. 文件结构

```text
T-Closeness-code
├── data
│   ├── adult.data          # 原始数据
│   ├── adult_x.txt         # 类别属性的泛化层次（每行格式：子节点;父节点）
│   └── anonymized.data     # 输出的匿名后数据
├── utils
│   ├── read_adult_data.py  # 读取数据
│   ├── utility.py          # 比较函数
│   ├── gentree.py          # 树结构
│   └── numrange.py         # 数值范围
├── mondrian_tcloseness.py  # （待补充）满足 t-closeness 的 Mondrian
└── anonymizer.py           # 不同参数的使用方法
```

### 5. 关键变量

- `QI_LEN`：定义的准标识符个数。
- `QI_ORDER[i]`：第 $i$ 个准标识符对应的排序后属性值序列。
- `QI_RANGE[i]`：第 $i$ 个准标识符的属性值范围大小。
- `QI_DICT[i][qi_value]`：第 $i$ 个准标识符的属性值 `qi_value` 在 `QI_ORDER[i]` 中的索引，即 `QI_ORDER[i][index] = qi_value` 时 `QI_DICT[i][qi_value] = index`。
- `RESULT`：存储每个划分段。
- `T_THRESHOLD`：全局变量 $t$ 值。
- `GLOBAL_SA_DIST`：数据集的原始（全局）分布。

### 6. 运行方法

- 全局变量 `T_THRESHOLD` 在代码中修改。
- 命令：`python anonymizer.py [r|s] [k|qi|data|k-value]`
  - `r`：宽松划分模式；`s`：严格划分模式；
  - `k`：改变参数 $k$ 进行实验；`qi`：改变准标识符个数进行实验；`data`：改变数据集大小进行实验；
  - `k-value`（具体数值）：使用指定的 $k$ 值运行。

---

## 五、附录：本次修订说明

1. EMD 定义补充了 $\min_{F}$，并将流量守恒式第二项的下标由 $\sum_{i=1}^{m}f_{ji}$ 修正为 $\sum_{j=1}^{m}f_{ji}$；同时补充 $m$、$d_{ij}$ 的含义。
2. t-closeness 定义的记号由 $D(P_s,P_T)\le t$ 统一为 $D(P,Q)\le t$，并说明 $P$ 为等价类内分布、$Q$ 为全局分布。
3. EMD 简式 $\frac{1}{2}\sum_i|pv_i-qv_i|$ 补充了适用前提（类别型敏感属性、ground distance 取 1）；有序或数值型敏感属性给出累积分布形式。
4. 原“Generalization（泛化）”小节实际描述的是划分（partitioning），标题改为“划分（Partitioning）”；松弛划分补充“并列值平均分配”的说明。
5. NCP 记号由 $NCP_{A_i}\{t\}$ 改为 $NCP_{A_i}(t)$，并补全 $y_i$、$z_i$、$|A_i|$ 的含义；GCP 中明确 $n$ 为准标识符属性个数。
6. “合法切分条件”明确为左右两个子集分别满足记录数不少于 $k$ 且各自与全集的 EMD 距离不超过 $t$。
7. Adult 数据集属性描述更正为 6 个数值属性 + 8 个类别属性 + 1 个类别标签。
8. “全集分布”统一为“全局分布”，与 `GLOBAL_SA_DIST` 的命名保持一致。
9. 严格划分算法补充了 “no allowable multidimensional cut” 判定条件的注释。
