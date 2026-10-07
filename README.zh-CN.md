# ovasm

**用测序数据组装植物线粒体和叶绿体基因组，给出的不只是一条序列，还有证据。**

[English](README.md) | [简体中文](README.zh-CN.md)

[OrganelleVerse](https://github.com/forageseed/organelleverse) 的一部分。

<p align="center"><img src="docs/workflow.svg" alt="ovasm 流程" width="780"></p>

## 为什么用 ovasm

- **原始数据只扫一遍。** 一次挑出细胞器 reads，之后的计算都在这一小份数据上做。
- **输出带证据的组装图，而不只是序列。** 每个连接点都有支持它的 reads 数；`molecules.fasta` 是代表分子。
  reads 分不出结构时，`decisive` 会写 `false`，不会悄悄猜一个。
- **一个小小的 Rust 程序。** 不需要 conda，不需要容器，也不需要系统库（gzip 是纯 Rust 实现），
  所以 Windows、Linux、macOS 都能编译。
- **直接可用于泛基因组。** 图同时写成 OV-GFA，`ovasm pan` 可以把多个样本合并。

## 安装

```bash
git clone https://github.com/forageseed/ovasm.git
cd ovasm
cargo build --release            # 需要 Rust 1.80 或更高版本
./target/release/ovasm --version
```

## 使用

需要全基因组测序数据（HiFi；Illumina 加 `--read-type sr`），以及近缘物种的细胞器基因组作为种子。
两个细胞器的种子一起给，这样各自招募各自的 reads：

```bash
ovasm run --reads sample.hifi.fastq.gz --organelle both \
  --seeds mitochondrion=mt_reference.fasta --seeds plastid=pt_reference.fasta \
  --out out/ --threads 8
```

结果写在 `out/mitochondrion/` 和 `out/plastid/`（只拼一个细胞器时，`--organelle mitochondrion` 或
`plastid` 直接写进 `out/`）：

| 文件 | 内容 |
|---|---|
| `assembly.gfa` | 组装图，片段带深度 |
| `evidence.json` | 每个连接和分叉有多少 reads 支持；哪些片段是重复 |
| `molecules.fasta` | 代表分子，沿 reads 支持的路径解开得到 |
| `linearization.json` | 代表序列怎么选的，reads 能不能分出结构（`decisive`） |
| `organelle.gfa` | OV-GFA 格式的图，可直接给 `ovasm pan` |
| `summary.json` | 运行摘要：用的 k、`k_accepted`、分子长度、连接数、文件列表 |

用结果之前先看 `summary.json`（`k_accepted`、`decisive`）。

没有参考序列时，`--recruit discover` 会尝试按 k-mer 深度找出细胞器 reads（只支持 HiFi）。
它要求细胞器 reads 的深度明显高于核背景，可能失败；有种子时请给种子。
[OrganelleVerse](https://github.com/forageseed/organelleverse) 的 Python 包自带一套现成的种子库。

## 命令

`ovasm run` 把下面这些阶段串起来；每个阶段也都是单独的命令（`ovasm <命令> --help`）。

| 命令 | 作用 |
|---|---|
| `run` | 完整流程：招募、组装、证据、解环、统一格式 |
| `recruit` | 用种子从全基因组数据里招募细胞器 reads |
| `discover` | 不用种子，按 k-mer 深度找细胞器 reads |
| `panel` | 按保守蛋白判断一份 reads 来自哪个细胞器 |
| `correct` | 含噪长读长（ONT、CLR）自纠错 |
| `assemble` | 建组装图 |
| `evidence` | 统计图上每个连接的 reads 支持 |
| `linearize` | 选出代表序列 |
| `unify` | 把图改写成 OV-GFA |
| `pan` | 把多个样本的 OV-GFA 合并成一个泛基因组图 |
| `components` | 只保留属于目标细胞器的连通部分 |
| `identify` | 找组装结果最接近的参考基因组 |
| `label` | 按基因组成给图的片段标注线粒体或质体 |

## 状态

已在金标准基因组（日本晴、拟南芥 Col-0、丹参）上测试。更多物种的测试仍在进行，
所以在新物种上的结果请先核对。

## 在 OrganelleVerse 里使用

```python
ov.assembly.assemble(data, organelle="mitochondrion", method="ovasm")
```

OrganelleVerse 通过环境变量 `ORG_VERSE_OVASM_BIN` 或 `PATH` 找到这个程序。

## 许可证

MIT
