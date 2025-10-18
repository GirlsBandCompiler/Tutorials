# 等式约束消除

在多面体分析中, 由于大部分循环变量都只取整数，我们会经常用到只含有整数的矩阵.

我们知道, 多面体的定义是 $P = \{ \mathbf{x} \mid A\mathbf{x} \leq \mathbf{b} \}$, 这也被称作多面体的 [H-表示](polyhedron.md#def-1). 这里我们认为 $\mathbf{x}$ 的取值范围是 $\mathbb{N}^k$, 也就是全为整数的向量.

但是, 在进行[依赖分析](../analysis/pluto.md#计算循环中的依赖)的时候, 我们产生了一些等式. 也就是说, 我们的集合变为了 $P = \{ \mathbf{x} \mid A\mathbf{x} = \mathbf{b}, C\mathbf{x}\leq \mathbf{d} \}$. 为了使用 [Barvinok 算法](barvinok.md), 我们需要消除 $A\mathbf{x} = \mathbf{b}$ 的这些约束, 将 $P$ 转化为只有不等式的形态.

为了方便, 我们设 $A$ 是 $m\times n$ 的矩阵, 也就是 $\mathbf{x}$ 有 $m$ 个元素, 一共有 $n$ 个等式.

你可能会想问, 为什么不能直接把 $A\mathbf{x} = \mathbf{b}$ 的每一行等式 $\mathbf{a}\cdot \mathbf{x} = \mathbf{b}$ 转化为两个不等式 $\mathbf{a}\cdot \mathbf{x} \leq \mathbf{b}$ 和 $\mathbf{a}\cdot \mathbf{x} \geq \mathbf{b}$ 呢?

这是因为, $\mathbf{x}$ 在等式约束下, 实际上只有 $m - \mathrm{rank}(A)$ 个自由的维度. 如果将等式改为两个不等式, $\mathbf{x}$ "看起来"就有了 $m$ 个自由维度, 这会对 Barvinok 算法产生干扰. Barvinok 算法有时会将原本的多面体拆成几个小的多面体, 然后将其中的整点个数相加或相减; 额外的自由维度会导致这些"小"多面体中含有无穷多个整点, 相减时就会出现问题.

我们知道, 如果 $\{\mathbf{v}_1, \dots, \mathbf{v}_n\}$ 是 $A$ 零空间的一组基, 而 $\mathbf{x}_0$ 是 $A\mathbf{x} = \mathbf{b}$ 的一个解, 那么 $A\mathbf{x} = \mathbf{b}$ 的所有解都可以写成 $\mathbf{x} = V\mathbf{y} + \mathbf{x}_0$ 的形式, 其中 $V$ 是以 $\mathbf{v}_i$ 为列向量构成的矩阵, 而 $\mathbf{y}$ 则是任意取定的参数. 接下来, 只需要将它代入剩下的不等式, 就完全消除了等式约束:

$$P = \{\mathbf{x} \mid C(V\mathbf{y}+\mathbf{x}_0) \leq \mathbf{d}\}$$

现在唯一的问题, 就是我们需要保证 $V$, $\mathbf{y}$ 和 $\mathbf{x}_0$ 都是整数, 才能使用 Barvinok 算法.

## 幺模矩阵

**定义 1.** 对于整数方阵 $A$, 如果 $\det(A) = \pm 1$, 那么 $A$ 是**幺模矩阵** (*unimodular matrix*).

设 $A$ 是幺模的. 根据克拉默法则, $A\mathbf{x} = \mathbf{b}$ 的解是

$$x_i = \frac{\det(A_i)}{\det(A)}$$

其中 $A_i$ 是将 $A$ 的第 $i$ 列换成 $b$ 的结果.

当 $A$ 和 $\mathbf{b}$ 都只包含整数时, 由于 $\det(A) = \pm 1$, 我们有 $x_i$ 也必定为整数. 这说明方程 $A\mathbf{x} = \mathbf{b}$ 必有整数解. 这正是我们想要的性质.

不过, 一般而言, 我们遇到的约束 $A\mathbf{x} = \mathbf{b}$ 中, $A$ 都不是幺模的; $A$ 甚至大概率不是方阵. 所以, 我们需要一种方法将 $A$ 分解为幺模矩阵.

## Hermite 标准型

对于矩阵 $A$, 我们可以找到一个幺模矩阵 $U$, 使得 $H = UA$, 而且 $H$ 是行阶梯型的上三角阵. 这时, $H$ 被称为 Hermite 标准型.

多面体库应当可以自动计算 Hermite 标准型. 如果不行的话, 一种简单的实现方法是高斯消元, 不过消元时需要通过欧几里得算法保证所有数据都是整数.

高斯消元的过程只包含"交换行/列"或是"将一列乘以一个数并加到另一列"两种操作. 这些操作的行列式都是 $\pm 1$, 所以最终产生的矩阵 $U$ 是幺模的.

**性质 1.** 设 $A^T$ 的 Hermite 标准型是 $H = UA^T$, 那么 $U$ 的最后 $n - \mathrm{rank}(H)$ 行就是 $A$ 零空间的基.

**证明.** $UA^T = H$, 那么 $AU^T = H^T$. $H$ 是行阶梯的, 所以 $H^T$ 的最后 $n - \mathrm{rank}(H)$ 列都是零, 所以 $U^T$ 的最后 $n - \mathrm{rank}(H)$ 列在 $A$ 的零空间中. 考虑到 $U$ 是满秩的 (因为 $\det(U)\neq 0$), 这些列向量线性无关. 又考虑到 $A$ 的零空间一共就 $n - \mathrm{rank}(H)$ 维 (因为 $\mathrm{rank}(H) = \mathrm{rank}(A)$), 所以这些列向量是一个极大线性无关组, 也就是一组基.

现在, 只要找到一组特解 $A\mathbf{x}_0 = \mathbf{b}$, 就可以完成了.

## Smith 标准型

在 Hermite 标准型上更进一步, 我们可以找到两个幺模矩阵 $U, V$, 使得 $D = UAV$, 其中 $D$ 是对角阵. 此时 $D$ 被称为 Smith 标准型.

构造 Smith 标准型的方法和 Hermite 标准型是类似的, 也是通过高斯消元完成.

对于方程 $A\mathbf{x}_0 = \mathbf{b}$, 可以将其直接写作 $V^{-1}DU^{-1}\mathbf{x}_0 = \mathbf{b}$. 从而 $D(U^{-1}\mathbf{x}_0) = V\mathbf{b}$. 考虑到 $D$ 是对角阵, 我们很容易解出 $U^{-1}\mathbf{x}_0$. 当然, 这一步可能产生非整数解, 比如约束是 $2x_1 + 2x_2 = 5$; 如果出现了这种情况, 那么原来的多面体就没有任何整点.

如果 $U^{-1}\mathbf{x}_0 = \mathbf{y}$ 是整数向量, 那么我们就获得了 $\mathbf{x}_0 = U\mathbf{y}$ 这个特解. 至此, 我们就消除了所有的等式约束.
