# Flexibility-score-tool-of-amino-acids
A tool which can score residues in protein based on R.
参考不依赖三维结构信息，被广泛应用于蛋白质序列分析、无序区域预测以及结构功能研究中的Karplus-Schulz氨基酸柔性评价体系。将20种标准氨基酸划分为刚性（A，L，H，V，Y，I，F，C，W，M）和柔性（K，S，G，P，D，E，Q，T，N，R）两种，通过统计目标残基两侧刚性氨基酸的邻接数目（0/1/2），结合高分辨率X射线晶体结构数据库中的归一化主链温度因子（B-factor）参数（BNORM0/BNORM1/BNORM2），采用滑动窗口加权算法（权重系数：0.25，0.50，0.75，1.00，0.75，0.50，0.25），逐位计算序列中各氨基酸的柔性值。然而，经典方法存在局限性：序列首尾端残基因邻域信息的缺失导致柔性预测存在偏差；缺乏对序列整体柔性水平的量化评估。本研究编写交互式分析工具，剔除序列首尾各3个氨基酸的局部柔性值，消除边界效应干扰，引入有效残基（N ≥ 7时，第4至N - 3位）柔性值的算术均值：
Mean Flexibility Score = 1/(N-6) ∑_(i=4)^(N-3)▒〖F_i （N≥7 ）〗

Citation：Karplus PA, Schulz GE. Prediction of chain flexibility in proteins[J]. Naturwissenschaften. 1985, 72(4):212-213.
