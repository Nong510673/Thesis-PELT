為方便描述與分析 UASB 系統長期運行的 COD 去除效率變化趨勢，本研究採用 R 語言套件 changepoint（Killick & Eckley, 2014）之 Pruned Exact Linear Time（PELT）演算法進行變點分析。PELT 方法基於最大概似估計（maximum likelihood estimation）框架，透過動態規劃（dynamic programming）並結合修剪（pruning）策略，有效降低計算複雜度至線性時間，適用於長時間序列之變點偵測。

本研究以均值與變異數同時變化的檢定統計量（test statistic = "Normal"）為模型假設，懲罰項（penalty）採用 Modified Bayesian Information Criterion（MBIC），以降低過度偵測之風險，並設定最小區段長度（minseglen）為 15 天，以避免過短的區段造成偵測結果不穩定。
