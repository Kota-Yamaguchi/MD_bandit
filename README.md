# MD_bandit

[環境構築]
conda create -n [仮想環境名] python=3.6 
source activate [仮想環境名]

conda install -c conda-forge mdanalysis

仮想環境のMDAnalysisのanalysisの中にpca_mod.pyを設置する。

・MDAnalysisの場所の確認方法
  import MDAnalysis
  
  print(MDAnalysis.__file__)
  


[使用方法]

python bandit.py -s {gro} -t {xtc} -r {gro}　-s {bandit reaction coordinate}
-s トポロジーファイルにあたるもの,GRO形式で
-t 解析したいトラジェクトリ
-r ターゲットとなるGROファイル
-s バンディットスコアにしたい反応座標を指定する


[説明]

PCA,RMSD、ContactMapの３つの反応座標のどれが一番良いかをバンディットアルゴリズムで探索する

PaCSスコアを計算　→ preバンディットスコアを計算　→ バンディットスコアを計算

バンディットスコア　＝ preバンディットスコア - 1サイクル前のpreバンディットスコア　+ 忘却率 * 今までのバンディットスコア
