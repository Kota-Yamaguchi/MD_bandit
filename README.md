# MD_bandit

[環境構築]
conda create -n [仮想環境名] python=3.6 
source activate [仮想環境名]

conda install -c conda-forge mdanalysis

[使用方法]

python bandit.py -s {gro} -t {xtc} -r {gro}
-s トポロジーファイルにあたるもの,GRO形式で
-t 解析したいトラジェクトリ
-r ターゲットとなるGROファイル



[説明]

PCA,RMSD、ContactMapの３つの反応座標のどれが一番良いかをバンディットアルゴリズムで探索する

現在PCAを引くとエラーが発生する。
