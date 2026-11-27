<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type='text/x-mathjax-config'>
	MathJax.Hub.Config({
		tex2jax: {
			inlineMath: [['$', '$']],
			dsplayMath: [ ['$$','$$'],["\\[","\\]"] ]
		}
	});
</script>

# TMTGEM v1.3.1 マニュアル（日本語版）
## Takuto Minami’s Tsunami-Generated Electro-Magnetic simulation code

#### 2026年11月27日 時点

南 拓人 (Takuto MINAMI)
tminami@port.kobe-u.ac.jp
（神戸大学 理学研究科）

※ 本書は英語版マニュアル TMTGEM_manual_v1_3.md の日本語訳です。内容に差異がある場合は英語版を優先してください。

### 目次

1. [はじめに](#Introduction)
2. [重要な注意事項](#ImportantNotes)
3. [必要な環境](#Required_environments)
4. [サンプルシミュレーションの実行](#Run_sample_simulations)  
 - サンプル1: [東北](#Tohoku_small)
 - サンプル2: [イースター島](#Easter)
 - サンプル3: [チリ](#Chile)
 - サンプル4: [地中海](#Mediterranean)

5. [コードの詳細](#Details_of_code)  
5.1. [メッシュ生成部](#Meshgeneration)  
5.2.  [時間領域シミュレーション部](#Time-domain-simulation)
6. [TMTGEM の出力ファイル](#Output_files)
7. [リリースノート](#Release_notes)
8. [質問・コメント歓迎！！](#Questions)
9. [参考文献](#References)



## 1 はじめに <a id="Introduction"></a>
TMTGEM は、海洋中の起電力を介して津波が生成する電磁場を計算するシミュレーションコードです。TMTGEM は Fortran で書かれており、スパース直接法ソルバー PARDISO を利用するために Intel Math Kernel Library (MKL) を使用します。TMTGEM は PARDISO による OpenMP 並列計算が可能で、1回のシミュレーションで行列分解が1回だけで済むため、時間領域の計算に適しています。
TMTGEM は 1: メッシュ生成部、2: シミュレーション部 の2つから構成されています。TMTGEM は現在のところ [COMCOT 津波シミュレーションコード](http://223.4.213.26/archive/tsunami/cornell/comcot_down.htm) の出力を用いることを前提としています。サンプルコードには、[Kawashima and Toh (2016)](#KT2016) による修正を加えた COMCOT ver 1.7 が含まれています。
TMTGEM はメッシュ生成に [Gmsh](https://www.soest.hawaii.edu/gmt/) と、[Oishi et al. (2013)](#O2013) に記載の押し出し (extrude) アルゴリズムを用いています。押し出し用サブルーチンの大部分は、オープンソースの流体シミュレーションコード [Fluidity](http://fluidityproject.github.io/) に由来します。
	

## 2 重要な注意事項 <a id="ImportantNotes"></a>
座標系: TMTGEM では X: 東向き、Y: 北向き、Z: 上向き
解いている方程式: Minami et al. (2017, JGR) を参照してください
TMTGEM の地形データ: 経度 [deg]、緯度 [deg]、標高 [m, 下向き正] を仮定

## 3 必要な環境 <a id="Required_environments"></a>
以下のパッケージを PC にインストールしてください。
メッシュ生成用の [Gmsh](http://gmsh.info/)
MKL ライブラリ付きの Intel Fortran コンパイラ（“ifort –mkl=parallel ***” が使えることを確認してください）
シミュレーション結果の表示用に GMT ([generic mapping tool](https://www.soest.hawaii.edu/gmt/))
ghostscript（“gv” コマンド用のみ。例: TMTGEM/Tohoku/flow の plot_z.sh）

### 必要なデータセット
海底地形データ: https://topex.ucsd.edu/WWW_html/mar_topo.html から xyz データをダウンロードしてください。
ただし、COMCOT と TMTGEM では下向き正の水深データを仮定しています。
標高の符号を変換するには、例えば TMTGEM/Tohoku/topo/mktopo.sh を利用してください。
（任意）電気伝導度データ: 海水の電気伝導度データをダウンロードして利用できます
（任意）背景地球磁場データ

## 4 サンプルコードの実行 <a id="Run_sample_simulations"></a>
とにかく TMTGEM を動かしてみましょう！！
TMTGEM/ のホームフォルダにある Tohoku/、Chile/、Mediterranean/、Easter/ の4つのフォルダがサンプルシミュレーション用です。おおよその実行時間を表1に示します。

表1. サンプルシミュレーションの概要と所要時間

TMTGEM/

|コマンド	| Tohoku/    |Tohoku/  | Tohoku/ | Easter/ | Chile/ | Mediterranean/ | 
|:----------:|:-----:|:----:|:----:|:----:|:----:|:----:|
|run_comcot.sh <br> 時間 |flow/ <br> 4m30s| - | - | flow/ <br> 5m30s| flow/ <br> 8m46s  | flow/ <br> 8m7s|
|tetmeshgen.sh <br> 時間 | mesh/ <br> 2m15s  | - | -  | mesh/ <br> 1m43s | mesh/ <br> 5m10s | mesh/ <br> 3m21s |
|run.sh <br>  時間	     | em/ <br> 7m52s |	em_IGRF/ <br> 8m38s	| em_woa/ <br> 8m32s	| em/ <br> 7m50s	| em/ <br> 7m57s| em/ <br> 15m14s |

（“時間” は計算時間で、2.8 GHz Intel Core i7・16 GB RAM の MacBook Pro で計測しました。./run.sh は 4 スレッドの OpenMP 並列計算で実行しました。“NX*NY” は comcot 計算の格子構成です。“# of lines” は生成された 3-D メッシュ em3d.msh の辺の数を百万 (M) 単位で示したものです。）

## 4.1 サンプル1: 東北 <a id="Tohoku_small"></a>

### 4.1.1 Tohoku/em/ でサンプルシミュレーションを実行するための準備 <a id="Preparation"></a>
津波電磁場シミュレーションには、準備として津波の流速データとシミュレーション用メッシュが必要です。ターミナルで以下のコマンドを実行して準備してください。

#### 津波シミュレーションと TMTGEM メッシュ生成用の *.xyz 地形ファイルの準備
	$ cd TMTGEM/Tohoku/topo
	$ ./mktopo.sh                       (topo.xyz が生成される)
#### comcot の実行
	$ cd TMTGEM/Tohoku/flow
    $ ./run_comcot.sh                   (津波流速データを生成)
#### TMTGEM 用メッシュの生成
	$ ./plot_z.sh 000600                (津波波高を確認、図1a 参照)
	$ cd TMTGEM/Tohoku/mesh
	$ ./tetmeshgen.sh                   (メッシュ生成)
#### 生成したメッシュの確認
	$ gmsh polygonki.msh                 (生成した 2 次元メッシュを確認、図1b 参照)
	$ gmsh em3d.msh                      (生成した 3 次元メッシュを確認、図1c 参照)

準備が成功すると、以下の図を確認できます。
![Fig1](./images/Fig1.png)
図1. (a) Tohoku/flow/ で “./plot_z.sh 000600” により描いた津波波高。(b) Tohoku/mesh/ で “gmsh polygonki.msh” により描いた polygonki.msh の 2 次元メッシュ。(c) Tohoku/mesh/ で “gmsh em3d.msh” により描いた最終的な 3 次元メッシュ。

### 4.1.2 東北サンプルでの TMTGEM の実行
比較的小さな東北メッシュ “Tohoku/mesh/em3d.msh” を用いたサンプルシミュレーションが3種類あります。

##### Tohoku/em/
 一様な背景磁場と、与えられた海水電気伝導度を使用
##### Tohoku/em_IGRF/
 IGRF の背景磁場と、一様な海水電気伝導度を使用
##### Tohoku/em_woa/
 一様な背景磁場と、不均質な海水電気伝導度を使用

上記3つのサンプルは、t = 5 分という短い計算期間ではほぼ同じ結果になります。いくつかの種類の制御ファイル（例: Tohoku/em/tmtgem.ctl）で TMTGEM が正しく動くことを確認してください。各フォルダで以下のコマンドを実行してください。

	$ cd TMTGEM/Tohoku/em             (Tohoku/em の場合)
	$ ./run.sh                        (TMTGEM シミュレーションを実行)
### シミュレーションの実行
	$ ./plotgeomag.sh                 (使用した背景磁場を確認、図2)
	$ ./plotcond.sh                   (使用した海水電気伝導度を確認、図3)
	$ ./bxyzplot.sh 000060 3          (it = 60 の第3成分 (Bz) を確認、図4)
	
サンプルシミュレーションが正常に終了すると、以下の図を確認できます。

 
![Fig2a](./images/Fig2a.png)
図2a. 使用した地磁気主磁場 (Fx, Fy, Fz) = (0,0,-35000) [nT]（シミュレーション後に Tohoku/em/ または Tohoku/em_woa で ./plotgeomag.sh により描いた geomag.pdf）


 
![Fig2b](./images/Fig2b.png)
図2b. IGRF 主磁場を用いた場合（Tohoku/em_IGRF で ./plotgeomag.sh により描いた geomag.pdf）。

![Fig3a](./images/Fig3a.png)
図3a. 海面に接する四面体要素 (surfcond.out) と海底に接する四面体要素 (bottcond.out) における海水電気伝導度の分布図。Tohoku/em/ または Tohoku/em_IGRF で ./plotcond.sh により描いたもの。
 
![Fig3b](./images/Fig3b.png)
図3b. 海面に接する四面体要素（上段）と海底に接する四面体要素（下段）における海水電気伝導度の分布図。シミュレーション後に Tohoku/em_woa で ./plotcond.sh により描いたもの。左列は、元の 1 度 × 1 度格子の海水電気伝導度データ woa13_decav81B0_C00an01.csv (Tyler et al., 2017) から計算した分布です（“_wocomp” は “w/o complementation（補間なし）” の意味）。右列は TMTGEM シミュレーション用に補間した海水電気伝導度分布です。左列では、データのないビンは既定の海水電気伝導度 3.3 S/m で埋めています。

![Fig4](./images/Fig4.png)

図4. シミュレーション後に Tohoku/em で “./bxyzplot.sh 000060 3” を実行した結果。黒のコンターは水平粒子速度の等振幅線です。it = 000060 では、Tohoku/em_IGRF と Tohoku/em_woa でもほぼ同じ結果になります。

[](![Fig6](./images/Fig6.png))

<!--(Fig. 6 Vertical component of tsunami-generated magnetic field, drawn by “./bxyzplot.sh 000060 3” (left) and the surface conductivity drawn by ./plotcond.sh in Tohoku/em_large.) -->

[]( ![Fig7](./images/Fig7.png) ) 

<!--(Figure 7. The upward component of IGRF magnetic field, drawn by ./plotgeomag.sh in Tohoku/em_large. Note that IGRF fields are calculated only in the area for COMCOT calculation.)-->

## 4.2 サンプル2: イースター島 <a id="Easter"></a>
TMTGEM は外洋の孤立した島も扱うことができます。サンプルシミュレーションの実行手順は 4.1 節と同じです。以下の3つのコマンドを実行してください。
TMTGEM/Easter/flow/ で

	$./run_comcot.sh

TMTGEM/Easter/mesh/ で

	$./tetmeshgen.sh

TMTGEM/Easter/em/ で

	$./run.sh

結果は以下のとおりです。

![Fig8](./images/Fig8.png)
図8. tetmeshgen.sh で生成した polygonki.msh。Easter/mesh/ で $gmsh polygonki.msh により描いたもの。

![Fig9](./images/Fig9.png)
図9. comcot シミュレーションによる海面変位の結果。Easter/flow で $./plot_z.sh 000600 により描いたもの。
 
![Fig10](./images/Fig10.png)
図10. (a) “./bxyaplot.sh 000060 3” により描いた津波磁場の鉛直成分。黒のコンターは |vh|、緑のコンターは水深 2800 m を示します。(b) “./plotgeomag.sh” により描いた IGRF モデルの Fz。(c) “./ixyhplot_xy.sh 000060” により描いた、海洋中の鉛直積分した全電流。黒のコンターは (a) と同じで、矢印は電流の向きを示します。(d) “./plotcond.sh” により描いた海面の電気伝導度。(a)–(d) の図は TMTGEM/Easter/em で作成しました。


## 4.3 サンプル3: チリ <a id="Chile"></a>
TMTGEM は、地形領域の右上の角が陸上にある場合も扱えるようになりました。これまでのサンプルと同様に、以下の手順で結果を得ることができます。

TMTGEM/Chile/flow/ で

	$./run_comcot.sh

TMTGEM/Chile/mesh/ で

	$./tetmeshgen.sh

TMTGEM/Chile/em/ で

	$./run.sh

結果は以下のとおりです。

![Fig11](./images/Fig11.png) 

図11. tetmeshgen.sh で生成した polygonki.msh。Chile/mesh/ で “gmsh polygonki.msh” と “gmsh em3d.msh” により描いたもの。

![Fig12](./images/Fig12.png)

 図12. comcot シミュレーションの結果。Chile/flow で
 $./plot_z.sh 000300
 および
 $./plot_z.sh 000600 により描いたもの。

![Fig13](./images/Fig13.png)
 
図13. 図10 と同様、ただし Chile/em でのシミュレーション。(a)〜(d) は、“./run.sh” の終了後にそれぞれ “./bxyzplot.sh 000060 3”、“./plotgeomag.sh”、“./ixyhplot_xy.sh 000060”、“plotcond.sh” を実行してください。


## 4.4 サンプル4: 地中海 <a id="Mediterranean"></a>
	改良された TMTGEM のメッシュ生成部は、多くの島や半島が計算境界で切られる地中海のような非常に複雑な海底地形の領域でも安定して動作します。

TMTGEM/Mediterranean/flow/ で

	$./run_comcot.sh

TMTGEM/Mediterranean/mesh/ で

	$./tetmeshgen.sh

TMTGEM/Mediterranean/em/ で

	$./run.sh

結果は以下のとおりです。

![Fig14](./images/Fig14.png) 

図14. ./tetmeshgen.sh で生成したメッシュ。Mediterranean/mesh で以下により描いたもの。

	$gmsh polygonki.msh 

および

	$gmsh em3d.msh
 
![Fig15](./images/Fig15.png)

図15. Mediterranean/flow/ で “./plot_z.sh 000600” により描いた comcot シミュレーションの結果。
 
![Fig16](./images/Fig16.png)

図16. Mediterranean/em の結果。図の形式は図10 と同じです。(b) と (d) では、comcot の計算領域に対応する黒い矩形内の値のみが妥当である点に注意してください。
 “./bxyaplot.sh 000060 3” により描いた津波起因磁場の鉛直成分（左）、以下により描いた IGRF の Fz

	$./plotgeomag.sh （中央）
	
および、Mediterranean/em で以下により描いた海面の電気伝導度

	$./plotcond.sh

## 5. コードの詳細 <a id="Details_of_code"></a>
## 5.1 メッシュ生成部 <a id="Meshgeneration"></a>
TMTGEM で最も重要かつ難しいのは、津波電磁場シミュレーションに適した四面体メッシュを生成する部分です。メッシュ生成のおおまかな流れは以下のとおりで、すべて各メッシュ生成フォルダ（例: Tohoku/mesh）の tetmeshgen.sh に含まれています。

2 次元三角形メッシュを生成する
海底地形データを読み込み、三角形メッシュに反映する
海面の三角形を海底まで押し出して、海洋の 3 次元メッシュ (mesh1) を生成する
大気と固体地球の 3 次元四面体メッシュ (mesh2) を生成する
mesh1 と mesh2 を結合して、3 次元電磁場シミュレーション用のメッシュとする

流れの理解には Minami et al. (2017) の Fig. 1 を参照してください。

< tetmeshgen.sh による主な出力ファイル >
*.geo : ジオメトリファイル（gmsh で表示）
polygonki.msh : 地表面/海面の 2 次元三角形メッシュ（gmsh で表示）
ocean.msh : 海洋層のみの 3 次元四面体メッシュ（gmsh で表示）
pre3d.msh : 大気と固体地球のみの 3 次元四面体メッシュ（gmsh で表示）
em3d.msh : 津波電磁場シミュレーション用の最終的な 3 次元四面体メッシュ（gmsh で表示）
lineinfo.dat : em3d.msh の辺の情報

出力メッシュの表示（各メッシュファイルは gmsh で表示してください）

![Fig17](./images/Fig17.png)

図17. tetmeshgen.sh の出力例。

メッシュ生成の制御パラメータはすべて制御ファイル（例: Tohoku/mesh/mesh.ctl）に含まれています。mesh.ctl の各パラメータの意味は以下を参照してください。
### < mesh.ctl >
    # of header        |3
    ## This control file is generated on Sep 7, 2016, (see m_param_mesh.f90)
    ## It is confirmed that this control file works well
    --------10--------20--------------------------------
    topofile           |../topo/W120E176S20N60_1min.xyz 入力地形ファイル
    west bound  [deg]  | 140.0   図8 参照
    east bound  [deg]  | 148.0   図8 参照
    south bound [deg]  | 34.0    図8 参照
    north bound [deg]  | 42.0    図8 参照
    nmcal lower bnd[km]| -500.0  図10 参照
    nmcal upper bnd[km]| 1000.0  図10 参照
    sizein  [km]       |20.0     図8 参照
    sizeout [km]       |100.0    図8 参照
    sizecoastratio [km]|0.3      図9 参照
    lenout  [km]       |300.0    図8 参照
    # of observatories |2        細かいメッシュが必要な観測点の数
    name of obs        |ESA      観測点名
    lon,lat[deg]alt[km]|141.35472       39.236944      0.0（Alt は使用しない）
    name of obs        |b14
    lon,lat[deg]alt[km]|144.807595      39.0582       -5.380（書式 g15.7）
    reso at obs [km]   |1.0      図8 参照
    sigma forGauss [km]|50.0     図8 参照
    0:none,1:vres<hres |0        とりあえず 0 に設定
    # of max layers    |3        図11 参照
    nlayer=1 up to [km]|-0.3     図11 参照
    nlayer=2 up to [km]|-2.0     図11 参照
    header             |polygon  変更可能だが “polygon” が安全
    bgmesh: upzin  [km]|20.0
    bgmesh: downzin[km]|-20.0
    bgmesh:sizein3d[km]|14.0        
    以下は bgmesh3d.pos 用のパラメータで、現在は使用していません。

そのまま設定してください。

## < 2 次元三角形メッシュの制御方法 >
 
polygonki.msh
図8. メッシュ生成パラメータの意味。
	海岸線に沿ったメッシュサイズの制御方法（“coastratio” を使用）
 
polygonki.msh
図9. 海岸線に沿った解像度に関するパラメータの意味。


	em3d.msh の上端・下端境界の設定方法（“nmcal upper bnd” と “nmcal lower bnd” を使用）

 
pre3d.geo
図10. 3 次元メッシュの構造。



	水深に応じた海洋中の層数の制御方法（“# of max layers” と “nlayer=* up to [km]” を使用）
 
< em3d.msh の鉛直断面 >
図11. 海洋中の層数を制御するパラメータ。

注: 背景メッシュ “bgmesh3d.pos” は現在使用できません。

	pre3d.msh において地表面・海面・海底の三角形はそのまま保持される（v1.3.1 以降）

gmsh で pre3d.geo（大気と固体地球）の 3 次元メッシュを生成する際、これまでは bgmesh3d.pos に従って海面・海底・地表面の三角形が細分化されていました。この細分化で追加された節点は ocean.msh には存在しないため、em3d.msh が海洋層と大気・固体地球との境界面で不整合（non-conforming）になっていました。v1.3.1 以降、tmtgem/mesh/mk3dgeo.f90 は pre3d.geo に以下の transfinite 制約を追加します。

	Transfinite Line{1:nline} = 2;
	Transfinite Surface{1:ntriall};

ここで nline は 2 次元三角形の辺の数、ntriall は 2 次元三角形の数（pre3d.geo の面 1〜ntriall）です。これらの制約により、gmsh はこれらの辺や三角形の内部に節点を追加しないため、pre3d.msh の面上の三角形は ocean.msh のものと一致します。外側の箱の 16 本の辺と 10 枚の面には制約をかけていないため、従来どおり gmsh によってメッシュが切られます。



## 5.2 時間領域シミュレーション部 <a id="Time-domain-simulation"></a>
	シミュレーション部の制御パラメータ（例: Tohoku/em_large/tohoku.ctl）は以下のとおりです。

< TMTGEM/Tohoku/em_large/tohoku.ctl >

	--------10--------20--------------------
	1:err,2:fwd,3:both |2  1: div 誤差計算, 2: 通常の順計算
	3d em mesh file    |../mesh_large/em3d.msh       <- メッシュ情報
	line infro file    |../mesh_large/lineinfo.dat   <- メッシュ情報
	mesh control file  |../mesh_large/mesh.ctl       <- メッシュ情報
	ocean mesh file    |../mesh_large/ocean.msh      <- メッシュ情報
	start time [sec]   |0.0                          <- 開始時刻
	time interval [sec]|5.0                          <- 電磁場シミュレーションの時間間隔
	end time   [sec]   |300.0                        <- 電磁場シミュレーションの終了時刻
	1:COMCOT,2:analytic|1                            <- 入力流速場が COMCOT か解析解か
	xgrid file         |../flow/layer01_x.dat        <- COMCOT の経度グリッドファイル
	ygrid file         |../flow/layer01_y.dat        <- COMCOT の緯度グリッドファイル
	depth grid  file   |../flow/layer01.dat          <- COMCOT の水深グリッドファイル
	header for m file  |../flow/m_01_                <- COMCOT 入力ファイルのヘッダ
	header for n file  |../flow/n_01_                <- COMCOT 入力ファイルのヘッダ
	t int COMCOT [sec] |5.0                          <- COMCOT 入力の時間間隔
	1:vh,2:vz,3:vh+vz  |1                            <- 入力流速に用いる成分
	0:IGRF,1:FIL,2:cons|0                            <- 背景磁場の選択
	Decimal year       |2011.0                       <- IGRF の場合、時刻を10進年で指定
	xgrid file         |../flow/layer01_x.dat        <- IGRF を与える経度グリッド
	ygrid file         |../flow/layer01_y.dat        <- IGRF を与える緯度グリッド
	out folder  (bxyz) |./bxyz/                      <- 磁場の出力フォルダ
	out folder  (exyz) |./exyz/                      <- 電場の出力フォルダ
	out folder  (vxyz) |./vxyz/                      <- 流速場の出力フォルダ
	# of observatories |2                            <- 時系列を出力する観測点の数
	1:lonlath,2:xyz[km]|1                            <- 経度緯度か xyz 座標かの選択
	obs name           |ESA                          <- 観測点名
	obs position [km]  |141.35472      39.236944        0.0           <- 位置 3g15.7
	obs name           |b14
	obs position [km]  |144.807595     39.0582          -5.380        <- 位置 3g15.7
	0:no,1:nxny,2:seafl|1   平面図オプション 1: 格子, 2: 全三角形
	nx, ny             |400         400                               <- x, y 方向の格子点数
	time interval [sec]|60.0                                          <- 平面図出力の時間間隔
	1:same dpth,2:seafl|2                            <- 平面図の z 1: 同一深さ, 2: 海底
	nx,ny-> kimsh file |../mesh_large/polygonki.msh  <- 平面図用のメッシュ情報
	nx,ny->ki23dptr fl |../mesh_large/ki23dptr.dat   <- 平面図用のメッシュ情報
	# of cond from top |3      
	air cond   [S/m]   |1.e-8 大気の電気伝導度
	ocean cond [S/m]   |3.3   海水の電気伝導度（既定値）
	crust cond [S/m]   |0.01  地殻/マントルの電気伝導度
	1:no,2:ocecond file|2     1: 一様, 2: ファイルから読み込み
	need woa csv file  |../../woa/woa13_decav81B0_C00an01.csv 海水の電気伝導度


一部のオプションは動作しない可能性があります。挙動がおかしい場合はソースプログラムを確認してください。

# 5.3 TMTGEM の出力ファイル <a id="Output_files"></a>
1)  ### "bxyz/[Site Name]_bxyz_ts.dat"
	“# of observatories” > 0 とし、“obs name” と観測点位置を設定すると、既定の出力フォルダ設定では磁場と電場の時系列がそれぞれ ./bxyz/ と ./exyz/ に出力されます。磁場時系列ファイルは以下の4列形式です。
	Time [sec],  bx [nT], by [nT], bz [nT] 
	ここで bx, by, bz はそれぞれ津波起因磁場の東向き、北向き、上向き成分です。

2) ### "exyz/[Site Name]_exyz_ts.dat"
	bxyz/[Site Name]_bxyz_ts.dat と同様に、TMTGEM の実行により各観測点の電場も得られます。10列形式は以下のとおりです。

	Time [sec], [e の x,y,z 成分], [v×F の x,y,z 成分], [e+v×F の x,y,z 成分]
 
	ここで e は誘導電場、v×F は誘導起電力、e+v×F は観測機器で観測される全電場です。値の単位はすべて [μV/m] です。x,y,z 成分の向きは [Site Name]_bxyz_ts.dat と同じです。

3) ### "bxyz/bxyz_xy2D******.dat"  (coord_xy2D.dat とともに)  
	“plan view option” を 1 または 2 にすると、計算した磁場変動の平面図用の出力ファイルが得られます。bxyz_xy2D*****.dat という名前のファイルが “time interval [sec]” ごとに bxyz/ フォルダに保存されます。****** は反復ステップを表し、t = ****** × “time interval [sec]” のスナップショットに対応します。ファイルは以下の3列形式です。

	bx [nT], by [nT], bz[nT]

	各行の位置はファイル “bxyz/coord_xy2D.dat” に対応しており、その形式は以下のとおりです。

	格子点数 (nx * ny)

	X [km], Y [km], Z [km] 

	ここで X, Y, Z は東向き距離、北向き距離、および海面を基準とした標高です。

4) ### "exyz/ixyh_xy2D******.dat" (coord_xy2D.dat とともに)
	海洋層で鉛直積分した全電流密度 $$i_H h=∫σ[e_H+(v×F)_H]dz$$ の出力ファイルです。値の単位は [μA\/m^2×km]=[A/km] です。2列形式は以下のとおりです。

	$ih_x$ (A/km), $ih_y$ (A/km)

	ファイルの行数は、bxyz/bxyz_xy2D******.dat と同様に nx と ny に基づき nx * ny です。ih_x と ih_y はそれぞれ i_H h の東向き成分と北向き成分です。



# 6. リリースノート <a id="Release_notes"></a>
### Version 1.3.1 (2026年11月27日)
tmtgem/mesh/mk3dgeo.f90 を修正し、gmsh で pre3d.geo の 3 次元メッシュを生成する際に、海面・海底・地表面の三角形がそれ以上細分化されないようにしました。pre3d.geo に "Transfinite Line{1:nline} = 2;" と "Transfinite Surface{1:ntriall};" を追加することで、ocean.msh に存在しない節点が gmsh により追加されるのを防ぎ、em3d.msh が海洋との境界面で整合（conforming）するようにしました。5.1 節の「pre3d.msh において地表面・海面・海底の三角形はそのまま保持される」を参照してください。

### Version 1.3 (2021年12月10日)
solver/m_param.f90 に新しいパラメータ “nalyer” を追加しました。1回のシミュレーションで comcot の複数の層を使用できます。制御ファイルの12行目で層の数を変更することで、comcot シミュレーションのネスト格子の層を増やすことができます。サンプルシミュレーションは Tohoku/em_2layer を参照してください。

### Version 1.2.3 (2021年7月26日)
pre3d.msh の 2 つ以上のボリュームを反映して em3d.msh を生成できるよう、tmtgem/mesh/combine3d.f90 を修正しました。

### Version 1.2.2 (2021年5月31日)
mk3dgeo.f90 のバグを修正し、[2] make n3d と [3] calzkiinland を適切な順序に並べ替えました。Tohoku/ のシェルスクリプトの一部を gmt 6 に対応させました。

### Version 1.2.1 (2021年5月24日)
Easter、Chile、Mediterranean のサンプルフォルダにある gmt を含むシェルスクリプトを GMT6 用に更新しました。

### Version 1.2 (2020年8月5日)
2020年8月5日、TIARES の例において tmtgem/mesh/m_intersection.f90 と m_coastline.f90 のバグを修正しました。

### Version 1.1
ver 1.0 から ver 1.1 では、図22 に示すような交差を修正する新しいモジュールを組み込みました。現在、この種の交差は tmtgem/mesh/coastline.exe が交差モジュール（mesh/m_intersection.f90 参照）により自動的に回避します。
ver 1.0 から ver 1.1 で、ソルバー部 tmtgem/solver/* に大きな変更はありません。

![Fig22](./images/Fig22.png)

図22. polygonki.geo における交差の発生例（“gmsh polygonki.geo” で表示）。

# 7 質問・コメント歓迎！！<a id="Questions"></a>
TMTGEM の実行中に問題が生じた場合は、takuto.minami126@gmail.com までご連絡ください。まだ多くのバグが残っていると思われ、本マニュアルも十分ではありません。皆様からの質問やコメントを反映してマニュアルを改訂することで、他のユーザーの助けにもなります。

# 8 参考文献 <a id="References"> </a>

- Kawashima, I., & Toh, H. (2016). Tsunami-generated magnetic fields may constrain focal mechanisms of earthquakes. Scientific reports, 6, 28603.<a id="KT2016"></a>

- Minami, T., Toh, H., Ichihara, H., & Kawashima, I. (2017). Three‐Dimensional Time Domain Simulation of Tsunami‐Generated Electromagnetic Fields: Application to the 2011 Tohoku Earthquake Tsunami. Journal of Geophysical Research: Solid Earth, 122(12), 9559-9579. <a id="M2017"></a>

- Oishi, Y., Piggott, M. D., Maeda, T., Kramer, S. C., Collins, G. S., Tsushima, H., & Furumura, T. (2013). Three‐dimensional tsunami propagation simulations using an unstructured mesh finite element model. Journal of Geophysical Research: Solid Earth, 118(6), 2998-3018.<a id="O2013"></a>

- Tyler, R. H., Boyer, T. P., Minami, T., Zweng, M. M., & Reagan, J. R. (2017). Electrical conductivity of the global ocean. Earth, Planets and Space, 69(1), 156. <a id="T2017"></a>
 

# 付録A: 支配方程式の離散化形式の導出

支配方程式

$$ \nabla \times \nabla \times \mathbf{A} + \mu \sigma \frac{\partial \mathbf{A}}{\partial t} = \mu \sigma (\mathbf{v} \times \mathbf{B}_0)
$$ 

に対し、ガラーキン法による有限要素法では、次式を数値的に解きます。

$$ \int_\Omega \left[ w\cdot \left( \nabla \times \nabla \times \mathbf{A} + \mu \sigma \frac{\partial \mathbf{A}}{\partial t}\right)-\mu \sigma (\mathbf{v}\times \mathbf{B}_0)\right]dV=0 $$
ここで w は任意の試験関数で、ガラーキン法では w に辺要素の基底関数を採用します。Ω は計算領域です。
関係式 $$\nabla \cdot ( \mathbf{A}\times \mathbf{B} ) = \mathbf{B}⋅(\nabla \times \mathbf{A} ) -\mathbf{A}\cdot \left( \nabla \times \mathbf{B} \right)$$を用いると、式 (1) の左辺第1項は次のように変形されます。
$$ \int_\Omega { \mathbf{w} \cdot ( ∇ × ∇ × \mathbf{A} ) dV } = \int_\Omega { ∇ \cdot \left( ∇ \times \mathbf{A} \right) × \mathbf{w} } + \left( ∇ × \mathbf{w} \right) \cdot (∇×A) dV $$
$$ = \int_\Omega {∇⋅((∇×A)×w)+(∇×w)⋅(∇×A)}dV$$
$$=\int_{\partial \Omega } \left[ ((∇×\mathbf{A})×\mathbf{w})\cdot \mathbf{n}dS \right] + \int_\Omega {(∇ × \mathbf{w}) \cdot (∇ × \mathbf{A}) }dV $$
$$= \int_{\partial \Omega} [(B×w)⋅ndS] +  ∫_\Omega  ( \nabla × \mathbf{w}) \cdot  (\nabla \times \mathbf{A} ) dV ⋅⋅⋅(2) $$
ここで $B=∇×A $ は磁場、$∂Ω$ は計算領域の境界面、$n$ は境界面 ∂Ω の外向き単位法線ベクトルです。∂Ω において
$n⋅(B×w)=w⋅(n×B)=0$ とできることを考慮します。
ここで $∇⋅(n×B)=B⋅(∇×n)-n⋅(∇×B)=-n⋅(μj)=0$ は境界面を貫く電流がないことを示します。すると、式 (2) は次のようになります。
$$ ∫_Ω[w⋅(∇×∇×A)dV]=∫_Ω{(∇×w)⋅(∇×A)}dV⋅⋅⋅(3)$$
式 (3) を式 (1) に代入すると、次式を得ます。
$$ \int_\Omega {(\nabla \times \mathbf{w})\cdot (\nabla \times \mathbf{A})}dV+∫_Ω[μσw⋅∂A/∂t dV]=μσ\int_\Omega\left[ w⋅(v\times B_0)dV \right].⋅⋅⋅(4)$$
有限要素法では、計算領域を有限個の四面体に分割し、各四面体要素について体積積分を計算します。以下では e 番目の要素 Ω_e の体積積分のみを考えると、式 (4) は e 番目の要素について次のように書き直されます。
$$ ∫_{Ω_e} {(∇×\mathbf{w}_i )⋅(∇×A_e )}dV+∫_{Ω_e} \left[μσ\mathbf{w}_i⋅\frac{\partial A_e}{\partial t} dV \right] =μσ∫_{Ω_e}[\mathbf{w}_i⋅(v×B_0)dV]⋅⋅⋅(5)$$
ここで $\mathbf{w}_i  (i=1,⋅⋅⋅,6)$ は e 番目の四面体の各辺に対する辺要素基底関数です。e 番目の要素内のベクトルポテンシャルを次のように表します。
$$ A_e=∑_{j=1}^6 \left[ \mathbf{w}_j  [Al]_j \right]⋅⋅⋅⋅(6)$$
ここで
$$ \mathbf{w}_j=w_{kl}=λ_k ∇λ_l-λ_l ∇λ_k $$
$λ_k$ は節点基底関数です。$j$ と $k$, $l$ の関係を以下の表に示します。
 
#### 表A1. 辺番号 $i$, $j$ に対する節点番号 
|i, j (辺番号) |	k|	l|	m|	n |
|--|:--:|--  |--| -- |
|1|	1|	2|	3|	4|
|2|	2|	3|	1|	4|
|3|	1|	3|	4|	2|
|4|	1|	4|	2|	3|
|5|	2|	4|	3|	1|
|6|	3|	4|	1|	2|
||||||||

表A1
式 (6) を用いると、式 (5) は次のように変形できます。

$$ ∑_{j=1}^6 \left[ ∫_{Ω_e} {(∇×\mathbf{w}_i )⋅(∇×\mathbf{w}_j )}dV [Al]_j \right] + ∑_{j=1}^6 \left[ ∫_{Ω_e}[μσ\mathbf{w}_i⋅\mathbf{w}_j] dV \frac{d}{dt} [Al]_j \right] = μσ_e ∫_{Ω_e} [\mathbf{w}_i⋅(v×B_0)dV]⋅⋅⋅(7)$$

これを行列形式で書き直すと、

$$[M_e ]+μσ_e [N_e ]  \frac{d}{dt} [Al]_e =s_e ⋅⋅⋅(8)$$

$$ [Al]_e =[[Al]_1  [Al]_2  [Al]_3  [Al]_4  [Al]_5  [Al]_6 ]^T⋅⋅⋅(9)$$
$$[M_e ]_{ij}=∫_{Ω_e}{(∇×\mathbf{w}_i )⋅(∇×\mathbf{w}_j )}dV=1/3v x_mn⋅1/3v x_{m'n'}×v=1/9v x_mn⋅x_{m'n'}⋅⋅⋅(10)$$
$$[N_e ]_{ij}=∫_{Ω_e}[\mathbf{w}_i⋅\mathbf{w}_j dV]=∫_{Ω_e}[(λ_k ∇λ_l-λ_l ∇λ_k )⋅(λ_{k'} ∇λ_{l'}-λ_{l'} ∇λ_{k'} )dV]$$
$$ =∫_{Ω_e}(λ_k λ_{k'} ∇λ_l⋅∇λ_{l'}-λ_k λ_{l'} ∇λ_l⋅∇λ_{k'}-λ_l λ_{k'} ∇λ_k⋅∇λ_{l'}+λ_l λ_{l'} ∇λ_k⋅∇λ_{k'} )dV$$
$$=∇λ_l⋅∇λ_{l'} ∫_{Ω_e}[λ_k λ_{k'} dV]-∇λ_l⋅∇λ_{k'} ∫_{Ω_e}[λ_k λ_{l'} dV]$$
$$-∇λ_k⋅∇λ_{l'} ∫_{Ω_e}[λ_l λ_{k'} dV]+∇λ_k⋅∇λ_{k'} ∫_{Ω_e}[λ_l λ_{l'} dV]⋅⋅⋅(11)$$
$$[s_e ]_i=μσ_e ∫_{Ω_e}[\mathbf{w}_i⋅(v×B_0 )dV]=μσ_e ∫_{Ω_e} [\mathbf{w}_i⋅(∑_{j=1}^4[λ_k [v×B_0]_j ])dV] =μσ_e ∑_{m=1}^4[(∫_{Ω_e}[\mathbf{w}_i λ_m dV])⋅[v×B_0 ]_m ] $$
$$=μσ_e ∑_{m=1}^4 [(∫_{Ω_e}[(λ_k ∇λ_l-λ_l ∇λ_k ) λ_m dV]) [⋅v×B_0 ]_m ] (12)$$
ここで、時間発展に後退オイラー法（2次の後退差分）を導入します。
$$\frac{\partial u_{n+2}}{\partial t}=\frac{1}{2Δt} (3u_{n+2}-4u_{n+1}+u_n ).⋅⋅⋅(13)$$
式 (13) を式 (8) に適用すると、次式を得ます。
$$ \left[ [M_e ]+\frac{3μσ_e}{2Δt} [N_e ] \right] [Al_e ]_{n+2}=\frac{1}{2Δt} [N_e ](4[Al_e]_{n+1}-[Al_e]_n )+[s_e]_{n+2}⋅⋅⋅(14)$$

# 付録B: 支配方程式のスケーリング
式 (14) より、各時間ステップで解く方程式は次のように書き直されます。

$$∑_{j=1}^6 \left[ ∫_{Ω_e}(\nabla \times \mathbf{w}_i)\cdot (\nabla \times \mathbf{w}_j )dV+\frac{3μσ_e}{2Δt} ∫_{Ω_e}\mathbf{w}_i⋅\mathbf{w}_j dV \right][Al_j]_{n+2} $$

$$=∑_{j=1}^6\frac{1}{2Δt} \left[ ∫_{Ω_e}\mathbf{w}_i \cdot \mathbf{w}_j dV\right](4[Al_j]_{n+1}-[Al_j]_n )+\mu \sigma_e ∑_{m=1}^4 \left[∫_{Ω_e}\mathbf{w}_i λ_m dV \right]⋅[\mathbf{v}×\mathbf{B}_0 ]_m⋅⋅⋅(15)$$

TMTGEM では、スケーリング長 L（TMTGEM では 1000 m に設定）を用いてスケーリングした式 (15) を解きます。$∇=\frac{1}{L} ∇'$,$w=\frac{1}{L} \mathbf{w}'$,$dV=L^3 dV'$, $Al[s⋅V/m⋅m]=∫A⋅dl=\frac{1}{L} Al'$（ここで $Al'$ の単位は [s⋅mV/km⋅km]、$E={\partial A}/{\partial t}$）を用いて式 (17) を書き直すことを考えます。さらに、$v=\frac{1}{L}⋅v'$ [mm/s]、$B_0=\frac{1}{L}^3⋅B_0'$ [nT] と仮定します。これらの置き換えは、w と ∇ が km 単位の座標系で w’ と ∇’ として構成されていることに対応します。

$$∑_{j=1}^6 \left[\int_{Ω_e}\frac{1}{L^2}\left( ∇'× \mathbf{w}'_i \right)\cdot \frac{1}{L^2} \left( ∇'×\mathbf{w}_j'\right) L^3 dV'+\frac{3μσ_e}{2Δt} ∫_{Ω_e}\frac{1}{L^2} \mathbf{w}_i'\cdot \mathbf{w}'_j L^3 dV' \right] \frac{1}{L} [Al'_j]_{n+2}$$

$$=∑_{j=1}^6 \frac{1}{2Δt} \left[∫_{Ω_e}\frac{1}{L^2} \mathbf{w}_i'\cdot \mathbf{w}'_j L^3 dV'\right]\frac{1}{L}\left(4[Al'_j]_{n+1}-[Al'_j]_n \right)$$

$$+μσ_e ∑_{m=1}^4 \left[\int_{Ω_e}\frac{1}{L} \mathbf{w}_i' λ_m L^3 dV'\right]\cdot \left[\left(\frac{1}{L} \mathbf{v}'\times\frac{1}{L^3}  \mathbf{B}'_0 \right)_m \right]_{n+2} $$

$$⟺∑_{j=1}^6\left[∫_{Ω_e}(∇'×\mathbf{w}_i' )⋅(∇'×\mathbf{w}_j')dV'+\frac{3μσ_e L^2}{2Δt}  ∫_{Ω_e}\mathbf{w}_i'⋅\mathbf{w}'_j dV'\right] [Al'_j]_{n+2}$$

$$=∑_{j=1}^6\frac{L^2}{2Δt} \left[∫_{Ω_e}\mathbf{w}_i'\cdot \mathbf{w}'_j dV'\right](4[Al'_j]_{n+1}-[Al'_j]_n )+ \mu \sigma_e ∑_{m=1}^4 \left[\int_{Ω_e}\mathbf{w}_i' λ_m dV'\right]⋅[[\mathbf{v}'×\mathbf{B}'_0]_m]_{n+2} $$


