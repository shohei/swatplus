# SWAT+ 安定同位体トレーサー実装ガイド

## 概要

J2000水文モデルの同位体コンポーネント（Watson et al., 2022; Horita & Wesolowski, 1994）をSWAT+（Fortran）に移植した実装の解説です。降水・土壌水・地下水・河川の各コンパートメントを流れる水の**安定同位体比（δ¹⁸O、δD）**をトレースし、蒸発分別と水文過程ベースのハイドログラフ分離を計算します。

---

## 1. 背景：なぜ同位体トレーサーが必要か

水の安定同位体比（δ¹⁸O、δD）は降水イベントごとに固有の「指紋」を持ちます。土壌・地下水・河川に至るまでその指紋が保存されるため、

- どの流量成分（地表流出・中間流・基底流）が河川流量に何割寄与しているか（**ハイドログラフ分離**）
- 水が集水域内でどの経路をたどったか（**トランジットタイム推定**）
- 蒸発によって土壌水がどれだけ濃縮されるか（**蒸発分別**）

を定量化するための強力なツールになります。

---

## 2. J2000 の 6 コンポーネントと対応する Fortran ファイル

| Java クラス | Fortran ファイル | 役割 |
|---|---|---|
| `IsotopeMixer_` | `isotope_mixing.f90` : `iso_mix_binary` | 2 貯水池の体積加重混合 |
| `IsotopeMixer` | `isotope_mixing.f90` : `iso_mix_array` | 配列版双方向混合（混合比率付き） |
| `IsotopeMultiMixer` | `isotope_mixing.f90` : `iso_mix_multi` | 1 対多の体積比率分配混合 |
| `Isotope_fractionation` | `isotope_fractionation.f90` | 蒸発時の液相‐気相分別 |
| `Binary_isotope_mixing` | `isotope_hydsep.f90` | 2 成分ハイドログラフ分離 |
| `Tertiary_Isotope_mixing` | `isotope_hydsep.f90` | 3 成分比例分配 |

---

## 3. ファイル構成

```
src/
├── isotope_module.f90         # データ構造・グローバル変数（モジュール）
├── isotope_init.f90           # 配列確保・precip.iso 読み込み
├── isotope_mixing.f90         # 混合ユーティリティ（3 関数/サブルーチン）
├── isotope_rain.f90           # 降水δを土壌第 1 層に混合
├── isotope_lch.f90            # 土壌層・地下水へのδ値伝播
├── isotope_fractionation.f90  # 蒸発同位体分別
└── isotope_hydsep.f90         # ハイドログラフ分離（2 成分・3 成分）
```

変更ファイル：

```
src/hru_control.f90   # 毎タイムステップの HRU ループに呼び出しを追加
src/proc_hru.f90      # 初期化 (iso_init) を追加
```

入力ファイル（実行ディレクトリに配置）：

```
precip.iso            # 気象観測所ごとの月別降水δ値
```

---

## 4. データ構造（`isotope_module.f90`）

### 主な配列

```fortran
integer :: iso_on  = 0   ! 0=無効, 1=シミュレーション実行
integer :: num_iso = 1   ! 同位体の種類（1=δ¹⁸O のみ, 2=δ¹⁸O+δD）

real :: iso_soil(ihru, layer, iso)   ! 土壌各層の δ 値 (‰)
real :: iso_aqu(ihru, iso)           ! 浅層帯水層の δ 値 (‰)
real :: iso_precip(iwst, iso)        ! 気象観測所の日降水 δ 値 (‰)
real :: iso_precip_mo(iwst, 12, iso) ! 月別降水 δ 値（入力値）
```

### フラックス出力配列

```fortran
real :: iso_d_surq(ihru, iso)  ! 地表流出の δ 値 (‰)
real :: iso_d_latq(ihru, iso)  ! 中間流の δ 値 (‰)
real :: iso_d_tile(ihru, iso)  ! 暗渠排水の δ 値 (‰)
real :: iso_d_perc(ihru, iso)  ! 最下層浸透水の δ 値 (‰)
real :: iso_d_evap(ihru, iso)  ! 蒸発水蒸気の δ 値 (‰)
```

### ハイドログラフ分離結果

```fortran
real :: iso_comp_rain(ihru)    ! 降水寄与割合（2 成分）
real :: iso_comp_gw(ihru)      ! 地下水寄与割合（2 成分）
real :: iso_comp_rain_n(ihru)  ! 正規化降水寄与割合
real :: iso_comp_gw_n(ihru)    ! 正規化地下水寄与割合
real :: iso_comp_a(ihru)       ! 実測同位体流量（3 成分）
real :: iso_comp_b(ihru)       ! シミュレーション同位体流量（3 成分）
```

### グローバルパラメータ

| 変数 | デフォルト値 | 意味 |
|---|---|---|
| `iso_k` | 1.0 | 大気蒸気 δ_A 計算の季節係数 |
| `iso_x` | 0.9 | 土壌水の交換率（蒸発分別に使用） |
| `iso_min_comp_rain` | 0.0 | 2 成分分離の降水割合下限（キャリブレーション） |
| `iso_min_comp_gw` | 0.0 | 2 成分分離の地下水割合下限 |

---

## 5. 混合アルゴリズム（`isotope_mixing.f90`）

同位体の混合はすべて**体積加重平均（mass-flux mixing）**で行います。

### 5.1 二成分逐次混合（`iso_mix_binary`）

J2000 の `IsotopeMixer_` に相当します。

$$
\delta_{\text{混合}} = \frac{\delta_A \cdot V_A + \delta_B \cdot V_B}{V_A + V_B}
$$

```fortran
real function iso_mix_binary(delta_a, vol_a, delta_b, vol_b)
```

**用途：** 降水が土壌第 1 層に加わる場面、浸透水が下層に加わる場面など、水の合流が起きるあらゆる箇所で使用。

### 5.2 配列版双方向混合（`iso_mix_array`）

J2000 の `IsotopeMixer` に相当します。

$$
x_i = \frac{\delta_{A,i} \cdot V_{A,i} + \delta_{B,i} \cdot (V_{B,i} / p)}{V_{A,i} + V_{B,i} / p}
$$

ここで $p$ は混合比率（0 < p ≤ 1）。`bidir=.true.` のとき A と B の両方を $x_i$ に更新します。

```fortran
subroutine iso_mix_array(n, vol_a, conc_a, vol_b, conc_b, bidir, prop)
```

### 5.3 一対多混合（`iso_mix_multi`）

J2000 の `IsotopeMultiMixer` に相当します。1 つの供給源 A を複数の宛先 B(i) に体積比率で分配します。

$$
w_i = \frac{V_{B,i}}{\sum_j V_{B,j}}, \quad
x_i = \frac{\delta_A \cdot V_A w_i + \delta_{B,i} \cdot V_{B,i}}{V_A w_i + V_{B,i}}, \quad
\delta_A^{\text{new}} = \sum_i x_i \cdot w_i
$$

```fortran
subroutine iso_mix_multi(vol_a, conc_a, n, vol_b, conc_b)
```

---

## 6. 降水入力（`isotope_rain.f90`）

毎タイムステップ、気象観測所の月別テーブルから当日の降水 δ 値を取得し、`iso_mix_binary` で土壌第 1 層の δ 値を更新します。

```
δ_soil1(new) = iso_mix_binary(δ_soil1(old), V_soil1, δ_precip, V_precip)
```

- `V_precip` = `w%precip`（mm）
- `V_soil1` = `soil(j)%phys(1)%st`（mm）

---

## 7. 土壌層内の δ 値伝播（`isotope_lch.f90`）

`cs_lch.f90` の水フラックス計算と並行して動作します。各土壌層について**完全混合（well-mixed reservoir）仮定**を適用します。

```
各層の仮定：
  出力フラックスの δ = その層の現在の δ 値

層 jj の処理順：
  1. 上層からの浸透水を iso_mix_binary で混合
  2. 地表流出（第 1 層のみ）の δ = iso_d_surq
  3. 中間流の δ = iso_d_latq（全層の体積加重平均）
  4. 暗渠排水の δ = iso_d_tile
  5. 下向き浸透水の δ = iso_soil(layer, jj)（次層へ）
```

最下層からの浸透水は浅層帯水層に混合されます：

```
δ_aqu(new) = iso_mix_binary(δ_aqu(old), V_aqu_proxy, δ_perc, V_perc)
```

---

## 8. 蒸発同位体分別（`isotope_fractionation.f90`）

J2000 の `Isotope_fractionation` を移植。土壌表面からの蒸発（`es_day > 0` のとき）に対して、以下の手順で土壌第 1 層の δ 値を更新します。

### ステップ 1：平衡分別係数 α（Horita & Wesolowski, 1994）

δ¹⁸O に対する液相-気相平衡分別係数（*T* はケルビン）：

$$
\alpha^+ = \exp\!\left(\frac{1}{1000}\left[
\frac{1158.8 \, T^3}{10^9} - \frac{1620.1 \, T^2}{10^6} + \frac{794.84 \, T}{10^3} - 161.04 + \frac{2.9992 \times 10^9}{T^3}
\right]\right)
$$

$$
\varepsilon_{\text{mas}} = (\alpha^+ - 1) \times 1000 \quad [\text{‰}]
$$

### ステップ 2：動力学的分別係数 ε_k（Merlivat, 1978）

拡散率比 *D*(H₂¹⁸O)/*D*(H₂¹⁶O) ≈ 0.9755 を用いて：

$$
\varepsilon_k = 0.9755 \times (1 - 0.9755) \times 1000 \times (1 - h)
$$

ここで $h$ は相対湿度（0～1）。SWAT+ の `w%rhum` は分数（0～1）として格納されています。

### ステップ 3：濃縮勾配 *m* と大気蒸気 δ_A（Gibson et al., 2016）

$$
m = \frac{h - 10^{-3}(\varepsilon_k + \varepsilon_{\text{mas}}/\alpha^+)}{1 - h + 10^{-3}\varepsilon_k}
$$

$$
\delta_A = \frac{\delta_P - k \cdot \varepsilon_{\text{mas}}}{1 + \varepsilon_{\text{mas}} \times 10^{-3}}
$$

ここで $k$ は季節係数（デフォルト 1）、$\delta_P$ は当日の降水 δ 値。

### ステップ 4：定常濃縮極限値 δ*（Gonfiantini, 1986）

$$
\delta^* = \frac{h \cdot \delta_A + \varepsilon_k + \varepsilon_{\text{mas}}/\alpha^+}{h - 10^{-3}(\varepsilon_k + \varepsilon_{\text{mas}}/\alpha^+)}
$$

### ステップ 5：残存土壌水の δ（Craig & Gordon, 1965）

$$
\delta_S = \delta_{S_0} - \delta^*(1 - x)^m + \delta^*
$$

ここで $x$ は交換率パラメータ（デフォルト 0.9）。

### ステップ 6：蒸発水蒸気の δ_E

$$
\delta_E = \frac{(\delta_S - \varepsilon_{\text{mas}})/\alpha^+ - h \cdot \delta_A - \varepsilon_k}{1 - h + 10^{-3}\varepsilon_k}
$$

---

## 9. ハイドログラフ分離（`isotope_hydsep.f90`）

### 9.1 2 成分分離（Binary IHS）

降水と地下水を端成分として、Sklash & Farvolden (1979) の標準式を適用します。

$$
f_{\text{rain}} = \frac{\delta_{\text{stream}} - \delta_{\text{gw}}}{\delta_{\text{rain}} - \delta_{\text{gw}}}, \qquad
f_{\text{gw}}   = \frac{\delta_{\text{stream}} - \delta_{\text{rain}}}{\delta_{\text{gw}} - \delta_{\text{rain}}}
$$

> **注意：** J2000 の `Binary_isotope_mixing.java` には演算子優先順位のバグが含まれています。本実装では水文学の標準式を使用しています。

キャリブレーション下限（`iso_min_comp_rain`、`iso_min_comp_gw`）を適用後、中間流成分を除いた割合で正規化します。

### 9.2 3 成分比例分配（Tertiary）

$$
\text{comp\_A} = \delta_{\text{stream}} \times Q_{\text{total}}
$$

降水 δ データが利用可能（`iso_d_rain ≠ -99`）な場合：

$$
\text{comp\_B} = \delta_{\text{rain}} \times Q_{\text{surf}}
              + \delta_{\text{gw}}   \times Q_{\text{base}}
              + \delta_{\text{sw}}   \times Q_{\text{lat}}
$$

降水 δ データがない場合：

$$
\text{comp\_B} = \delta_{\text{gw}} \times Q_{\text{base}}
              + \delta_{\text{sw}}  \times Q_{\text{lat}}
$$

`comp_A ≈ comp_B` となれば、各成分の割合推定が観測された同位体比と整合しています。

> **注意：** J2000 の `Tertiary_Isotope_mixing.java` は `if(isotopeRain != -99)` の分岐ロジックが反転しています。本実装では科学的に正しい解釈を採用しています。

---

## 10. SWAT+ への組み込み方法

### 10.1 hru_control.f90 での呼び出しフロー

毎日・毎 HRU のシミュレーションループで `cs_lch` の直後に実行されます：

```fortran
! isotope tracking (stable water isotopes: delta-18O, delta-D)
if (iso_on == 1) then
  if (iso_atmo == "y") then
    call iso_rain    ! 降水δを土壌第1層に混合
  end if
  call iso_lch       ! 土壌層を通じたδ値伝播
  call iso_frac      ! 土壌第1層の蒸発分別
  call iso_hydsep    ! ハイドログラフ分離診断
end if
```

### 10.2 初期化（proc_hru.f90）

HRU 関連の初期化の最後に呼び出されます：

```fortran
call iso_init  ! stable water isotope initialization
```

---

## 11. 入力ファイル（`precip.iso`）

実行ディレクトリに配置します。

```
precip.iso  （1行目：タイトル）
iso_on  num_iso  k  x  min_comp_rain  min_comp_gw
station_name  d18O_jan  d18O_feb  ...  d18O_dec
[station_name  dD_jan  dD_feb  ...  dD_dec]
```

### 最小構成例（δ¹⁸O のみ、1 観測所）

```
precip.iso - Monthly precipitation delta-18O values
1  1  1.0  0.9  0.0  0.0
sta01  -5.5  -6.1  -7.2  -6.8  -5.0  -3.5  -3.2  -3.8  -4.5  -5.8  -6.5  -5.9
```

| 項目 | 値の例 | 説明 |
|---|---|---|
| `iso_on` | 1 | 同位体シミュレーションを有効化 |
| `num_iso` | 1 | δ¹⁸O のみ（2 にすると δD も追加） |
| `k` | 1.0 | 季節係数（通常は 1.0） |
| `x` | 0.9 | 土壌水交換率（0〜1） |
| `min_comp_rain` | 0.0 | 2 成分分離の降水割合下限 |
| `min_comp_gw` | 0.0 | 2 成分分離の地下水割合下限 |

ファイルが存在しない場合は `iso_on = 0`（無効）のままで動作に影響しません。

---

## 12. 主な制限事項と今後の拡張

| 項目 | 現状 | 拡張案 |
|---|---|---|
| 同位体種類 | δ¹⁸O のみ（`num_iso=1`）| `num_iso=2` でδD を追加、Majoube (1971) の係数を使用 |
| 帯水層モデル | HRU 単位の簡易代理変数 | `aqu_d%stor` を使った正確な体積に切り替え |
| 降水 δ 時系列 | 月平均値 | 日別時系列ファイルへの拡張 |
| 出力 | 配列変数として保持 | `hru_iso_output.f90` を追加して CSV 出力 |
| 蒸散分別 | 非分別と仮定（Evaristo et al., 2015） | 植物水利用の同位体分別を追加 |
| チャネルルーティング | 未実装 | `ch_water` の `iso` フィールドを追加 |

---

## 13. 参考文献

- Watson A, Vystavna Y, Kralisch S, Helmschrot J, van Rooyen J, Miller J (2022). Development of an isotope-enabled rainfall-runoff model. *Hydrological Processes*.
- Watson A, Birkel C, Kralisch S (2023). J2000-ISO isotope mixing components. FSU Jena / Stellenbosch University.
- Horita J, Wesolowski DJ (1994). Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature. *Geochimica et Cosmochimica Acta*, 58(16), 3425–3437.
- Gibson JJ et al. (2016). Stable isotope mass balance of lakes: a contemporary perspective. *Quaternary Science Reviews*, 131, 316–328.
- Merlivat L (1978). Molecular diffusivities of H₂¹⁶O, HD¹⁶O, and H₂¹⁸O in gases. *Journal of Chemical Physics*, 69(6), 2864–2871.
- Craig H, Gordon LI (1965). Deuterium and oxygen-18 variations in the ocean and marine atmosphere. *Proceedings of a Conference on Stable Isotopes in Oceanographic Studies*.
- Sklash MG, Farvolden RN (1979). The role of groundwater in storm runoff. *Journal of Hydrology*, 43, 45–65.
- Gonfiantini R (1986). Environmental isotopes in lake studies. *Handbook of Environmental Isotope Geochemistry*, 2, 113–168.
- Evaristo J, Jasechko S, McDonnell JJ (2015). Global separation of plant transpiration from groundwater and streamflow. *Nature*, 525, 91–94.
