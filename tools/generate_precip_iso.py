#!/usr/bin/env python3
"""
generate_precip_iso.py — SWAT+ precip.iso generator

TxtInOut ディレクトリを引数として渡すと、weather-sta.cli と weather-wgn.cli から
気象局名・座標を自動抽出し、Open-Meteo ERA5 の 2001-2020 気候値と経験式を使って
δ¹⁸O・δD を推定して precip.iso を生成します。

精度の目安: δ¹⁸O ±2-4 ‰ (RMSE 対 GNIP データ)
参考文献:
  - Clark & Fritz (1997) Environmental Isotopes in Hydrogeology
  - Bowen & Revenaugh (2003) JGR doi:10.1029/2003JD003560
  - Craig (1961) GMWL: δD = 8·δ¹⁸O + 10

使用例:
  # TxtInOut フォルダから気象局を自動検出
  python3 tools/generate_precip_iso.py /path/to/TxtInOut

  # 出力先を指定（デフォルトは TxtInOut/precip.iso）
  python3 tools/generate_precip_iso.py /path/to/TxtInOut --out /path/to/precip.iso

  # δ¹⁸O のみ
  python3 tools/generate_precip_iso.py /path/to/TxtInOut --num_iso 1
"""

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# ---------------------------------------------------------------------------
# weather-sta.cli / weather-wgn.cli パーサー
# ---------------------------------------------------------------------------

def parse_weather_sta(path: str) -> list[tuple[str, str]]:
    """
    weather-sta.cli を読み込み、[(sta_name, wgn_name), ...] を返す。
    sta_name: SWAT+ オブジェクト名（precip.iso の局ラベルに使用）
    wgn_name: weather-wgn.cli 内の WGN ステーション名
    """
    stations = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    # 行 1: ファイルヘッダ、行 2: 列ラベル → スキップ
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 2:
            stations.append((parts[0], parts[1]))
    return stations


def parse_weather_wgn(path: str) -> dict[str, tuple[float, float, float]]:
    """
    weather-wgn.cli を読み込み、{wgn_name: (lat, lon, elev_m)} を返す。

    フォーマット:
      行 1: ファイルヘッダ
      各局の先頭行: name lat lon elev num_years  (先頭が空白でない)
      続く行:  列ヘッダ + 12 行の月別値  (先頭が空白)
    """
    wgn = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0].isspace():
                continue                      # 列ヘッダ・月別値行はスキップ
            parts = line.split()
            if len(parts) < 4:
                continue                      # ファイルヘッダ行などはスキップ
            name = parts[0]
            if ":" in name:
                continue                      # "weather-wgn.cli:" ヘッダ行
            try:
                lat  = float(parts[1])
                lon  = float(parts[2])
                elev = float(parts[3])
                wgn[name] = (lat, lon, elev)
            except ValueError:
                pass
    return wgn


def get_stations_from_txtinout(txtinout: str) -> list[dict]:
    """
    TxtInOut フォルダから気象局リストを組み立てて返す。
    各要素: {"name": sta_name, "lat": float, "lon": float, "elev": float}
    """
    sta_path = os.path.join(txtinout, "weather-sta.cli")
    wgn_path = os.path.join(txtinout, "weather-wgn.cli")

    for p in (sta_path, wgn_path):
        if not os.path.isfile(p):
            sys.exit(f"[ERROR] 見つかりません: {p}")

    sta_list = parse_weather_sta(sta_path)
    wgn_dict = parse_weather_wgn(wgn_path)

    stations = []
    missing  = []
    for sta_name, wgn_name in sta_list:
        if wgn_name not in wgn_dict:
            missing.append(f"{sta_name} → wgn={wgn_name}")
            continue
        lat, lon, elev = wgn_dict[wgn_name]
        stations.append({"name": sta_name, "lat": lat, "lon": lon, "elev": elev})

    if missing:
        print(f"[WARN] weather-wgn.cli に対応する WGN が見つからない局:",
              file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)

    if not stations:
        sys.exit("[ERROR] 有効な気象局が 1 件も見つかりませんでした。")

    return stations


# ---------------------------------------------------------------------------
# Open-Meteo ERA5 archive API（日別 → 月次気候値）
# ---------------------------------------------------------------------------

def fetch_monthly_era5(lat: float, lon: float,
                       start_year: int = 2001, end_year: int = 2020) -> dict:
    """
    ERA5 の日別データを取得し、{elevation, T_monthly[12], P_monthly[12]} を返す。
    T_monthly: 月平均気温 [°C]
    P_monthly: 月合計降水量の年平均 [mm]
    """
    params = {
        "latitude":   lat,
        "longitude":  lon,
        "start_date": f"{start_year}-01-01",
        "end_date":   f"{end_year}-12-31",
        "daily":      "temperature_2m_mean,precipitation_sum",
        "timezone":   "UTC",
    }
    url = ("https://archive-api.open-meteo.com/v1/archive?"
           + urllib.parse.urlencode(params))

    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = json.load(resp)
    except urllib.error.HTTPError as e:
        sys.exit(f"[ERROR] Open-Meteo API HTTP {e.code}: {e.reason}\n  URL: {url}")
    except urllib.error.URLError as e:
        sys.exit(f"[ERROR] Open-Meteo に接続できません: {e.reason}")

    elev  = float(data.get("elevation", 0.0))
    dates = data["daily"]["time"]
    T_raw = data["daily"]["temperature_2m_mean"]
    P_raw = data["daily"]["precipitation_sum"]

    T_mo = _daily_to_monthly(dates, T_raw, agg="mean")
    P_mo = _daily_to_monthly(dates, P_raw, agg="sum")

    return {"elevation": elev, "T_monthly": T_mo, "P_monthly": P_mo}


def _daily_to_monthly(dates: list, values: list, agg: str = "mean") -> list:
    """
    日別値を 12 ヶ月の気候値に集計する。
    agg="mean": 月平均（気温用）
    agg="sum":  月合計の年数平均（降水量用）
    """
    sums        = [0.0] * 12
    day_counts  = [0]   * 12
    year_set    = [set() for _ in range(12)]

    for d, v in zip(dates, values):
        if v is None:
            continue
        mo   = int(d[5:7]) - 1
        year = int(d[:4])
        sums[mo]       += v
        day_counts[mo] += 1
        year_set[mo].add(year)

    if agg == "mean":
        return [sums[i] / day_counts[i] if day_counts[i] > 0 else float("nan")
                for i in range(12)]
    else:
        n_years = [len(year_set[i]) for i in range(12)]
        return [sums[i] / n_years[i] if n_years[i] > 0 else float("nan")
                for i in range(12)]


# ---------------------------------------------------------------------------
# δ¹⁸O 推定（経験式）
# ---------------------------------------------------------------------------

def estimate_d18O(T_monthly: list, lat: float, alt_m: float) -> list:
    """
    月別 δ¹⁸O (‰ vs VSMOW) を推定する。

    [年均値]
      δ¹⁸O_ann = 0.521·MAT − 0.006·|lat| − 0.002·alt_m − 12.0
      (Bowen & Revenaugh 2003 簡略版; 全球 RMSE ≈ 2-3 ‰)

    [月別スケール]
      δ¹⁸O_m = δ¹⁸O_ann + slope_T · (T_m − MAT)
      slope_T は緯度で補間（大陸性中緯度に最適化）:
        |lat| ≤ 30°:  0.20 ‰/°C  （熱帯〜亜熱帯）
        |lat| ≈ 50°:  0.28 ‰/°C  （温帯大陸性）
        |lat| ≥ 65°:  0.40 ‰/°C  （高緯度; 温度効果が卓越）
    """
    MAT = sum(T_monthly) / 12.0

    d18O_ann = (0.521 * MAT
                - 0.006 * abs(lat)
                - 0.002 * alt_m
                - 12.0)
    d18O_ann = max(-25.0, min(-2.0, d18O_ann))

    alat = abs(lat)
    if alat <= 30.0:
        slope_T = 0.20
    elif alat >= 65.0:
        slope_T = 0.40
    else:
        slope_T = 0.20 + (alat - 30.0) / (65.0 - 30.0) * (0.40 - 0.20)

    return [round(d18O_ann + slope_T * (T - MAT), 2) for T in T_monthly]


def d18O_to_dD(d18O_mo: list) -> list:
    """GMWL: δD = 8·δ¹⁸O + 10 (Craig 1961)"""
    return [round(8.0 * d + 10.0, 1) for d in d18O_mo]


# ---------------------------------------------------------------------------
# precip.iso ライター
# ---------------------------------------------------------------------------

def write_precip_iso(path: str, stations: list,
                     iso_on: int, num_iso: int,
                     iso_k: float, iso_x: float,
                     min_comp_rain: float, min_comp_gw: float):
    with open(path, "w") as f:
        names = "+".join(s["name"] for s in stations)
        f.write(f"precip.iso: {names}  "
                f"[iso_on num_iso iso_k iso_x min_comp_rain min_comp_gw]\n")
        f.write(f"{iso_on}  {num_iso}  {iso_k}  {iso_x}  "
                f"{min_comp_rain}  {min_comp_gw}\n")
        for s in stations:
            d18O_str = "  ".join(f"{v:7.2f}" for v in s["d18O"])
            f.write(f"{s['name']}  {d18O_str}\n")
            if num_iso >= 2:
                dD_str = "  ".join(f"{v:7.1f}" for v in s["dD"])
                f.write(f"{s['name']}  {dD_str}\n")


# ---------------------------------------------------------------------------
# 診断サマリー出力
# ---------------------------------------------------------------------------

def print_summary(s: dict, num_iso: int):
    print(f"\n  [{s['name']}]  lat={s['lat']:.3f}°  lon={s['lon']:.3f}°  "
          f"elev={s['elev']:.0f} m")
    print(f"    MAT={s['MAT']:.1f} °C   MAP={s['MAP']:.0f} mm/yr")
    hdr   = "    " + "".join(f"  {m:>6}" for m in MONTHS)
    row_T = "    T(°C)" + "".join(f"  {v:6.1f}" for v in s["T_monthly"])
    row_d = "    d18O " + "".join(f"  {v:6.2f}" for v in s["d18O"])
    print(hdr);  print(row_T);  print(row_d)
    if num_iso >= 2:
        row_D = "    dD   " + "".join(f"  {v:6.1f}" for v in s["dD"])
        print(row_D)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SWAT+ precip.iso generator (Open-Meteo ERA5 + 経験式)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "txtinout",
        help="SWAT+ TxtInOut フォルダのパス "
             "(weather-sta.cli / weather-wgn.cli から気象局を自動検出)",
    )
    p.add_argument(
        "--out", default=None,
        help="出力ファイルパス（デフォルト: <txtinout>/precip.iso）",
    )
    p.add_argument("--start_year", type=int, default=2001,
                   help="ERA5 集計開始年（デフォルト: 2001）")
    p.add_argument("--end_year",   type=int, default=2020,
                   help="ERA5 集計終了年（デフォルト: 2020）")
    p.add_argument("--iso_on",  type=int,   default=1)
    p.add_argument("--num_iso", type=int,   default=2, choices=[1, 2],
                   help="1=δ¹⁸O のみ, 2=δ¹⁸O+δD（デフォルト: 2）")
    p.add_argument("--iso_k",   type=float, default=1.0,
                   help="季節性係数（デフォルト: 1.0）")
    p.add_argument("--iso_x",   type=float, default=0.9,
                   help="交換係数（デフォルト: 0.9）")
    p.add_argument("--min_comp_rain", type=float, default=0.0)
    p.add_argument("--min_comp_gw",   type=float, default=0.0)
    return p


def main():
    args = build_parser().parse_args()

    txtinout = os.path.abspath(args.txtinout)
    if not os.path.isdir(txtinout):
        sys.exit(f"[ERROR] ディレクトリが見つかりません: {txtinout}")

    out_path = args.out or os.path.join(txtinout, "precip.iso")

    # --- 気象局リストを TxtInOut から抽出 ---
    stations = get_stations_from_txtinout(txtinout)

    print(f"気象局を {len(stations)} 件検出:", file=sys.stderr)
    for s in stations:
        print(f"  {s['name']:30s}  lat={s['lat']:.5f}  lon={s['lon']:.5f}  "
              f"elev={s['elev']:.1f} m", file=sys.stderr)

    print(f"\nOpen-Meteo ERA5 {args.start_year}-{args.end_year} 気候値を取得中...",
          file=sys.stderr)

    # --- 各局の ERA5 データ取得 → δ推定 ---
    result = []
    for i, s in enumerate(stations):
        if i > 0:
            time.sleep(0.5)      # API レート制限を避けるための小休止

        print(f"  ({i+1}/{len(stations)}) {s['name']} ...",
              file=sys.stderr, end="", flush=True)

        era5 = fetch_monthly_era5(s["lat"], s["lon"],
                                  args.start_year, args.end_year)
        print(" OK", file=sys.stderr)

        T_mo = era5["T_monthly"]
        P_mo = era5["P_monthly"]
        elev = era5["elevation"]    # ERA5 グリッドの標高（より正確）

        d18O_mo = estimate_d18O(T_mo, s["lat"], elev)
        dD_mo   = d18O_to_dD(d18O_mo)

        result.append({
            **s,
            "elev":      elev,
            "MAT":       sum(T_mo) / 12.0,
            "MAP":       sum(P_mo),
            "T_monthly": T_mo,
            "d18O":      d18O_mo,
            "dD":        dD_mo,
        })
        print_summary(result[-1], args.num_iso)

    # --- precip.iso 書き込み ---
    write_precip_iso(
        out_path, result,
        iso_on=args.iso_on, num_iso=args.num_iso,
        iso_k=args.iso_k,   iso_x=args.iso_x,
        min_comp_rain=args.min_comp_rain,
        min_comp_gw=args.min_comp_gw,
    )

    print(f"\n出力: {out_path}")
    print("注意: 精度の目安は RMSE ≈ 2-4 ‰ です。"
          "実測 GNIP データがある場合はそちらを優先してください。")


if __name__ == "__main__":
    main()
