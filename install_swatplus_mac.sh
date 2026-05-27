#!/usr/bin/env bash
# install_swatplus_mac.sh
#
# SWAT+ Editor に自コンパイルした ARM64 バイナリを組み込む。
#
# 使い方:
#   ./install_swatplus_mac.sh [オプション] [バイナリパス]
#
# オプション:
#   --restore   バックアップから元の x86_64 バイナリに戻す
#
# 引数省略時:
#   バイナリパス: スクリプトと同じディレクトリで最新の arm64-Rel バイナリを自動検索
#
# 例:
#   ./install_swatplus_mac.sh
#   ./install_swatplus_mac.sh ./build/release/swatplus-61.0.2.61-...-gnu-mac_arm64-Rel
#   ./install_swatplus_mac.sh --restore

set -euo pipefail

# ---------------------------------------------------------------------------
# 定数
# ---------------------------------------------------------------------------
EDITOR_NAME="SWATPlusEditor.app"
SWAT_EXE_SUBDIR="Contents/Resources/app.asar.unpacked/static/swat_exe"
TARGET_FILENAME="rev60.5.7_64rel_mac"
BACKUP_SUFFIX=".orig_x86_64"

# ---------------------------------------------------------------------------
# 引数解析
# ---------------------------------------------------------------------------
DO_RESTORE=0
ARM64_BIN=""

for arg in "$@"; do
    case "$arg" in
        --restore) DO_RESTORE=1 ;;
        *)         ARM64_BIN="$arg" ;;
    esac
done

# ---------------------------------------------------------------------------
# SWAT+ Editor の検索
# ---------------------------------------------------------------------------
EDITOR_APP=""
for candidate in \
    "$HOME/SWATPlus/SWATPlusEditor/$EDITOR_NAME" \
    "/Applications/$EDITOR_NAME" \
    "$HOME/Applications/$EDITOR_NAME"
do
    if [[ -d "$candidate" ]]; then
        EDITOR_APP="$candidate"
        break
    fi
done

if [[ -z "$EDITOR_APP" ]]; then
    echo "ERROR: SWAT+ Editor が見つかりません。" >&2
    echo "  検索先: ~/SWATPlus/SWATPlusEditor/, /Applications/, ~/Applications/" >&2
    exit 1
fi

SWAT_EXE_DIR="$EDITOR_APP/$SWAT_EXE_SUBDIR"
TARGET="$SWAT_EXE_DIR/$TARGET_FILENAME"
BACKUP="$TARGET$BACKUP_SUFFIX"

echo ">>> SWAT+ Editor: $EDITOR_APP"
echo ">>> 対象バイナリ:  $TARGET"
echo ""

# ---------------------------------------------------------------------------
# --restore モード
# ---------------------------------------------------------------------------
if [[ $DO_RESTORE -eq 1 ]]; then
    if [[ ! -f "$BACKUP" ]]; then
        echo "ERROR: バックアップが見つかりません: $BACKUP" >&2
        exit 1
    fi
    cp "$BACKUP" "$TARGET"
    echo "=== 復元完了 ==="
    file "$TARGET"
    exit 0
fi

# ---------------------------------------------------------------------------
# ARM64 バイナリの検索（引数省略時）
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "$ARM64_BIN" ]]; then
    # スクリプトと同じリポジトリの build/release/ から arm64-Rel を探す
    ARM64_BIN=$(find "$SCRIPT_DIR" -path "*/build/release/*arm64*Rel" \
                    ! -name "*.orig_*" 2>/dev/null \
                | sort -V | tail -1)
    if [[ -z "$ARM64_BIN" ]]; then
        echo "ERROR: ARM64 バイナリが見つかりません。" >&2
        echo "  先に cmake --build build/release を実行するか、パスを引数で指定してください。" >&2
        echo "  使い方: $0 <バイナリパス>" >&2
        exit 1
    fi
    echo ">>> 自動検出: $ARM64_BIN"
fi

if [[ ! -x "$ARM64_BIN" ]]; then
    echo "ERROR: バイナリが見つからないか実行権限がありません: $ARM64_BIN" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# アーキテクチャ確認
# ---------------------------------------------------------------------------
BIN_ARCH=$(file "$ARM64_BIN" | grep -o "arm64\|x86_64\|universal" | head -1)
if [[ "$BIN_ARCH" != "arm64" && "$BIN_ARCH" != "universal" ]]; then
    echo "ERROR: 指定バイナリが arm64 ではありません (検出: ${BIN_ARCH:-不明})。" >&2
    echo "  arm64 または universal バイナリを指定してください。" >&2
    exit 1
fi
echo ">>> アーキテクチャ: $BIN_ARCH  ✓"

# ---------------------------------------------------------------------------
# SWAT+ バージョン確認
# ---------------------------------------------------------------------------
SWAT_VER=$(cd /tmp && "$ARM64_BIN" 2>&1 | grep -i "Revision" | head -1 | tr -d ' ' || true)
echo ">>> バージョン: ${SWAT_VER:-不明}"
echo ""

# ---------------------------------------------------------------------------
# バックアップ
# ---------------------------------------------------------------------------
if [[ -f "$TARGET" ]]; then
    EXISTING_ARCH=$(file "$TARGET" | grep -o "arm64\|x86_64\|universal" | head -1)
    if [[ "$EXISTING_ARCH" == "arm64" && -f "$BACKUP" ]]; then
        echo ">>> 既に ARM64 バイナリがインストールされています（バックアップ済み）。"
        echo "    上書きしてよい場合は続行してください。"
    fi
    if [[ ! -f "$BACKUP" ]]; then
        cp "$TARGET" "$BACKUP"
        echo ">>> バックアップ: $BACKUP"
    else
        echo ">>> バックアップ既存のためスキップ: $BACKUP"
    fi
fi

# ---------------------------------------------------------------------------
# インストール
# ---------------------------------------------------------------------------
cp "$ARM64_BIN" "$TARGET"
chmod +x "$TARGET"

# macOS Quarantine 属性を除去（ダウンロードしたファイルの場合）
xattr -d com.apple.quarantine "$TARGET" 2>/dev/null || true

# ---------------------------------------------------------------------------
# 確認
# ---------------------------------------------------------------------------
echo ""
echo "=== インストール完了 ==="
file "$TARGET"
ls -lh "$TARGET"
echo ""
echo "元に戻すには:"
echo "  $0 --restore"
