"""98_relocate_paths.py

把 `code/` 下 46 个脚本里硬编码的工作根目录（默认 `D:/AVF_GSE_2611`）改写到新的
仓库路径，使代码在别人的机器上可直接运行。

为什么需要它
------------
流水线是在 Windows 上以固定路径开发的，脚本里大量使用绝对路径（约 70 处）。
与其逐个人工修改（易漏、易改坏已验证过的脚本），不如用一次机械改写，
并让改写本身可回溯、可逆转。

用法
----
    # 预览将要改动的行数与文件数（不写盘）
    python code/98_relocate_paths.py --to /home/me/avf --dry-run

    # 实际改写
    python code/98_relocate_paths.py --to /home/me/avf

    # 改回原路径（例如仍在本机原目录运行）
    python code/98_relocate_paths.py --to "D:/AVF_GSE_2611"

注意
----
* 只改 `code/*.py` 中的字面量，不动任何数据。
* 斜杠与反斜杠两种写法都会被替换。
* 改写后请跑 `python code/90_verify_numbers.py` 确认数字溯源仍全部通过。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

CODE = Path(__file__).resolve().parent
OLD_DEFAULT = "D:/AVF_GSE_2611"


def variants(root: str) -> list[str]:
    """返回同一根目录的常见字面量写法。"""
    r = root.rstrip("/\\")
    return [r, r.replace("/", "\\\\"), r.replace("/", "\\"), r.replace("\\", "/")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--to", required=True, help="新的仓库根目录，如 /home/me/avf")
    ap.add_argument("--from", dest="src", default=OLD_DEFAULT, help="当前脚本中的根目录")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    src_v = variants(a.src)
    dst = a.to.rstrip("/\\")

    n_files = n_hits = 0
    for p in sorted(CODE.glob("*.py")):
        if p.name == Path(__file__).name:
            continue
        s = p.read_text(encoding="utf-8")
        orig = s
        hits = 0
        for v in src_v:
            c = s.count(v)
            if c:
                hits += c
                s = s.replace(v, dst)
        if s != orig:
            n_files += 1
            n_hits += hits
            print(f"  {p.name:34s} {hits:3d} 处")
            if not a.dry_run:
                p.write_text(s, encoding="utf-8")

    verb = "将改写" if a.dry_run else "已改写"
    print(f"\n{verb} {n_files} 个文件、{n_hits} 处路径 → {dst}")
    if a.dry_run:
        print("（dry-run，未写盘；去掉 --dry-run 执行）")


if __name__ == "__main__":
    main()
