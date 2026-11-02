"""一次性填好仓库身份信息：CITATION.cff、.zenodo.json、git 提交作者。

注意：本脚本把 creators 重写为**单一作者**。仓库现在是两位作者（Tang, Yuewu + Hu, Daijun），
CITATION.cff / .zenodo.json 为权威版本；除非确实要退回单作者，不要重新运行本脚本。

为什么需要它：Zenodo 建记录时读的是 CITATION.cff / .zenodo.json 里的 creators，
而 GitHub 页面上显示的是 git 提交里的 author。三处不一致，就会出现"DOI 记录写着
REPLACE-FAMILY-NAME，提交历史挂着 AVF Analysis"这种尴尬局面。

用法（先预览，确认无误再去掉 --dry-run）：

    python code/97_set_identity.py \\
        --family "Hu" --given "Xiaoming" \\
        --affiliation "Department of Nephrology, XXX Hospital" \\
        --email "huxm@example.edu.cn" \\
        --orcid 0000-0002-1825-0097 \\
        --fix-git-author --dry-run

ORCID 是可选的，但**不能填假号**——Zenodo 会做 ISO 7064 校验，`0000-0000-0000-0000`
这类占位会直接让归档失败（本仓库第一次失败就是这个原因）。本脚本会先替你校验。
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFF = ROOT / "CITATION.cff"
ZEN = ROOT / ".zenodo.json"


def orcid_checksum_ok(orcid: str) -> bool:
    """按 ISO 7064 MOD 11-2 校验 ORCID 末位校验码。"""
    body = re.sub(r"[^0-9X]", "", orcid.upper())
    if len(body) != 16:
        return False
    total = 0
    for ch in body[:15]:
        total = (total + int(ch)) * 2
    remainder = total % 11
    expected = (12 - remainder) % 11
    return ("X" if expected == 10 else str(expected)) == body[15]


def normalise_orcid(orcid: str) -> str:
    digits = re.sub(r"[^0-9X]", "", orcid.upper())
    return "-".join(digits[i:i + 4] for i in range(0, 16, 4))


def patch_cff(text: str, family: str, given: str, affil: str, orcid: str | None) -> str:
    text = re.sub(r'family-names:\s*"[^"]*"', f'family-names: "{family}"', text)
    text = re.sub(r'given-names:\s*"[^"]*"', f'given-names: "{given}"', text)
    text = re.sub(r'affiliation:\s*"[^"]*"', f'affiliation: "{affil}"', text)
    if orcid:
        line = f'    orcid: "https://orcid.org/{orcid}"'
        if re.search(r"^\s*#\s*orcid:", text, flags=re.M):
            text = re.sub(r"^\s*#\s*orcid:.*$", line, text, count=1, flags=re.M)
        elif "orcid:" not in text:
            text = text.replace(
                f'affiliation: "{affil}"', f'affiliation: "{affil}"\n{line}', 1
            )
    return text


def patch_zenodo(text: str, family: str, given: str, affil: str, orcid: str | None) -> str:
    data = json.loads(text)
    creator = {"name": f"{family}, {given}", "affiliation": affil}
    if orcid:
        creator["orcid"] = orcid
    data["creators"] = [creator]
    return json.dumps(data, ensure_ascii=False, indent=2) + "\n"


def run(cmd: list[str]) -> tuple[int, str]:
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True, help='姓，如 "Hu"')
    ap.add_argument("--given", required=True, help='名，如 "Xiaoming"（不要连写）')
    ap.add_argument("--affiliation", required=True)
    ap.add_argument("--email", default=None, help="git 提交邮箱；不填则只改姓名")
    ap.add_argument("--orcid", default=None, help="完整 ORCID iD；无则省略")
    ap.add_argument("--fix-git-author", action="store_true",
                    help="把已有提交的作者重设为当前身份（仓库已推送时会提示 force push）")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    orcid = None
    if a.orcid:
        if not orcid_checksum_ok(a.orcid):
            sys.exit(f"ORCID 校验位不通过：{a.orcid}\n"
                     "请到 https://orcid.org 复制完整 iD；宁可不填，也不要填占位。")
        orcid = normalise_orcid(a.orcid)
        print(f"ORCID 校验通过：{orcid}")

    print(f"\n即将写入：{a.given} {a.family} / {a.affiliation}"
          f"{' / ' + orcid if orcid else ' / (无 ORCID)'}")

    new_cff = patch_cff(CFF.read_text(encoding="utf-8"), a.family, a.given, a.affiliation, orcid)
    new_zen = patch_zenodo(ZEN.read_text(encoding="utf-8"), a.family, a.given, a.affiliation, orcid)

    # 写盘前先确认仍然是合法 YAML / JSON，避免把仓库改坏
    import yaml
    yaml.safe_load(new_cff)
    json.loads(new_zen)

    touched = []
    for path, content in ((CFF, new_cff), (ZEN, new_zen)):
        if path.read_text(encoding="utf-8") != content:
            touched.append(path.name)
            if not a.dry_run:
                path.write_text(content, encoding="utf-8")
    print(f"  文件：{'、'.join(touched) if touched else '无需改动'}")

    if a.email:
        ident = ["-c", f"user.name={a.given} {a.family}", "-c", f"user.email={a.email}"]
        for args in (["config", "user.name", f"{a.given} {a.family}"],
                     ["config", "user.email", a.email]):
            run(["git", *args])
        print(f"  git 身份：{a.given} {a.family} <{a.email}>")

        if a.fix_git_author:
            if a.dry_run:
                print("  [dry-run] 将执行 git rebase --root 重设全部提交作者")
            else:
                code, out = run(ident + ["rebase", "--root", "--exec",
                                         "git commit --amend --no-edit --reset-author"])
                if code != 0:
                    run(["rebase", "--abort"])
                    print(f"  [!] 重设作者失败，已回滚：{out.splitlines()[-1] if out else ''}")
                else:
                    code, out = run(["rev-list", "--count", "HEAD"])
                    print(f"  提交作者已重设为当前身份（共 {out} 个提交）")
                    print("  仓库已推送过，需要覆盖远端：git push --force-with-lease origin main")

    if a.dry_run:
        print("\n[dry-run] 未写入任何文件。去掉 --dry-run 正式执行。")
        return

    print("\n下一步：git add -A && git commit -m \"Set repository identity\" && git push")


if __name__ == "__main__":
    main()
