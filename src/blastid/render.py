"""문제 행의 상위 hit 표를 PNG로 그린다."""
from __future__ import annotations

from pathlib import Path


def render_table(path: Path, title: str, hits: list[dict], n: int = 15) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    have = {f.name for f in font_manager.fontManager.ttflist}
    family = [
        f
        for f in (
            "Noto Sans CJK KR",
            "Noto Sans CJK SC",
            "WenQuanYi Micro Hei",
            "WenQuanYi Zen Hei",
            "NanumGothic",
            "Gothic A1",
            "DejaVu Sans",
        )
        if f in have
    ]
    plt.rcParams["font.family"] = family or ["DejaVu Sans"]
    cols = [
        "#", "Description", "Scientific Name", "Max\nScore", "Total\nScore",
        "Query\nCover", "E value", "Per.\nIdent", "Acc.\nLen", "Accession",
    ]
    rows = []
    for i, h in enumerate(hits[:n]):
        t = h.get("title") or ""
        t = t if len(t) <= 62 else t[:59] + "..."
        rows.append(
            [
                i + 1,
                t,
                (h.get("name") or "")[:34],
                h.get("max_s"),
                h.get("total_s"),
                h.get("qc_s"),
                h.get("e_s"),
                h.get("pi_s"),
                h.get("acc_len"),
                h.get("accession"),
            ]
        )
    if not rows:
        rows = [["", "No significant similarity found", "", "", "", "", "", "", "", ""]]
    fig = plt.figure(figsize=(18, 0.29 * (len(rows) + 2) + 0.5), dpi=90)
    ax = fig.add_axes((0.005, 0.0, 0.99, 0.93))
    ax.axis("off")
    fig.text(0.006, 0.985, title, ha="left", va="top", fontsize=12, fontweight="bold")
    tb = ax.table(
        cellText=rows,
        colLabels=cols,
        loc="upper left",
        cellLoc="left",
        colWidths=[0.025, 0.34, 0.16, 0.05, 0.05, 0.05, 0.06, 0.055, 0.045, 0.085],
    )
    tb.auto_set_font_size(False)
    tb.set_fontsize(9.5)
    tb.scale(1, 1.45)
    for (r, _c), cell in tb.get_celld().items():
        cell.set_edgecolor("#BFBFBF")
        if r == 0:
            cell.set_height(cell.get_height() * 1.9)
            cell.get_text().set_fontsize(9)
            cell.set_facecolor("#1F4E78")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        elif r == 1:
            cell.set_facecolor("#E2EFDA")
        elif r % 2 == 0:
            cell.set_facecolor("#F5F8FB")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
