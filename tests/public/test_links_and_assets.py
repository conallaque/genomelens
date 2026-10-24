"""Every published link and asset resolves.

Unglamorous, and the first thing a reader hits if it is wrong.
"""
import re

from conftest import ROOT


def test_readme_links_and_images_resolve():
    text = (ROOT / "README.md").read_text()
    targets = re.findall(r"\[[^\]]+\]\(([^)]+)\)", text)
    targets += re.findall(r'(?:src|href)="([^"]+)"', text)
    broken = []
    for t in targets:
        if t.startswith(("http", "#", "mailto:")):
            continue
        path, _, frag = t.partition("#")
        if path and not (ROOT / path).exists():
            broken.append(t)
    assert not broken, f"broken README targets: {broken}"


def test_all_markdown_cross_references_resolve():
    broken = []
    for md in ROOT.rglob("*.md"):
        if ".git" in md.parts:
            continue
        for t in re.findall(r"\[[^\]]+\]\(([^)]+)\)", md.read_text()):
            if t.startswith(("http", "#", "mailto:")):
                continue
            path, _, _frag = t.partition("#")
            if path and not (md.parent / path).exists():
                broken.append(f"{md.relative_to(ROOT)} -> {t}")
    assert not broken, f"broken cross-references: {broken}"


def test_committed_pngs_are_real_images():
    pngs = sorted(p for p in ROOT.rglob("*.png") if ".git" not in p.parts)
    assert pngs, "no PNG found to check"
    for p in pngs:
        data = p.read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{p.relative_to(ROOT)} is not a PNG"
        assert len(data) > 10_000, f"{p.relative_to(ROOT)} is too small to be readable"
