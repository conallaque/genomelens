"""Every published link and asset resolves.

Unglamorous, and the first thing a reader hits if it is wrong.
"""
import re
from conftest import BENCH, ROOT, SAMPLES

REQUIRED = ("summary.json", "report-public.html", "report-public.pdf", "preview.png")


def test_every_sample_ships_the_full_artifact_set(sample):
    for name in REQUIRED:
        p = BENCH / sample / name
        assert p.is_file(), f"missing {sample}/{name}"
        assert p.stat().st_size > 0, f"empty {sample}/{name}"


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


def test_preview_images_are_real_images(sample):
    data = (BENCH / sample / "preview.png").read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "preview is not a PNG"
    assert len(data) > 20_000, "preview is too small to be readable"
