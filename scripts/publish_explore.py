#!/usr/bin/env python3
"""Publishes a local build: the export (`data/explore/`) and the QA viewer (`logs/viewer/`) as
assets of a GitHub release of this repo, then starts the Pages workflow that deploys the viewer.

The build takes 30-40 min on a GitHub runner and ~5.5 min here, so it runs here and only its result
goes online. defileViz downloads `explore.zip` from the same release when it deploys (its
`scripts/fetch-explore.mjs`). Refused: an export marked as built from uncommitted code (here or in
defile-dataset), from another commit than HEAD or from a commit not pushed; a `--taxa` build since
the last full one; a viewer older than the export.

Usage:
    python scripts/publish_explore.py             # to the rolling `dev` pre-release
    python scripts/publish_explore.py --dry-run   # check and zip, upload nothing
    python scripts/publish_explore.py --tag v2025.1
"""

import argparse
import json
import os
import subprocess
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPORT_DIR = os.path.join(ROOT, "data", "explore")
VIEWER_DIR = os.path.join(ROOT, "logs", "viewer")
ZIP_DIR = os.path.join(ROOT, "logs", "publish")
DEV_TAG = "dev"  # rolling pre-release, until the first versioned release
EXPORT_ASSET = "explore.zip"  # the name defileViz downloads
VIEWER_ASSET = "viewer.zip"  # the name .github/workflows/pages.yml downloads
PAGES_WORKFLOW = "pages.yml"
UNCOMMITTED = "uncommitted"


def run(*args: str) -> str:
    return subprocess.run(
        args, capture_output=True, text=True, cwd=ROOT, check=True
    ).stdout.strip()


def check(export_dir: str, viewer_dir: str) -> dict:
    """The export's manifest, after the checks that make what goes online traceable to commits."""
    with open(os.path.join(export_dir, "manifest.json"), encoding="utf-8") as f:
        manifest = json.load(f)
    problems = []
    sha, dataset_sha = manifest.get("git_sha") or "", manifest.get("dataset_git_sha") or ""
    if UNCOMMITTED in sha:
        problems.append(f"built from uncommitted code ({sha}): commit, then rebuild")
    elif sha != run("git", "rev-parse", "--short", "HEAD"):
        problems.append(f"built at {sha}, HEAD is elsewhere: rebuild")
    elif not run("git", "branch", "-r", "--contains", "HEAD"):
        problems.append(f"{sha} is not pushed")
    if not dataset_sha or UNCOMMITTED in dataset_sha:
        problems.append(
            f"defile-dataset built from uncommitted code ({dataset_sha or 'unknown'}):"
            " commit there, rebuild it, then build here with --dataset"
        )
    catalogue = os.path.getmtime(os.path.join(export_dir, "taxa.json"))
    species = os.path.join(export_dir, "species")
    if any(os.path.getmtime(os.path.join(species, f)) > catalogue for f in os.listdir(species)):
        problems.append("species files rewritten since the last full build (--taxa): rebuild all")
    viewer = os.path.join(viewer_dir, "index.html")
    if not os.path.exists(viewer) or os.path.getmtime(viewer) < catalogue:
        problems.append("the viewer is older than the export: run scripts/explore_viewer.py")
    if problems:
        raise SystemExit("Not published:\n  " + "\n  ".join(problems))
    return manifest


def zip_folder(folder: str, path: str) -> str:
    """`folder`'s files at the zip's root."""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for base, _, files in os.walk(folder):
            for f in sorted(files):
                full = os.path.join(base, f)
                z.write(full, os.path.relpath(full, folder))
    print(f"  {os.path.relpath(path, ROOT)}: {os.path.getsize(path) / 1e6:.1f} MB")
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tag", default=DEV_TAG, help=f"release to publish to (default {DEV_TAG})")
    ap.add_argument("--dry-run", action="store_true", help="check and zip, upload nothing")
    args = ap.parse_args(argv)

    manifest = check(EXPORT_DIR, VIEWER_DIR)
    print(f"Export {manifest['git_sha']} of dataset {manifest['dataset_git_sha'][:7]}")
    os.makedirs(ZIP_DIR, exist_ok=True)
    assets = [
        zip_folder(EXPORT_DIR, os.path.join(ZIP_DIR, EXPORT_ASSET)),
        zip_folder(VIEWER_DIR, os.path.join(ZIP_DIR, VIEWER_ASSET)),
        os.path.join(EXPORT_DIR, "manifest.json"),
    ]
    if args.dry_run:
        print("Dry run: nothing uploaded")
        return 0

    if subprocess.run(
        ["gh", "release", "view", args.tag], capture_output=True, cwd=ROOT
    ).returncode:
        dev = args.tag == DEV_TAG
        run(
            "gh",
            "release",
            "create",
            args.tag,
            "--target",
            run("git", "rev-parse", "HEAD"),
            "--title",
            "Development build" if dev else args.tag,
            "--notes",
            "Local builds, replaced at each publish: manifest.json names the commits."
            if dev
            else f"Explore export {manifest['git_sha']}.",
            *(["--prerelease"] if dev else []),
        )
    run("gh", "release", "upload", args.tag, *assets, "--clobber")
    print(f"Uploaded to release {args.tag}")
    run("gh", "workflow", "run", PAGES_WORKFLOW, "-f", f"tag={args.tag}")
    print(f"Started {PAGES_WORKFLOW}: the viewer deploys in about a minute")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
