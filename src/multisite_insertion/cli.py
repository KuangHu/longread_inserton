"""Command-line interface for the assembly-based native-repeat route."""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from . import __version__
from .io import digest
from .pipeline import run
from .validation import validate

TOOLS = [
    "pals",
    "piler2",
    "muscle",
    "build_lmer_table",
    "RepeatScout",
    "blastn",
    "blastx",
    "makeblastdb",
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["run", "discover"]:
        p = sub.add_parser(
            name,
            help=(
                "Full analysis"
                if name == "run"
                else "Native outputs and coordinate inventory only"
            ),
        )
        p.add_argument("--manifest", type=Path, required=True)
        p.add_argument(
            "--out", type=Path, required=True, help="New or empty output directory"
        )
        p.add_argument(
            "--methods", nargs="+", choices=["PILER", "RS"], default=["PILER", "RS"]
        )
        p.add_argument(
            "--subset",
            choices=["tnp-supported", "all"],
            default="tnp-supported",
            help="Post-hoc boundary subset; discovery always retains every family",
        )
        p.add_argument(
            "--threads",
            type=int,
            default=2,
            help="BLASTX threads; native discovery defaults remain unchanged",
        )
        for tool in TOOLS:
            p.add_argument("--" + tool, default=tool)
    p = sub.add_parser(
        "validate",
        help="Verify source-coordinate/sequence consistency and native hashes",
    )
    p.add_argument("--run", type=Path, required=True)
    p = sub.add_parser(
        "doctor", help="Inspect installed tools and record executable hashes"
    )
    for tool in TOOLS:
        p.add_argument("--" + tool, default=tool)
    args = parser.parse_args(argv)
    try:
        if args.command == "doctor":
            report = {}
            for tool in TOOLS:
                path = shutil.which(getattr(args, tool))
                if not path:
                    report[tool] = {"available": False}
                    continue
                path = str(Path(path).resolve())
                flags = (
                    ["-version"]
                    if tool in ["blastn", "blastx", "makeblastdb", "muscle"]
                    else []
                )
                result = subprocess.run(
                    [path, *flags], capture_output=True, text=True, timeout=15
                )
                report[tool] = dict(
                    available=True,
                    path=path,
                    sha256=digest(path),
                    banner=(result.stdout + result.stderr)[:2000],
                    returncode=result.returncode,
                )
            print(json.dumps(report, indent=2))
            return 0 if all(r["available"] for r in report.values()) else 1
        if args.command == "validate":
            print(json.dumps(validate(args.run), indent=2))
            return 0
        if len(set(args.methods)) != len(args.methods):
            raise ValueError("Do not repeat a method in --methods")
        state = run(
            args.manifest,
            args.out,
            methods=args.methods,
            subset=args.subset,
            bins={t: getattr(args, t) for t in TOOLS},
            threads=args.threads,
            discovery_only=args.command == "discover",
        )
        print(json.dumps(state, indent=2))
        return 0
    except (
        ValueError,
        OSError,
        RuntimeError,
        KeyError,
        subprocess.TimeoutExpired,
    ) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
