"""Make a `conda env export` file safe to publish.

Replaces the environment name (an absolute cluster path) with a plain name, drops the prefix line,
and replaces the channel list with the public channels used to build the environment. Cluster
channel entries can contain access tokens (for example https://conda.anaconda.org/t/<token>/...)
and must not appear in a public repository.

Usage:
  python sanitize_env_export.py IN.yml OUT.yml [--name bvas]
"""
import argparse
import re

PUBLIC_CHANNELS = ["conda-forge", "bioconda"]


def sanitize(lines, name):
    out, in_channels = [], False
    for line in lines:
        if line.startswith("name:"):
            out.append(f"name: {name}\n")
            continue
        if line.startswith("prefix:"):
            continue
        if line.startswith("channels:"):
            out.append("channels:\n")
            out.extend(f"  - {c}\n" for c in PUBLIC_CHANNELS)
            in_channels = True
            continue
        if in_channels and line.startswith("  - "):
            continue
        in_channels = False
        out.append(line)
    text = "".join(out)
    if re.search(r"/t/[^/\s]+/", text):
        raise SystemExit("token-like URL still present after sanitizing")
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--name", default="bvas")
    a = ap.parse_args()
    with open(a.infile, encoding="utf-8") as fh:
        text = sanitize(fh.readlines(), a.name)
    with open(a.outfile, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


if __name__ == "__main__":
    main()
