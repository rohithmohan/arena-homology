#!/usr/bin/env bash
# End-to-end: fetch bundle -> extract -> RCSB sequences -> MMseqs2 clustering -> re-scoring -> report.
# Requires: node, python3 with numpy, scipy and anarci (pip install -r requirements.txt), and HMMER's
# hmmscan on PATH (used by ANARCI; e.g. `apt-get install hmmer` or `conda install -c bioconda hmmer`).
# MMseqs2 is taken from $MMSEQS, then PATH, else the release 18 static binary is downloaded into tools/.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/raw data/derived

if [[ -z "${MMSEQS:-}" ]] && ! command -v mmseqs >/dev/null && [[ ! -x tools/mmseqs/bin/mmseqs ]]; then
  echo "downloading MMseqs2 release 18 (linux avx2)"
  mkdir -p tools
  curl -sSL https://github.com/soedinglab/MMseqs2/releases/download/18-8cc5c/mmseqs-linux-avx2.tar.gz | tar xz -C tools
fi

curl -sS https://cofoldarena.ai/ -o data/raw/index.html
BUNDLE=$(grep -o 'assets/index-[^"]*\.js' data/raw/index.html | head -1)
echo "bundle: $BUNDLE"
curl -sS "https://cofoldarena.ai/$BUNDLE" -o data/raw/bundle.js
node scripts/01_extract.js data/raw/bundle.js data/derived
python3 -I scripts/02_fetch_sequences.py      # caches RCSB responses in data/raw/rcsb/
python3 -I scripts/03_cluster.py
python3 -I scripts/04_rescore.py              # asserts the published leaderboard numbers first
python3 -I scripts/05_report.py
