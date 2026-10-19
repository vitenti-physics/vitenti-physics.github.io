#!/usr/bin/env bash
# Measures, for every slide of a rendered reveal.js deck, how far its content
# extends past the right and bottom edges (in slide pixels), after MathJax has
# typeset. Serves _site on 127.0.0.1:8765 with server.js, which injects
# measure.js into <deck>-check.html and delays the page load event so headless
# Firefox screenshots only after the measurement has been posted back.
#   usage: mede_slides.sh <path under _site, e.g. courses/x/slides-aula1.html> [width height]
# Run from the site root after `quarto render ... --profile drafts`.
set -e
H=$(cd "$(dirname "$0")" && pwd)
W=${2:-1920}; HT=${3:-1080}
T=$(mktemp -d); mkdir -p "$T/ffprof"
node "$H/server.js" "$PWD/_site" 8765 "$T/report.json" & SRV=$!
trap 'kill $SRV' EXIT
sleep 1
timeout 90 firefox --headless --no-remote --profile "$T/ffprof" --window-size=$W,$HT \
  --screenshot "$T/shot.png" "http://127.0.0.1:8765/${1%.html}-check.html?measure" 2>&1 | grep -v headless || true
python -c "
import json
for r in json.load(open('$T/report.json')):
    print(f\"{r['i']:2d} right+{r['overflowRight']:4d} bottom+{r['overflowBottom']:4d} sec={r['secW']}x{r['secH']} {r['worst']:20s} {r['title']}\")"
