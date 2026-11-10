#!/bin/bash
# Downloads the MESSENGER MDIS 8-color MDR mosaic (PDS msgrmds_5001, 64 ppd, 55 tiles + labels), approved by Michael 3 Nov 2026.
B=https://planetarydata.jpl.nasa.gov/img/data/messenger
cd "$(dirname "$0")"
while read p h f; do
  for ext in LBL IMG; do
    g="${f%.IMG}.$ext"
    if [ -f "$g.done" ]; then continue; fi
    curl -sS --retry 5 --retry-delay 10 -C - -o "$g" "$B/$p/$h/$g" && touch "$g.done" || echo "FAILED $g" >> download_errors.log
  done
done < files.txt
echo "ALL DONE $(date)" >> download.log
