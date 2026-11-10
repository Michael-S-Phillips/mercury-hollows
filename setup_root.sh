#!/bin/bash
# Link this repository at ~/mnt/hollows-emissivity, the path under which the analyses ran (see README, "Running the scripts").
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"; LINK="$HOME/mnt/hollows-emissivity"
if [ -e "$LINK" ] || [ -L "$LINK" ]; then
  echo "$LINK already exists ($(readlink "$LINK" || echo 'not a link')); remove or rename it first."; exit 1
fi
mkdir -p "$HOME/mnt" && ln -s "$HERE" "$LINK" && echo "linked $LINK -> $HERE"
