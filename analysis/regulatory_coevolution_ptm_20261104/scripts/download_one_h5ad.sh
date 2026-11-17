#!/bin/zsh
set -euo pipefail

if [[ $# -ne 3 ]]; then
  print -u2 "usage: $0 URL OUTPUT SHA256_OUTPUT"
  exit 2
fi

url="$1"
output="$2"
checksum_output="$3"
partial="${output}.partial"

mkdir -p "${output:h}" "${checksum_output:h}"
curl --fail --location --retry 3 --continue-at - --output "$partial" "$url"
h5dump -H "$partial" >/dev/null
mv "$partial" "$output"
shasum -a 256 "$output" >"${checksum_output}.partial"
mv "${checksum_output}.partial" "$checksum_output"
