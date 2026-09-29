#!/usr/bin/env bash
set -euo pipefail

PATH=/home/diriano/MERQURY.FK/:$PATH
export PATH
eval "$(mamba shell hook --shell bash)"
mamba activate genomescope2

species=(
    officinarum
    spontaneum
    robustum
    barberi
)

mkdir -p KatComp_results

for ((i = 0; i < ${#species[@]}; i++)); do
    SP1="${species[$i]}"
    KatComp -w10 -h10 -pdf -T10 "$SP1" SP803280.19.ktab "KatComp_results/KatComp__${SP1}__vs__SP803280"
    for ((j = i + 1; j < ${#species[@]}; j++)); do
        SP2="${species[$j]}"
        OUT="KatComp_results/KatComp__${SP1}__vs__${SP2}"

        # Confirm the expected FastK tables exist before starting.
        [[ -f "${SP1}.ktab" ]] || { echo "Missing: ${SP1}.ktab" >&2; exit 1; }
        [[ -f "${SP2}.ktab" ]] || { echo "Missing: ${SP2}.ktab" >&2; exit 1; }

        echo "Running: ${SP1} vs ${SP2}"
        KatComp -w10 -h10 -pdf -T10 "$SP1" "$SP2" "$OUT"
    done
done
