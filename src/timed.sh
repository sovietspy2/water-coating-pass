#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <command...> [runs]"
  echo "Example: $0 ./myprog input.dat 20"
  exit 1
fi

runs=10
last="${!#}"
if [[ "$last" =~ ^[0-9]+$ ]]; then
  runs="$last"
  set -- "${@:1:$(($#-1))}"
fi

cmd=( "$@" )

# Run once, keep and print stdout (before stats)
out="$("${cmd[@]}")"
printf "%s\n" "$out"

total_ms=0
min_ms=0
max_ms=0

for ((i=1; i<=runs; i++)); do
  start_ns=$(date +%s%N)
  "${cmd[@]}" > /dev/null
  end_ns=$(date +%s%N)

  dur_ms=$(((end_ns - start_ns) / 1000000))
  total_ms=$((total_ms + dur_ms))

  if (( i == 1 || dur_ms < min_ms )); then min_ms=$dur_ms; fi
  if (( i == 1 || dur_ms > max_ms )); then max_ms=$dur_ms; fi

  printf "run %d: %d ms\n" "$i" "$dur_ms"
done

avg_ms=$((total_ms / runs))

fmt_mmssms () {
  local ms=$1
  local m=$((ms / 60000))
  local s=$(((ms % 60000) / 1000))
  local r=$((ms % 1000))
  printf "%d:%02d.%03d" "$m" "$s" "$r"
}

printf "\nStats:\n"
printf "  total: %d ms (%s)\n" "$total_ms" "$(fmt_mmssms "$total_ms")"
printf "  avg  : %d ms (%s)\n" "$avg_ms"   "$(fmt_mmssms "$avg_ms")"
printf "  min  : %d ms (%s)\n" "$min_ms"   "$(fmt_mmssms "$min_ms")"
printf "  max  : %d ms (%s)\n" "$max_ms"   "$(fmt_mmssms "$max_ms")"
