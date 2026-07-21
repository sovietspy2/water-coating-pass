#!/usr/bin/env python3
"""Extract the scalar results an engine reports, and show their spread across runs.

Structures tell you *that* two runs diverged; the scalars tell you how far the
optimiser or integrator actually travelled. TINKER-DETERMINISM.md section 2 measured
a 26 kcal/mol energy spread and an iteration count varying 888/881/761 between
identical `minimize` runs -- both invisible in an RMSD number, both decisive.

Reads one captured stdout/log per replicate and reports each scalar's min, max,
spread (max-min) and standard deviation.
"""

import argparse
import json
import math
import os
import re
import sys

# --------------------------------------------------------------------------
# Tinker
# --------------------------------------------------------------------------

# "  Final Function Value :          -3597.9345"
_TINKER_FINAL = re.compile(r"Final\s+(Function Value|RMS Gradient|Gradient Norm)\s*:\s*(-?[\d.]+(?:[eE][-+]?\d+)?)")
# Iteration table rows: "     1234     -3597.9345    0.0098   ..."
_TINKER_ITER_ROW = re.compile(r"^\s*(\d+)\s+-?\d+\.\d+\s+\d")
# "  Total Energy              -1234.5678 Kcal/mole"
_TINKER_MD_VALUE = re.compile(
    r"^\s*(Total Energy|Potential Energy|Kinetic Energy|Temperature)\s+(-?[\d.]+(?:[eE][-+]?\d+)?)"
)


def extract_tinker_minimize(text):
    values = {}
    for key, raw in _TINKER_FINAL.findall(text):
        values[key.lower().replace(" ", "_")] = float(raw)

    iterations = [int(m.group(1)) for m in (_TINKER_ITER_ROW.match(l) for l in text.splitlines()) if m]
    if iterations:
        values["iterations"] = float(max(iterations))
    return values


def extract_tinker_dynamic(text):
    """`dynamic` prints an averages block per save interval; the last one is the
    end state of the trajectory."""
    values = {}
    for line in text.splitlines():
        match = _TINKER_MD_VALUE.match(line)
        if match:
            values[match.group(1).lower().replace(" ", "_")] = float(match.group(2))
    return values


# --------------------------------------------------------------------------
# GROMACS
# --------------------------------------------------------------------------

# "Steepest Descents converged to Fmax < 200 in 1234 steps"
# "Polak-Ribiere Conjugate Gradients converged to Fmax < 750 in 5678 steps"
_GMX_CONVERGED = re.compile(r"converged to \S+ < \S+ in (\d+) steps")
# "Potential Energy  = -1.23456789e+04"
_GMX_SCALAR = re.compile(r"^\s*(Potential Energy|Maximum force|Norm of force)\s*=\s*(-?[\d.]+(?:[eE][-+]?\d+)?)")


def extract_gromacs_minimize(text):
    values = {}
    converged = _GMX_CONVERGED.findall(text)
    if converged:
        values["steps"] = float(converged[-1])
    for line in text.splitlines():
        match = _GMX_SCALAR.match(line)
        if match:
            values[match.group(1).lower().replace(" ", "_")] = float(match.group(2))
    return values


def extract_gromacs_md(text):
    """Parse the final `Energies (kJ/mol)` block from md.log.

    The block is fixed-width: name and value columns are 15 characters wide, with
    header rows and value rows alternating until a blank line.
    """
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.strip().startswith("Energies (kJ/mol)")]
    if not starts:
        return {}

    values = {}
    cursor = starts[-1] + 1
    while cursor + 1 < len(lines):
        header, row = lines[cursor], lines[cursor + 1]
        if not header.strip() or not row.strip():
            break
        names = [header[i:i + 15].strip() for i in range(0, len(header), 15)]
        numbers = [row[i:i + 15].strip() for i in range(0, len(row), 15)]
        for name, number in zip(names, numbers):
            if not name or not number:
                continue
            try:
                values[name.lower().replace(" ", "_").replace(".", "")] = float(number)
            except ValueError:
                pass
        cursor += 2
    return values


EXTRACTORS = {
    "tinker-minimize": extract_tinker_minimize,
    "tinker-dynamic": extract_tinker_dynamic,
    "gromacs-minimize": extract_gromacs_minimize,
    "gromacs-md": extract_gromacs_md,
}


# --------------------------------------------------------------------------
# Aggregation / reporting
# --------------------------------------------------------------------------


def aggregate(per_replicate):
    """Per-scalar min/max/spread/sd, keeping only scalars every replicate reported."""
    keys = set(per_replicate[0])
    for values in per_replicate[1:]:
        keys &= set(values)

    summary = {}
    for key in sorted(keys):
        series = [values[key] for values in per_replicate]
        mean = sum(series) / len(series)
        if len(set(series)) == 1:
            # Identical values must report exactly zero; the naive variance formula
            # leaves ~1e-12 rounding noise on large energies, which reads as signal.
            sd = 0.0
        elif len(series) > 1:
            variance = sum((v - mean) ** 2 for v in series) / (len(series) - 1)
            sd = math.sqrt(variance)
        else:
            sd = 0.0
        summary[key] = {
            "values": series,
            "min": min(series),
            "max": max(series),
            "spread": max(series) - min(series),
            "mean": mean,
            "sd": sd,
            "all_equal": len(set(series)) == 1,
        }
    return summary


def fmt(value):
    # 10 significant digits, not 6: Tinker energies are ~-2095.54 and the run-to-run
    # spread is ~0.002, so 6 digits renders min and max as the same string and hides
    # exactly the difference this report exists to show.
    if value == 0:
        return "0"
    if abs(value) < 1e-4 or abs(value) >= 1e9:
        return f"{value:.6e}"
    text = f"{value:.10g}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def fmt_stat(value):
    """Derived statistics (SD) do not need the full width raw values do -- printing a
    standard deviation to 10 significant digits implies precision that is not there."""
    if value == 0:
        return "0"
    if abs(value) < 1e-4 or abs(value) >= 1e9:
        return f"{value:.3e}"
    return f"{value:.6g}"


def write_markdown(summary, label, path, missing):
    out = [f"#### Engine scalars — {label}", ""]

    if not summary:
        out += ["_No scalars could be parsed from the engine logs._", ""]
    else:
        out += ["| scalar | all runs equal | min | max | spread (max-min) | SD |",
                "|---|---|---|---|---|---|"]
        for key, stat in summary.items():
            out.append(
                f"| `{key}` | {'yes' if stat['all_equal'] else '**no**'} "
                f"| {fmt(stat['min'])} | {fmt(stat['max'])} "
                f"| {fmt(stat['spread'])} | {fmt_stat(stat['sd'])} |"
            )
        out += ["", "Per-run values:", "",
                "| scalar | " + " | ".join(f"run {i + 1}" for i in range(len(next(iter(summary.values()))["values"]))) + " |",
                "|---" * (len(next(iter(summary.values()))["values"]) + 1) + "|"]
        for key, stat in summary.items():
            out.append(f"| `{key}` | " + " | ".join(fmt(v) for v in stat["values"]) + " |")
        out.append("")

    if missing:
        out += [f"> Some replicates reported scalars the others did not: {', '.join(sorted(missing))}", ""]

    with open(path, "a") as handle:
        handle.write("\n".join(out) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+", help="one captured engine log per replicate")
    parser.add_argument("--kind", required=True, choices=sorted(EXTRACTORS))
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-md")
    parser.add_argument("--out-json")
    return parser.parse_args()


def main():
    args = parse_args()

    missing_files = [p for p in args.logs if not os.path.isfile(p)]
    if missing_files:
        print("ERROR: missing log file(s): " + ", ".join(missing_files), file=sys.stderr)
        return 2

    extractor = EXTRACTORS[args.kind]
    per_replicate = []
    for path in args.logs:
        with open(path, errors="replace") as handle:
            per_replicate.append(extractor(handle.read()))

    all_keys = set().union(*per_replicate) if per_replicate else set()
    common_keys = set(per_replicate[0]).intersection(*per_replicate[1:]) if len(per_replicate) > 1 else all_keys
    missing = all_keys - common_keys

    summary = aggregate(per_replicate) if per_replicate else {}

    if args.out_json:
        with open(args.out_json, "w") as handle:
            json.dump({"label": args.label, "kind": args.kind, "scalars": summary}, handle, indent=2)

    if args.out_md:
        write_markdown(summary, args.label, args.out_md, missing)

    varying = [k for k, s in summary.items() if not s["all_equal"]]
    if varying:
        detail = ", ".join(f"{k} spread {fmt(summary[k]['spread'])}" for k in varying)
        print(f"scalars varying across runs: {detail}")
    elif summary:
        print(f"all {len(summary)} parsed scalar(s) identical across runs")
    else:
        print("no scalars parsed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
