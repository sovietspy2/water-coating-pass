#!/bin/sh
# Reduce a multi-model PDB to its first model only, safely.
# - Creates: input.pdb.multimodel.orig
# - Rewrites: input.pdb
# - Leaves files with 0 or 1 MODEL records unchanged
# - Preserves header before first MODEL, except stale bookkeeping records
# - Preserves only first MODEL ... ENDMDL block
# - Preserves trailing non-model, non-coordinate records after all models
# - Writes exactly one final END record
# - Refuses malformed MODEL/ENDMDL structure

set -eu
umask 077
LC_ALL=C
export LC_ALL

usage() {
    printf 'Usage: %s <file.pdb>\n' "$0" >&2
    exit 1
}

fail() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

cleanup() {
    [ -n "${tmp-}" ] && [ -e "$tmp" ] && rm -f "$tmp"
}

make_tmp() {
    dir=$1
    base=$2

    if command -v mktemp >/dev/null 2>&1; then
        mktemp "$dir/.${base}.tmp.XXXXXX"
        return
    fi

    i=0
    while [ "$i" -lt 100 ]; do
        cand=$dir/.${base}.tmp.$$.${i}
        if ( set -C; : > "$cand" ) 2>/dev/null; then
            printf '%s\n' "$cand"
            return
        fi
        i=$((i + 1))
    done

    return 1
}

count_records() {
    file=$1
    name=$2
    awk -v want="$name" '
        {
            sub(/\r$/, "", $0)
            rec = substr($0, 1, 6)
            sub(/[[:space:]]+$/, "", rec)
            if (rec == want) c++
        }
        END { print c + 0 }
    ' "$file"
}

[ "$#" -eq 1 ] || usage

pdb=$1
[ -f "$pdb" ] || fail "File not found: $pdb"
[ -r "$pdb" ] || fail "File not readable: $pdb"
[ -w "$pdb" ] || fail "File not writable: $pdb"

model_count=$(count_records "$pdb" MODEL) || fail "Could not scan MODEL records"
endmdl_count=$(count_records "$pdb" ENDMDL) || fail "Could not scan ENDMDL records"

if [ "$model_count" -eq 0 ]; then
    printf 'No MODEL records found; file left unchanged: %s\n' "$pdb"
    exit 0
fi

if [ "$model_count" -eq 1 ]; then
    printf 'Only one MODEL record found; file left unchanged: %s\n' "$pdb"
    exit 0
fi

[ "$endmdl_count" -eq "$model_count" ] || fail \
    "Malformed PDB: MODEL/ENDMDL count mismatch ($model_count MODEL, $endmdl_count ENDMDL)"

backup="${pdb}.multimodel.orig"
[ ! -e "$backup" ] || fail "Backup already exists: $backup"

cp "$pdb" "$backup" || fail "Could not create backup: $backup"

dir=$(dirname "$pdb")
base=$(basename "$pdb")
tmp=$(make_tmp "$dir" "$base") || fail "Could not create temporary file"

trap cleanup EXIT HUP INT TERM

awk '
function recname(line, r) {
    r = substr(line, 1, 6)
    sub(/[[:space:]]+$/, "", r)
    return r
}

function is_coord_record(r) {
    return (r == "ATOM" || r == "HETATM" || r == "ANISOU" || r == "SIGATM" || r == "SIGUIJ")
}

function is_chain_term_record(r) {
    return (r == "TER")
}

function is_bookkeeping_record(r) {
    return (r == "MASTER" || r == "NUMMDL")
}

BEGIN {
    in_first_model = 0
    in_skipped_model = 0
    first_model_done = 0

    kept_models = 0
    kept_endmdls = 0
    kept_atoms = 0
}

{
    sub(/\r$/, "", $0)
    rec = recname($0)

    if (in_first_model) {
        if (rec == "MODEL") {
            print "Malformed PDB: nested MODEL inside first model" > "/dev/stderr"
            exit 2
        }

        if (rec == "ENDMDL") {
            print
            in_first_model = 0
            first_model_done = 1
            kept_endmdls++
            next
        }

        print
        if (rec == "ATOM" || rec == "HETATM") kept_atoms++
        next
    }

    if (in_skipped_model) {
        if (rec == "MODEL") {
            print "Malformed PDB: nested MODEL inside skipped model" > "/dev/stderr"
            exit 2
        }

        if (rec == "ENDMDL") {
            in_skipped_model = 0
            next
        }

        next
    }

    if (!first_model_done) {
        if (rec == "MODEL") {
            print
            in_first_model = 1
            kept_models++
            next
        }

        if (rec == "ENDMDL") {
            print "Malformed PDB: ENDMDL before first MODEL" > "/dev/stderr"
            exit 2
        }

        if (rec == "END") next
        if (is_bookkeeping_record(rec)) next

        print
        next
    }

    if (rec == "MODEL") {
        in_skipped_model = 1
        next
    }

    if (rec == "ENDMDL") {
        print "Malformed PDB: ENDMDL outside MODEL block" > "/dev/stderr"
        exit 2
    }

    if (rec == "END") next
    if (is_bookkeeping_record(rec)) next

    # After the retained first model, drop stray coordinate-section records that
    # would otherwise duplicate or corrupt the kept model/trailer structure.
    if (is_coord_record(rec) || is_chain_term_record(rec)) next

    print
}

END {
    if (in_first_model || in_skipped_model) {
        print "Malformed PDB: unterminated MODEL block" > "/dev/stderr"
        exit 2
    }

    if (kept_models != 1 || kept_endmdls != 1) {
        print "Internal error: output does not contain exactly one MODEL/ENDMDL pair" > "/dev/stderr"
        exit 2
    }

    if (kept_atoms == 0) {
        print "Malformed PDB: first model contains no ATOM/HETATM records" > "/dev/stderr"
        exit 2
    }

    printf "%-80s\n", "END"
}
' "$backup" > "$tmp" || fail "Failed while processing PDB"

# Final sanity check on the generated file.
awk '
function recname(line, r) {
    r = substr(line, 1, 6)
    sub(/[[:space:]]+$/, "", r)
    return r
}
function is_coordish(r) {
    return (r == "ATOM" || r == "HETATM" || r == "ANISOU" || r == "SIGATM" || r == "SIGUIJ" || r == "TER")
}
BEGIN {
    in_model = 0
    models = 0
    endmdls = 0
    ends = 0
    atoms = 0
    bad = 0
    last = ""
}
{
    sub(/\r$/, "", $0)
    rec = recname($0)
    last = rec

    if (rec == "MODEL") {
        models++
        if (in_model) bad = 1
        in_model = 1
        next
    }

    if (rec == "ENDMDL") {
        endmdls++
        if (!in_model) bad = 1
        in_model = 0
        next
    }

    if (rec == "END") {
        ends++
        next
    }

    if (in_model && (rec == "ATOM" || rec == "HETATM")) atoms++
    if (!in_model && is_coordish(rec)) bad = 1
}
END {
    if (in_model) bad = 1
    if (models != 1) bad = 1
    if (endmdls != 1) bad = 1
    if (ends != 1) bad = 1
    if (atoms == 0) bad = 1
    if (last != "END") bad = 1

    exit bad ? 1 : 0
}
' "$tmp" || fail "Generated PDB failed sanity checks"

mv "$tmp" "$pdb" || fail "Failed to replace original PDB"
trap - EXIT HUP INT TERM
cleanup

printf 'Backup created: %s\n' "$backup"
printf 'Updated file:   %s\n' "$pdb"