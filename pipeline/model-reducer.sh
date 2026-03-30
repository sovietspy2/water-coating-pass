#!/bin/sh
# sh script to reduce a multi-model PDB to its first model only.
# - Creates: input.pdb.multimodel.orig
# - Rewrites: input.pdb
# - Preserves header before first MODEL
# - Preserves only first model coordinates
# - Preserves trailing non-model metadata after the model section
# - Drops subsequent MODEL/ENDMDL blocks and their coordinate contents

set -eu

usage() {
    printf 'Usage: %s <file.pdb>\n' "$0" >&2
    exit 1
}

fail() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

[ "$#" -eq 1 ] || usage

pdb=$1
[ -f "$pdb" ] || fail "File not found: $pdb"
[ -r "$pdb" ] || fail "File not readable: $pdb"
[ -w "$pdb" ] || fail "File not writable: $pdb"

backup="${pdb}.multimodel.orig"
tmp="${pdb}.tmp.$$"

[ ! -e "$backup" ] || fail "Backup already exists: $backup"

cp -- "$pdb" "$backup" 2>/dev/null || cp "$pdb" "$backup" || fail "Could not create backup: $backup"

cleanup() {
    rm -f -- "$tmp" 2>/dev/null || rm -f "$tmp"
}
trap cleanup EXIT HUP INT TERM

awk '
BEGIN {
    in_model = 0
    first_model_seen = 0
    keep_model = 0
    after_models = 0
    saw_end = 0
}

{
    rec = substr($0, 1, 6)

    if (rec == "MODEL ") {
        if (first_model_seen == 0) {
            first_model_seen = 1
            in_model = 1
            keep_model = 1
            after_models = 0
            print
        } else {
            in_model = 1
            keep_model = 0
            after_models = 0
        }
        next
    }

    if (rec == "ENDMDL") {
        if (in_model && keep_model) {
            print
        }
        in_model = 0
        keep_model = 0
        after_models = first_model_seen ? 1 : 0
        next
    }

    if (in_model) {
        if (keep_model) {
            print
        }
        next
    }

    if (first_model_seen == 0) {
        if (rec == "END   ") {
            saw_end = 1
            next
        }
        print
        next
    }

    if (after_models) {
        if (rec == "END   ") {
            saw_end = 1
            next
        }
        print
        next
    }

    print
}

END {
    if (saw_end == 0) {
        print "END"
    } else {
        print "END"
    }
}
' "$backup" > "$tmp" || fail "Failed while processing PDB"

[ -s "$tmp" ] || fail "Processing produced empty output"

mv -- "$tmp" "$pdb" 2>/dev/null || mv "$tmp" "$pdb" || fail "Failed to replace original PDB"

trap - EXIT HUP INT TERM
cleanup

printf 'Backup created: %s\n' "$backup"
printf 'Updated file:   %s\n' "$pdb"