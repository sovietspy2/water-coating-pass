#!/bin/sh

WHITELIST="pass5x5.sh pass5x1.sh mm_minim.sh cg.mdp md.mdp clean_pdb.py st.mdp cleanup.sh"
DRY_RUN=0

# Build the ( keep-this OR keep-that OR *.pdb ) expression as find args
set -- -name "*.pdb"
for f in $WHITELIST; do
  set -- "$@" -o -name "$f"
done

if [ "$DRY_RUN" -eq 1 ]; then
  find . -maxdepth 1 -type f ! \( "$@" \) -print
else
  find . -maxdepth 1 -type f ! \( "$@" \) -delete
fi
