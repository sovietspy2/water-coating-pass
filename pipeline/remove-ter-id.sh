PDB="${1:-system_ref.pdb}" #defaulting to mobywat convention

echo "Fixing TER with ID in the reference PDB file. Also renumbering atoms if necessary."
awk '
BEGIN {n=1}
$1=="ATOM" || $1=="HETATM" {
    printf "%-6s%5d%s\n", substr($0,1,6), n, substr($0,12)
    n++
    next
}
$1=="TER" {
    print "TER"
    next
}
{print}
' ${PDB} > system_ref.pdb.tmp && mv system_ref.pdb.tmp ${PDB}