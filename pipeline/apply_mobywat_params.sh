source "$SCRIPT_DIR/pipeline_common.sh"

PDB="${1:-system_ref.pdb}" #defaulting to mobywat convention

log "Adding REMARK section to reference PDB file."
sed -i '1i\
REMARK mobywat_reference_target [A]\
REMARK mobywat_reference_waters Auto\
#REMARK mobywat_min_ctol   1.0\
#REMARK mobywat_num_ctol   4\
#REMARK mobywat_step_ctol  0.50\
#REMARK mobywat_min_ptol   2.50\
#REMARK mobywat_num_ptol   1\
#REMARK mobywat_step_ptol  0\
#REMARK mobywat_min_stol   1.500\
#REMARK mobywat_num_stol   1\
#REMARK mobywat_step_stol  0
' "${PDB}"

log "Fixing TER with ID in the reference PDB file."
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