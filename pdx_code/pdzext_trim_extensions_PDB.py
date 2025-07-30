from Bio import PDB
from Bio.Data.IUPACData import protein_letters_3to1
import freesasa
import re

def trim_structure(input_pdb, output_pdb, start_res, end_res, chain_id=None):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("struct", input_pdb)

    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            hetflag, resseq, icode = residue.id
            if hetflag != " ":  # Skip HETATM
                return False
            if chain_id and residue.get_parent().id != chain_id:
                return False
            # ✅ INCLUDE only residues within the user-defined range
            return start_res <= resseq <= end_res

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, ResidueSelect())


def calculate_rsa_with_biopython_names(pdb_file, output_tsv):
    # Parse structure with Biopython to map residue IDs to names
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("trimmed", pdb_file)
    res_lookup = {}  # (chain, resnum + icode) -> resname

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if hetflag != " ": continue
                key = (chain.id, f"{resseq}{icode.strip()}")
                res_lookup[key] = residue.get_resname()

    # Run FreeSASA
    sasa_structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(sasa_structure)
    areas = result.residueAreas()

    # Write output
    with open(output_tsv, 'w') as f:
        f.write("chain\tresnum\taa\trsa\n")
        for chain in areas:
            for res_id, res_info in areas[chain].items():
                rsa = res_info.relativeTotal
                key = (chain, str(res_id))  # match Biopython format
                aa_3 = res_lookup.get(key, "UNK")
                aa_1 = protein_letters_3to1.get(aa_3.capitalize(), "X")  # fallback to "X"
                f.write(f"{chain}\t{res_id}\t{aa_1}\t{rsa:.3f}\n")

# ---- USAGE ----
if __name__ == "__main__":
    input_pdb = "1be9.pdb"
    trimmed_pdb = "1be9_calculated_ddg_C.pdb"
    output_tsv = "1be9_calculated_ddg_C_rsa_output.tsv"

    # Define the region to keep (inclusive)
    start_res = 303
    end_res = 394
    chain_id = "A"

    trim_structure(input_pdb, trimmed_pdb, start_res, end_res, chain_id)
    calculate_rsa_with_biopython_names(trimmed_pdb, output_tsv)
    print(f"✅ Done. Output written to {output_tsv}")
