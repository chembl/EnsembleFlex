import os
import sys
import argparse
import glob

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Extract ligands from PDB files and create reference structures.")
    parser.add_argument("-i", "--input", dest="input_path", required=True,
                        help="Input dataset directory path", metavar="PATH")
    parser.add_argument("-o", "--output", dest="output_dir", default="outdir",
                        help="Output directory name", metavar="STRING")
    parser.add_argument("-p", "--protein", dest="protein_pdb", default=None,
                        help="Optional PDB file to use as the protein structure", metavar="FILE")
    return parser.parse_args()

def validate_directories(input_path, output_path):
    """Validate input and output directories."""
    if not os.path.isabs(input_path):
        input_path = os.path.abspath(input_path)

    if not os.path.isdir(input_path):
        print("Error: The input directory does not exist.")
        sys.exit(1)

    if not os.path.isabs(output_path):
        output_path = os.path.abspath(os.path.join(os.getcwd(), output_path))

    if not os.path.isdir(output_path):
        try:
            os.makedirs(output_path)
        except OSError as e:
            print(f"Error: Could not create output directory '{output_path}'. {e}")
            sys.exit(1)

    return input_path, output_path

def validate_protein_pdb(protein_pdb, input_path):
    """Validate the protein PDB file."""
    if protein_pdb:
        if not os.path.isabs(protein_pdb):
            protein_pdb = os.path.abspath(protein_pdb)

        if not os.path.isfile(protein_pdb):
            print(f"Error: The specified protein PDB file '{protein_pdb}' does not exist.")
            sys.exit(1)
    else:
        pdb_files = sorted(glob.glob(os.path.join(input_path, "*.pdb")))
        if not pdb_files:
            print("Error: No PDB files found in the input directory.")
            sys.exit(1)
        protein_pdb = pdb_files[0]  # Default protein structure
        print(f"Using first PDB in input directory as protein structure: {protein_pdb}")

    return protein_pdb

def extract_hetatm_lines(pdb_file):
    """Extracts all HETATM lines from a PDB file."""
    hetatm_lines = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("HETATM"):
                hetatm_lines.append(line)
    return hetatm_lines

def extract_protein_lines(pdb_file):
    """Extracts protein structure lines (ATOM, TER, HEADER, TITLE, REMARK) from a PDB file."""
    protein_lines = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "TER", "HEADER", "TITLE", "REMARK")):
                protein_lines.append(line)
    return protein_lines

def create_pdb_outputs(input_path, output_path, protein_pdb):
    """Creates two PDB files: one with protein and ligands, one with ligands only (as multi-models)."""
    pdb_files = sorted(glob.glob(os.path.join(input_path, "*.pdb")))

    if not pdb_files:
        print("Error: No PDB files found in the input directory.")
        sys.exit(1)

    # Extract protein structure from specified or default protein PDB
    protein_lines = extract_protein_lines(protein_pdb)

    # Define output file paths
    protein_with_ligands_pdb = os.path.join(output_path, "ligands_with_protein.pdb")
    ligands_only_pdb = os.path.join(output_path, "ligands_only.pdb")

    # Write protein + ligands in multi-model format
    with open(protein_with_ligands_pdb, 'w') as out_file:
        out_file.write(f"MODEL 1\n")
        out_file.writelines(protein_lines)  # Write protein structure
        out_file.write("ENDMDL\n")

        model_number = 2
        for pdb_file in pdb_files:
            ligand_lines = extract_hetatm_lines(pdb_file)
            if ligand_lines:
                out_file.write(f"MODEL {model_number}\n")
                out_file.writelines(ligand_lines)
                out_file.write("ENDMDL\n")
                model_number += 1

        out_file.write("END\n")  # Ensure proper termination

    # Write ligands-only PDB in multi-model format
    with open(ligands_only_pdb, 'w') as out_file:
        model_number = 1
        for pdb_file in pdb_files:
            ligand_lines = extract_hetatm_lines(pdb_file)
            if ligand_lines:
                out_file.write(f"MODEL {model_number}\n")
                out_file.writelines(ligand_lines)
                out_file.write("ENDMDL\n")
                model_number += 1

        out_file.write("END\n")  # Ensure proper termination

    print(f"Protein PDB with ligands (multi-model) created: {protein_with_ligands_pdb}")
    print(f"Ligands-only PDB (multi-model) created: {ligands_only_pdb}")


def main():
    """Main function to handle argument parsing and processing."""
    # Parse command-line arguments
    args = parse_arguments()

    # Validate input and output directories
    input_path, output_path = validate_directories(args.input_path, args.output_dir)

    # Validate protein PDB file or use the first PDB in the input directory
    protein_pdb = validate_protein_pdb(args.protein_pdb, input_path)

    # Change working directory to output_path
    os.chdir(output_path)

    # Create PDB outputs
    create_pdb_outputs(input_path, output_path, protein_pdb)

if __name__ == "__main__":
    main()
