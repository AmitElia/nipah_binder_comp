from pathlib import Path
from Bio.PDB import MMCIFParser, PPBuilder

# --- config ---
INPUT_DIR = Path("./boltzgen_out_1000/protein_vs_2vsm_protein-anything/final_ranked_designs/final_10_designs")       # directory containing your .cif files
OUTPUT_FASTA = Path("./boltzgen_out_1000/protein_vs_2vsm_protein-anything/fall_complexes.fasta")
# If you ONLY want specific chains (e.g. binders), list them here; otherwise set to None
CHAIN_FILTER = {"A"}        # e.g. {"B"} or {"A", "B"} or None for all chains
# ---------------

parser = MMCIFParser(QUIET=True)
ppb = PPBuilder()

fasta_lines = []
seen_ids = set()

for cif_path in sorted(INPUT_DIR.glob("*.cif")):
    structure_id = cif_path.stem
    print(f"Parsing {cif_path}")

    structure = parser.get_structure(structure_id, str(cif_path))

    for model in structure:
        for chain in model:
            chain_id = chain.id

            # Optional: keep only specific chains
            if CHAIN_FILTER is not None and chain_id not in CHAIN_FILTER:
                continue

            # Build polypeptides from this chain
            peptides = ppb.build_peptides(chain)
            if not peptides:
                continue

            # Some chains are broken into multiple polypeptides; join them
            seq_str = "".join(str(pp.get_sequence()) for pp in peptides)
            if not seq_str:
                continue

            # Make a unique FASTA header
            header = f"{structure_id}_chain{chain_id}"
            # In case of duplicates, append an index
            base_header = header
            idx = 1
            while header in seen_ids:
                idx += 1
                header = f"{base_header}_{idx}"
            seen_ids.add(header)

            fasta_lines.append(f">{header}")
            fasta_lines.append(seq_str)

# Write output
if fasta_lines:
    OUTPUT_FASTA.write_text("\n".join(fasta_lines) + "\n")
    print(f"Wrote {len(fasta_lines)//2} sequences to {OUTPUT_FASTA}")
else:
    print("No sequences found in the provided .cif files.")