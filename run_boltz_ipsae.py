#!/usr/bin/env python
"""
Batch Boltz-2 co-folding for binder–target complexes.

Given:
  - a FASTA with binder sequences (one or more records)
  - a target structure (e.g. 2VSM.cif) + chain id (e.g. A)

This script will:
  1. Extract the target amino-acid sequence from the CIF.
  2. Generate one Boltz-2 YAML per binder:
       - protein A = target (from 2VSM.cif chain A)
       - protein B = binder (from FASTA)
       - templates: use the CIF as a template for chain A
  3. Run `boltz predict` for each YAML.
  4. Optionally run ipSAE on the resulting PAE/CIF pairs.

YAML format is aligned with the official Boltz-2 schema:
  version: 1
  sequences:
    - protein:
        id: A
        sequence: <target_sequence>
    - protein:
        id: B
        sequence: <binder_sequence>
  templates:
    - cif: /abs/path/to/2VSM.cif
      chain_id: A

"""

import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Optional

import yaml

try:
    # Biopython is used only for extracting the sequence from the CIF
    from Bio.PDB import MMCIFParser, PPBuilder
except ImportError:
    MMCIFParser = None
    PPBuilder = None


# ---------------------------
# Small utilities
# ---------------------------

def log(msg: str, stream=sys.stdout) -> None:
    print(msg, file=stream)
    stream.flush()


def run_cmd(cmd: List[str], cwd: Path, log_path: Path) -> int:
    """
    Run a subprocess, tee stdout/stderr to a log file, and return exit code.
    """
    cwd = cwd.resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w") as f:
        f.write(f"$ {' '.join(cmd)}\n\n")
        f.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            f.write(line)
            f.flush()
        proc.wait()
        return proc.returncode


def short_hash(seq: str, n: int = 8) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:n]


def read_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    """
    Very simple FASTA reader -> yields (header, sequence).
    """
    header = None
    chunks: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and chunks:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None and chunks:
        yield header, "".join(chunks)


def infer_temperature_from_header(header: str) -> str:
    """
    Try to pull a temperature like 'T_0.1' or 'T=0.1' out of the header.
    Used only for naming; not required by Boltz.
    """
    m = re.search(r"[Tt][=_](\d+(\.\d+)?)", header)
    return m.group(1) if m else "NA"


# ---------------------------
# Target sequence from CIF
# ---------------------------

def get_target_sequence_from_cif(cif_path: Path, chain_id: str) -> str:
    """
    Extract a single-letter amino-acid sequence for a given chain from a mmCIF file.

    Requires Biopython:
        pip install biopython
    """
    if MMCIFParser is None or PPBuilder is None:
        raise RuntimeError(
            "Biopython is required for get_target_sequence_from_cif() "
            "(pip install biopython)."
        )

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("target", str(cif_path))

    ppb = PPBuilder()
    for model in structure:
        if chain_id not in model:
            continue
        chain = model[chain_id]
        peptides = ppb.build_peptides(chain)
        if not peptides:
            continue
        # If multiple fragments, just concatenate them
        seq = "".join(str(pp.get_sequence()) for pp in peptides)
        if seq:
            return seq

    return ""


# ---------------------------
# YAML builder (Boltz-2 schema)
# ---------------------------

def make_yaml_for_pair(
    target_seq: str,
    binder_seq: str,
    target_chain: str = "A",
    binder_chain: str = "B",
    template_cif: Optional[Path] = None,
) -> str:
    """
    Build a minimal, schema-correct Boltz-2 YAML for a protein–protein complex
    with an optional structural template for the target chain.

    Format based on Boltz-2 prediction docs. :contentReference[oaicite:1]{index=1}
    """
    data: Dict = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": target_chain,
                    "sequence": target_seq,
                }
            },
            {
                "protein": {
                    "id": binder_chain,
                    "sequence": binder_seq,
                }
            },
        ],
    }

    if template_cif is not None:
        data["templates"] = [
            {
                "cif": str(template_cif.resolve()),
                "chain_id": target_chain,
            }
        ]

    return yaml.safe_dump(data, sort_keys=False)


# ---------------------------
# ipSAE utilities
# ---------------------------

def find_first(root: Path, want) -> Optional[Path]:
    """
    Walk the tree under `root` and return the first path where want(path, name) is True.
    """
    for p in root.rglob("*"):
        name = p.name.lower()
        try:
            if want(p, name):
                return p
        except Exception:
            pass
    return None


# ---------------------------
# Main
# ---------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Batch Boltz-2 co-folding for binder–target complexes."
    )
    ap.add_argument(
        "--fasta",
        required=True,
        help="FASTA file containing binder sequences.",
    )
    ap.add_argument(
        "--target_cif",
        required=True,
        help="Target structure (e.g. 2VSM.cif).",
    )
    ap.add_argument(
        "--target_chain",
        required=True,
        help="Target chain ID in the CIF (e.g. A).",
    )
    ap.add_argument(
        "--yaml_out_dir",
        required=True,
        help="Directory where YAML files will be written.",
    )
    ap.add_argument(
        "--run_root_dir",
        required=True,
        help="Root directory where per-YAML Boltz runs will live.",
    )
    ap.add_argument(
        "--logs_dir",
        required=False,
        help="Optional directory for extra pipeline logs (not required).",
    )

    # Boltz flags
    ap.add_argument(
        "--use_msa_server",
        action="store_true",
        help="Pass --use_msa_server to `boltz predict`.",
    )
    ap.add_argument(
        "--no_kernels",
        action="store_true",
        help="Pass --no_kernels to `boltz predict`.",
    )
    ap.add_argument(
        "--override",
        action="store_true",
        help="Pass --override to `boltz predict`.",
    )

    # ipSAE
    ap.add_argument(
        "--ipsae_py",
        help="Path to ipsae.py (if provided, ipSAE will be run).",
    )
    ap.add_argument(
        "--ipsae_out_root",
        help="Output root for ipSAE results.",
    )

    O_DIR = Path("./boltzgen_out_1000/protein_vs_2vsm_protein-anything/")
    INPUT_FASTA = Path(f"{O_DIR}/fall_complexes.fasta")

    if len(sys.argv) == 1:
        sys.argv.extend([
            f"--fasta", str(INPUT_FASTA),
            "--target_cif", "./inputs/2VSM.cif",
            "--target_chain", "A",
            "--yaml_out_dir", f"{O_DIR}/yaml_batch",
            "--run_root_dir", f"{O_DIR}/runs",
            "--logs_dir", f"{O_DIR}/logs",
            "--use_msa_server",
            "--no_kernels",
            "--override",
            "--ipsae_py", "ipsae.py",
            "--ipsae_out_root", f"{O_DIR}/ipSAE",
        ])

    args = ap.parse_args()

    binder_fasta = Path(args.fasta)
    if not binder_fasta.exists():
        log(f"[err] binder FASTA not found: {binder_fasta}", stream=sys.stderr)
        sys.exit(2)

    target_cif = Path(args.target_cif)
    if not target_cif.exists():
        log(f"[err] target CIF not found: {target_cif}", stream=sys.stderr)
        sys.exit(2)

    # ---- target sequence from CIF ----
    target_seq = get_target_sequence_from_cif(target_cif, args.target_chain)
    if not target_seq:
        log(
            f"[err] could not extract sequence for chain {args.target_chain} from {target_cif}",
            stream=sys.stderr,
        )
        sys.exit(2)

    yaml_dir = Path(args.yaml_out_dir)
    yaml_dir.mkdir(parents=True, exist_ok=True)

    run_root = Path(args.run_root_dir)
    run_root.mkdir(parents=True, exist_ok=True)

    # ---- build YAMLs for each binder ----
    yaml_paths: List[Path] = []
    seen_seqs: set = set()

    base_name = re.sub(r"_T_[^_]+$", "", binder_fasta.stem)

    for hdr, seq in read_fasta(binder_fasta):
        if not seq:
            continue
        if seq in seen_seqs:
            # global dedup on exact sequence
            continue
        seen_seqs.add(seq)

        T = infer_temperature_from_header(hdr)
        uid = short_hash(seq)
        yaml_base = f"{base_name}__T_{T}_{uid}"

        yp = yaml_dir / f"{yaml_base}.yaml"
        yaml_text = make_yaml_for_pair(
            target_seq=target_seq,
            binder_seq=seq,
            target_chain=args.target_chain,
            binder_chain="B",
            template_cif=target_cif,
        )
        yp.write_text(yaml_text)
        yaml_paths.append(yp)

    if not yaml_paths:
        log("[err] No YAMLs were created; aborting.", stream=sys.stderr)
        sys.exit(2)

    log(f"[ok] Wrote {len(yaml_paths)} YAMLs to {yaml_dir}")

    # ---- run Boltz-2 ----
    failures = 0
    for yp in yaml_paths:
        base = yp.stem
        run_dir = run_root / base
        run_dir.mkdir(parents=True, exist_ok=True)
        run_log = run_dir / "boltz_run.log"

        cmd = ["boltz", "predict", str(yp.resolve()), "--write_full_pae"]
        if args.use_msa_server:
            cmd.append("--use_msa_server")
        if args.no_kernels:
            cmd.append("--no_kernels")
        if args.override:
            cmd.append("--override")

        log(f"[run] boltz predict {yp}", stream=sys.stdout)
        rc = run_cmd(cmd, cwd=run_dir, log_path=run_log)
        if rc != 0:
            log(f"[err] boltz failed for {base} (exit {rc})", stream=sys.stderr)
            failures += 1

    log(f"[done] Boltz runs complete. failures={failures}")

    # ---- ipSAE (optional) ----
    if args.ipsae_py and args.ipsae_out_root:
        ipsae_py = Path(args.ipsae_py)
        ipsae_root = Path(args.ipsae_out_root)
        ipsae_root.mkdir(parents=True, exist_ok=True)

        summary_rows = []

        for yp in yaml_paths:
            base = yp.stem
            run_dir = run_root / base

            # Boltz outputs under ./boltz_results_<basename>/predictions/<basename>/
            cif = find_first(
                run_dir,
                lambda p, n: (
                    p.is_file()
                    and (n.endswith(".cif") or n.endswith(".mmcif") or n.endswith(".pdb"))
                    and ("model_0" in n or "output_model_0" in n)
                ),
            )
            pae = find_first(
                run_dir,
                lambda p, n: (
                    p.is_file()
                    and (n.endswith(".npz") or n.endswith(".json"))
                    and ("pae_" in n and ("model_0" in n or "output_model_0" in n))
                ),
            )

            out_txt = ipsae_root / f"{base}.txt"

            if cif and pae:
                cmd = [
                    sys.executable,
                    str(ipsae_py),
                    str(pae),
                    str(cif),
                    "10",
                    "10",
                ]
                rc = run_cmd(cmd, cwd=ipsae_root, log_path=out_txt)
                summary_rows.append((base, str(cif), str(pae), rc))
            else:
                log(f"[ipsae] {base} (missing cif or pae npz/json)")

        tsv = ipsae_root / "ipsae_all.tsv"
        with tsv.open("w") as f:
            f.write("base\tcif\tpae\texit_code\n")
            for r in summary_rows:
                f.write("\t".join(map(str, r)) + "\n")
        log(f"[ok] ipSAE summary: {tsv}")


if __name__ == "__main__":
    main()
