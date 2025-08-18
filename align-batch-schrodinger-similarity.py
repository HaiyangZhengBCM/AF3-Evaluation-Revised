#!/usr/bin/env python3
import os
import sys
import argparse
import traceback
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter

import pandas as pd
from schrodinger.structutils.structalign2 import cealign
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils import rmsd as rmsd_mod
from schrodinger.structutils.analyze import evaluate_asl

from schrodinger.application.canvas.fingerprintgui import CanvasFingerprintGeneratorGUI
from schrodinger.application.canvas.similarity import CanvasFingerprintSimilarity

# ----------------------
# Canvas logger & globals
# ----------------------
_canvas_logger = logging.getLogger("canvas_fpsim")
if not _canvas_logger.handlers:
    _canvas_logger.addHandler(logging.NullHandler())
_canvas_logger.setLevel(logging.WARNING)  # change to INFO/DEBUG if you want chatter

_FP_GEN = CanvasFingerprintGeneratorGUI(logger=_canvas_logger)  # Default: Linear, 32-bit
_SIM    = CanvasFingerprintSimilarity(logger=_canvas_logger)    # Default metric: Tanimoto

# ----------------------
# Core utilities
# ----------------------

def align_and_save_complex(pdb_id, ref_path, mob_path, out_dir):
    """CE-align mobile structure onto reference (Cα atoms) and save both."""
    os.makedirs(out_dir, exist_ok=True)
    ref = next(StructureReader(ref_path))
    mob = next(StructureReader(mob_path))

    alignment = cealign(ref, mob, transform_mobile=True)
    rmsd_bb = alignment.rmsd
    print(f"[{pdb_id}] Cα RMSD = {rmsd_bb:.3f} Å")

    ref_out = os.path.join(out_dir, f"aligned_ref_{pdb_id}.pdb")
    mob_out = os.path.join(out_dir, f"aligned_mob_{pdb_id}.pdb")
    with StructureWriter(ref_out) as w_ref:
        w_ref.append(ref)
    with StructureWriter(mob_out) as w_mob:
        w_mob.append(mob)

    return ref_out, mob_out, rmsd_bb


def split_ligands_from_pdb(pdb_file, out_dir, pdb_code, resname_filter=None, min_heavy=3):
    """
    Extract each disconnected ligand component (excluding waters) into its own PDB.
    Splits by connectivity (bond graph), not residue labels. Optionally filter by residue name.
    """
    os.makedirs(out_dir, exist_ok=True)
    st = next(StructureReader(pdb_file))

    lig_idxs = set(evaluate_asl(st, 'ligand and not water'))
    if not lig_idxs:
        return []

    if resname_filter:
        keep = set()
        for a in st.atom:
            if a.index in lig_idxs:
                r = a.getResidue()
                rname = (getattr(r, 'pdbres', None) or getattr(r, 'name', '') or '').strip()
                if rname == resname_filter:
                    keep.add(a.index)
        lig_idxs = keep
        if not lig_idxs:
            return []

    neighbors = {i: set() for i in lig_idxs}
    for b in st.bond:
        a1, a2 = b.atom1.index, b.atom2.index
        if a1 in lig_idxs and a2 in lig_idxs:
            neighbors[a1].add(a2)
            neighbors[a2].add(a1)

    visited, components = set(), []
    for i in lig_idxs:
        if i in visited:
            continue
        stack = [i]
        visited.add(i)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in neighbors[u]:
                if v not in visited:
                    visited.add(v)
                    stack.append(v)
        components.append(sorted(comp))

    def _heavy_atom_count_in(struct, atom_indices):
        return sum(1 for idx in atom_indices if struct.atom[idx].element != 'H')

    ligand_files = []
    for k, atom_idxs in enumerate(components, start=1):
        if _heavy_atom_count_in(st, atom_idxs) < min_heavy:
            continue
        sub = st.extract(atom_idxs)
        # Majority labels (for filenames only)
        names, chains, resnums = [], [], []
        for a in sub.atom:
            r = a.getResidue()
            names.append((getattr(r, 'pdbres', None) or getattr(r, 'name', '') or '').strip())
            chains.append((getattr(r, 'chain', '') or '').strip())
            rn = getattr(r, 'resnum', None)
            if rn is not None:
                resnums.append(rn)

        def mode(xs, default):
            c = Counter(x for x in xs if x)
            return c.most_common(1)[0][0] if c else default

        resname = mode(names, 'LIG')
        chain = mode(chains, 'X')
        resnum = mode(resnums, k)

        out_path = os.path.join(out_dir, f"ligand_{resnum}_{resname}_{chain}_{pdb_code}_{k}.pdb")
        with StructureWriter(out_path) as w:
            w.append(sub)
        ligand_files.append(out_path)

    print(f"[{pdb_code}] Split {len(components)} ligand components; wrote {len(ligand_files)} files.")
    return ligand_files


def compute_rmsd(gt_file, pr_file):
    """Chemistry-aware, direct (in-place) heavy-atom RMSD."""
    ref = next(StructureReader(gt_file))
    test = next(StructureReader(pr_file))
    conf = rmsd_mod.ConformerRmsd(ref, test, asl_expr='NOT atom.element H', in_place=True)
    conf.use_symmetry = True
    conf.renumber_structures = True
    conf.use_heavy_atom_graph = True
    return float(conf.calculate())


def heavy_atom_count(path):
    st = next(StructureReader(path))
    return sum(1 for a in st.atom if a.element != 'H')


def find_best_ligand_match(gt_files, pred_file):
    """Return the GT ligand with the lowest RMSD to the predicted ligand."""
    best_val = float('inf')
    pred_hvy = heavy_atom_count(pred_file)
    for gt in gt_files:
        gt_hvy = heavy_atom_count(gt)
        if abs(gt_hvy - pred_hvy) > max(2, int(0.1 * max(pred_hvy, 1))):
            print(f"  Skipping {os.path.basename(gt)} (heavy atoms {gt_hvy} vs {pred_hvy})")
            continue
        try:
            rmsd = compute_rmsd(gt, pred_file)
            print(f"  {os.path.basename(gt)} RMSD = {rmsd:.3f} Å")
        except Exception as e:
            print(f"  Skipping {os.path.basename(gt)} due to mapping error: {e}")
            continue
        if rmsd < best_val:
            best_val = rmsd
    if best_val == float('inf'):
        raise RuntimeError("No comparable ligands found (mapping/size mismatch).")
    return best_val

# ----------------------
# Canvas fingerprint utilities
# ----------------------

def _canvas_fp_from_pdb(path):
    st = next(StructureReader(path))
    return _FP_GEN.generate(st)

def compute_tanimoto_scores_canvas(gt_files, pred_file):
    """
    Compute Canvas Tanimoto similarity between predicted ligand and each GT ligand.
    Returns: (best_score, per_gt_list)
    """
    pred_fp = _canvas_fp_from_pdb(pred_file)
    per_gt, best = [], 0.0
    for gt in gt_files:
        try:
            gt_fp = _canvas_fp_from_pdb(gt)
            sim = float(_SIM.calculateSimilarity(pred_fp, gt_fp))
            per_gt.append((gt, sim))
            if sim > best:
                best = sim
        except Exception:
            per_gt.append((gt, None))
    return best, per_gt

# ----------------------
# Batch / multiprocessing
# ----------------------

def process_pair(pdb_id, real_path, pred_path, work_root):
    """Do everything for one pair; return minimal dict for CSV."""
    work_dir = os.path.join(work_root, pdb_id)
    os.makedirs(work_dir, exist_ok=True)
    try:
        ref_aln, pred_aln, ca_rmsd = align_and_save_complex(pdb_id, real_path, pred_path, work_dir)

        gt_ligs = split_ligands_from_pdb(ref_aln, os.path.join(work_dir, 'gt_ligs'), pdb_id)
        pr_ligs = split_ligands_from_pdb(pred_aln, os.path.join(work_dir, 'pr_ligs'), pdb_id)

        ligand_rmsd = "N/A"
        match_found = False
        tanimoto_best = "N/A"
        tanimoto_per_gt = []

        if pr_ligs:
            try:
                ligand_rmsd = round(find_best_ligand_match(gt_ligs, pr_ligs[0]), 3)
                match_found = True
            except Exception as e:
                print(f"[{pdb_id}] Ligand comparison failed: {e}")
                try:
                    tbest, per_gt = compute_tanimoto_scores_canvas(gt_ligs, pr_ligs[0])
                    tanimoto_best = round(tbest, 3)
                    tanimoto_per_gt = [(os.path.basename(g), (None if s is None else round(s, 3)))
                                       for g, s in per_gt]
                except Exception as e2:
                    print(f"[{pdb_id}] Canvas Tanimoto calculation failed: {e2}")
        else:
            print(f"[{pdb_id}] No ligands in predicted structure.")

        return {
            "pdb_id": pdb_id,
            "protein_ca_rmsd": round(float(ca_rmsd), 3) if ca_rmsd is not None else "N/A",
            "ligand_rmsd": ligand_rmsd,
            "match_found": match_found,
            "gt_lig_count": len(gt_ligs),
            "pr_lig_count": len(pr_ligs),
            "tanimoto_best": tanimoto_best,
            "tanimoto_per_gt": tanimoto_per_gt
        }
    except Exception as e:
        print(f"[{pdb_id}] ERROR: {e}")
        print(traceback.format_exc())
        return {
            "pdb_id": pdb_id,
            "protein_ca_rmsd": "N/A",
            "ligand_rmsd": "N/A",
            "match_found": False,
            "gt_lig_count": 0,
            "pr_lig_count": 0,
            "tanimoto_best": "N/A",
            "tanimoto_per_gt": []
        }


def find_pairs(real_dir, pred_dir):
    """
    Match real *_prep.pdb with predicted *_model_prep.pdb by PDB code prefix.
    Returns list of (pdb_id, real_path, pred_path).
    """
    pairs = []
    real_files = [f for f in os.listdir(real_dir) if f.endswith("_prep.pdb") and not f.endswith("_model_prep.pdb")]
    for f in real_files:
        pdb_id = f[:-9]  # strip "_prep.pdb"
        real_path = os.path.join(real_dir, f)
        pred_path = os.path.join(pred_dir, f"{pdb_id}_model_prep.pdb")
        if os.path.exists(pred_path):
            pairs.append((pdb_id, real_path, pred_path))
        else:
            print(f"[WARN] Predicted file missing for {pdb_id}: {pred_path}")
    return pairs


def main():
    ap = argparse.ArgumentParser(description="Batch CE-align, ligand RMSD, and Canvas Tanimoto.")
    ap.add_argument("--real_dir", required=True)
    ap.add_argument("--pred_dir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tanimoto_csv", default=None, help="Optional: write per-GT Tanimoto scores for failures")
    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)

    pairs = find_pairs(args.real_dir, args.pred_dir)
    if not pairs:
        print("No matching pairs found.")
        sys.exit(1)

    print(f"Found {len(pairs)} pairs. Processing with {args.workers} workers...")

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = [ex.submit(process_pair, pid, r, p, args.work_dir) for (pid, r, p) in pairs]
        for fut in as_completed(futures):
            rows.append(fut.result())

    df = pd.DataFrame(rows, columns=[
        "pdb_id", "protein_ca_rmsd", "ligand_rmsd", "match_found",
        "gt_lig_count", "pr_lig_count", "tanimoto_best", "tanimoto_per_gt"
    ])
    df.to_csv(args.out_csv, index=False)

    # Percentage of pairs (with predicted ligands) where no RMSD match was found
    considered = df[df["pr_lig_count"] > 0]
    if len(considered) > 0:
        failures = considered[~considered["match_found"]]
        failure_pct = 100.0 * len(failures) / len(considered)
        print(f"\nPairs with predicted ligands: {len(considered)}; no-match cases: {len(failures)} ({failure_pct:.1f}%)")
    else:
        print("\nNo predicted ligands found in any pair.")

    # Optional Tanimoto CSV for failures
    if args.tanimoto_csv:
        tan_rows = []
        for _, r in df[df["match_found"] == False].iterrows():
            for entry in (r["tanimoto_per_gt"] or []):
                gt_name, sim = entry
                tan_rows.append({
                    "pdb_id": r["pdb_id"],
                    "tanimoto_best": r["tanimoto_best"],
                    "gt_ligand_file": gt_name,
                    "tanimoto_to_pred": sim
                })
        if tan_rows:
            tan_df = pd.DataFrame(tan_rows, columns=["pdb_id", "tanimoto_best", "gt_ligand_file", "tanimoto_to_pred"])
            tan_df.to_csv(args.tanimoto_csv, index=False)
            print(f"Wrote Tanimoto details for failures to {args.tanimoto_csv}")
        else:
            print("No Tanimoto details to write.")

    print(f"\nSaved results to {args.out_csv}")
    print(df[["pdb_id", "protein_ca_rmsd", "ligand_rmsd", "match_found", "tanimoto_best"]]
          .head(min(10, len(df))).to_string(index=False))


if __name__ == "__main__":
    main()
