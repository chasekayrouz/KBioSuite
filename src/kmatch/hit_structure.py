"""Reconcile public full-atom assemblies with the exact local Foldseek trace."""

import re

import gemmi
import numpy as np

from .models import Progress
from .pdb_data import assembly_file, pdb_identity

CONTACT_ANGSTROM = 4.0


def residue_id(residue):
    return {
        "name": residue.name,
        "number": residue.seqid.num,
        "insertion_code": residue.seqid.icode.strip(),
        "label_sequence_id": residue.label_seq,
    }


def major_conformer(residue, preferred=None):
    copy = residue.clone()
    labels = sorted({a.altloc for a in copy if a.has_altloc()})
    chosen = preferred or max(
        labels, key=lambda label: sum(a.occ for a in copy if a.altloc == label), default=None
    )
    for i in reversed(range(len(copy))):
        atom = copy[i]
        if atom.element.is_hydrogen or (atom.has_altloc() and atom.altloc != chosen):
            del copy[i]
        else:
            atom.altloc = "\0"
    return copy


def ca_residues(chain):
    return [r for r in chain if any(a.name == "CA" and a.element.name == "C" for a in r)]


def xyz(residues):
    return np.array([list(next(a for a in r if a.name == "CA").pos) for r in residues])


def sequence(residues):
    return "".join(
        gemmi.find_tabulated_residue(r.name).one_letter_code.upper().strip() or "X"
        for r in residues
    )


def trace_mapping(trace_sequence, full_sequence):
    """Allow extra observed full-atom residues, never substitute a database residue."""
    if trace_sequence == full_sequence:
        return list(range(len(trace_sequence)))
    alignment = gemmi.align_string_sequences(
        list(trace_sequence), list(full_sequence), [], gemmi.AlignmentScoring()
    )
    qi = ti = 0
    mapping = []
    for count, operation in re.findall(r"(\d+)([MID])", alignment.cigar_str()):
        count = int(count)
        if operation == "I":
            raise ValueError("Full-atom structure is missing residues from the database trace.")
        if operation == "D":
            ti += count
        else:
            if trace_sequence[qi : qi + count] != full_sequence[ti : ti + count]:
                raise ValueError("Full-atom and database sequences differ.")
            mapping.extend(range(ti, ti + count))
            qi += count
            ti += count
    if len(mapping) != len(trace_sequence):
        raise ValueError("Full-atom structure cannot cover the complete database trace.")
    return mapping


def rigid_registration(source, reference):
    if len(source) < 3 or np.linalg.matrix_rank(source - source.mean(0)) < 2:
        raise ValueError("Too few independent Cα positions to establish the source frame.")
    a, b = source.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((source - a).T @ (reference - b))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    row_rotation = u @ correction @ vt
    rotation, translation = row_rotation.T, b - a @ row_rotation
    distances = np.linalg.norm(source @ rotation.T + translation - reference, axis=1)
    return rotation, translation, float(np.sqrt(np.mean(distances**2))), float(max(distances))


def transform_chain(chain, rotation, translation):
    for residue in chain:
        for atom in residue:
            atom.pos = list(np.asarray(list(atom.pos)) @ rotation.T + translation)


def full_atom_superposition(config, combined, trace, hit, allow_network, cancel, progress):
    path, source = assembly_file(config, hit.target, allow_network, cancel, progress)
    progress(Progress("export", "Matching full-atom PDB coordinates to the local database frame…"))
    structure = gemmi.read_structure(str(path))
    _, _, requested_chain = pdb_identity(hit.target)
    trace_residues = ca_residues(trace[0][0])
    trace_seq, trace_xyz = sequence(trace_residues), xyz(trace_residues)
    candidates = []
    for model_index, model in enumerate(structure):
        for chain_index, chain in enumerate(model):
            cancel.check()
            polymer = gemmi.Chain(chain.name)
            for residue in chain:
                if residue.entity_type == gemmi.EntityType.Polymer:
                    polymer.add_residue(major_conformer(residue))
            residues = ca_residues(polymer)
            if not residues:
                continue
            # Author names take precedence; some generated copies use a suffix or label ID.
            labels = {r.subchain for r in residues}
            if not (
                chain.name == requested_chain
                or requested_chain in labels
                or re.sub(r"-\d+$", "", requested_chain) == chain.name
            ):
                continue
            try:
                mapping = trace_mapping(trace_seq, sequence(residues))
                positions = xyz(residues)[mapping]
                rotation, translation, rmsd, maximum = rigid_registration(positions, trace_xyz)
                # Foldseek may retain a different backbone alternate than the most occupied one.
                # Select the experimentally present Cα alternate closest to its database position,
                # then apply the same strict frame check to the complete resolved chain.
                source_ca = [
                    r
                    for r in chain
                    if r.entity_type == gemmi.EntityType.Polymer
                    and any(a.name == "CA" and a.element.name == "C" for a in r)
                ]
                preferred = {}
                for db_index, source_index in enumerate(mapping):
                    original = source_ca[source_index]
                    choices = [a for a in original if a.name == "CA" and a.element.name == "C"]
                    if len(choices) > 1:
                        nearest = min(
                            choices,
                            key=lambda a: np.linalg.norm(
                                np.asarray(list(a.pos)) @ rotation.T
                                + translation
                                - trace_xyz[db_index]
                            ),
                        )
                        preferred[(original.seqid.num, original.seqid.icode, original.name)] = (
                            nearest.altloc
                        )
                if preferred:
                    polymer = gemmi.Chain(chain.name)
                    for residue in chain:
                        if residue.entity_type == gemmi.EntityType.Polymer:
                            key = (residue.seqid.num, residue.seqid.icode, residue.name)
                            polymer.add_residue(major_conformer(residue, preferred.get(key)))
                    residues = ca_residues(polymer)
                    positions = xyz(residues)[mapping]
                    rotation, translation, rmsd, maximum = rigid_registration(positions, trace_xyz)
            except ValueError:
                continue
            if rmsd <= 0.15 and maximum <= 0.5:
                raw_rmsd = float(np.sqrt(np.mean(np.sum((positions - trace_xyz) ** 2, axis=1))))
                candidates.append(
                    (
                        rmsd,
                        raw_rmsd,
                        model_index,
                        chain_index,
                        polymer,
                        residues,
                        mapping,
                        rotation,
                        translation,
                        maximum,
                    )
                )
    if not candidates:
        raise RuntimeError(
            "The downloaded PDB chain does not reproduce the saved database trace "
            "(required Cα RMSD ≤0.15 Å and maximum deviation ≤0.5 Å). "
            "No full-atom alignment was exported. Use the local Cα option or "
            "update the database and run a new search."
        )
    # Quantization differences should not outweigh exact assembly placement for identical copies.
    candidates.sort(key=lambda c: (round(c[0], 2), c[1], c[2], c[3]))
    rmsd, _, model_index, chain_index, polymer, residues, mapping, reg_u, reg_t, maximum = (
        candidates[0]
    )
    model = structure[model_index]
    original_chain = model[chain_index].name
    # Contacts are measured in the experimental source frame, before either rigid transform.
    heavy_xyz, heavy_residue = [], []
    for i, residue in enumerate(residues):
        for atom in residue:
            if atom.occ > 0:
                heavy_xyz.append(list(atom.pos))
                heavy_residue.append(i)
    heavy_xyz, heavy_residue = np.array(heavy_xyz), np.array(heavy_residue)
    contact_ligands, ligand_chain, contacting = [], gemmi.Chain("L"), set()
    for source_chain in model:
        for residue in source_chain:
            cancel.check()
            # Excludes water, protein/nucleic-acid polymers and branched glycans.
            if residue.entity_type != gemmi.EntityType.NonPolymer or residue.is_water():
                continue
            ligand = major_conformer(residue)
            ligand_atoms = [a for a in ligand if a.occ > 0]
            contacts = set()
            nearest = float("inf")
            for atom in ligand_atoms:
                distance = np.linalg.norm(heavy_xyz - np.array(list(atom.pos)), axis=1)
                nearest = min(nearest, float(distance.min()))
                contacts.update(map(int, heavy_residue[distance <= CONTACT_ANGSTROM]))
            if not contacts:
                continue
            number = len(contact_ligands) + 1
            contact_ligands.append(
                {
                    "component_id": residue.name,
                    "source_chain": source_chain.name,
                    "source_residue": residue_id(residue),
                    "export_chain": "L",
                    "export_residue_number": number,
                    "heavy_atoms": len(ligand_atoms),
                    "minimum_protein_distance_angstrom": nearest,
                    "contacting_hit_residues": [residue_id(residues[i]) for i in sorted(contacts)],
                }
            )
            ligand.seqid = gemmi.SeqId(number, " ")
            ligand.subchain = "L"
            ligand.entity_id = ""
            ligand.label_seq = None
            ligand_chain.add_residue(ligand)
            contacting.update(contacts)
    # Map through saved Foldseek alignment columns, never through residue-number equality.
    from .alignment import aligned_pairs

    qchain = combined[0]["Q"]
    _, qseq, _ = ca_chain_from(qchain)
    pairs = aligned_pairs(hit, qseq, trace_seq)
    query_residues = ca_residues(qchain)
    mapped_contacts = []
    for qi, ti in pairs:
        full_index = mapping[int(ti)]
        if full_index in contacting:
            mapped_contacts.append(
                {
                    "query_position": int(qi) + 1,
                    "database_target_position": int(ti) + 1,
                    "query_residue": residue_id(query_residues[int(qi)]),
                    "hit_residue": residue_id(residues[full_index]),
                }
            )
    # Snapshot coordinates before mutating the Gemmi model container.
    q_positions = xyz(query_residues)[pairs[:, 0]]
    # Registration only restores the database frame; the retained query/hit alignment is unchanged.
    rotation = np.asarray(hit.rotation) @ reg_u
    translation = np.asarray(hit.rotation) @ reg_t + np.asarray(hit.translation)
    polymer.name = "T"
    for residue in polymer:
        residue.subchain = "T"
        residue.entity_id = ""
    transform_chain(polymer, rotation, translation)
    transform_chain(ligand_chain, rotation, translation)
    del combined[0]["T"]
    combined[0].add_chain(polymer)
    if len(ligand_chain):
        combined[0].add_chain(ligand_chain)
    combined.entities.clear()
    combined.setup_entities()
    combined.assign_label_seq_id()
    t_positions = xyz(ca_residues(combined[0]["T"]))[[mapping[int(ti)] for ti in pairs[:, 1]]]
    full_rmsd = float(np.sqrt(np.mean(np.sum((q_positions - t_positions) ** 2, axis=1))))
    return {
        "query_aligned_ca_rmsd_angstrom": full_rmsd,
        "source": source,
        "source_model": model.num,
        "source_chain": original_chain,
        "representation": "full-atom protein; contacting non-water small molecules and ions",
        "conformer_policy": (
            "database-matching Cα alternate; otherwise highest summed occupancy per residue; "
            "no hydrogens"
        ),
        "registration": {
            "rotation": reg_u.tolist(),
            "translation": reg_t.tolist(),
            "ca_rmsd_angstrom": rmsd,
            "ca_max_deviation_angstrom": maximum,
            "ca_pairs": len(mapping),
            "purpose": "source→database frame only",
        },
        "source_to_query": {"rotation": rotation.tolist(), "translation": translation.tolist()},
        "contact_cutoff_angstrom": CONTACT_ANGSTROM,
        "contact_definition": (
            "any protein/ligand heavy-atom distance ≤4 Å; proximity, not a bond assignment"
        ),
        "ligands": contact_ligands,
        "hit_contact_residues": [residue_id(residues[i]) for i in sorted(contacting)],
        "aligned_query_contacts": mapped_contacts,
        "unmapped_hit_contact_count": len(contacting - {mapping[int(ti)] for _, ti in pairs}),
    }


def ca_chain_from(chain):
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    model.add_chain(chain.clone())
    structure.add_model(model)
    from .alignment import ca_chain

    return ca_chain(structure, "Query")
