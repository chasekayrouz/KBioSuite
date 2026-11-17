"""Run by an external PyMOL installation with -cqk; no Kmatch/Python ABI coupling."""

import json
import math
import os
import re
from pathlib import Path

import pymol
from pymol import cmd


def coordinates(name):
    return [list(atom.coord) for atom in cmd.get_model(name).atom]


def main():
    request = json.loads(Path(os.environ["KMATCH_ALIGNMENT_REQUEST"]).read_text())
    provenance = request["provenance"]
    query, target = (provenance["pymol_objects"][key] for key in ("query", "hit"))
    cmd.reinitialize()
    cmd.load(request["coordinates"], "kmatch_input", quiet=1)
    cmd.create(query, "kmatch_input and chain Q")
    cmd.create(target, "kmatch_input and chain T")
    full = provenance.get("full_atom")
    objects = [query, target]
    ligand_objects = []
    for ligand in (full or {}).get("ligands", []):
        code = re.sub(r"[^A-Za-z0-9_]", "_", ligand["component_id"])
        number = ligand["export_residue_number"]
        name = f"ligand_{code}_{number}"
        cmd.create(name, f"kmatch_input and chain L and resi {number}")
        objects.append(name)
        ligand_objects.append((name, ligand["heavy_atoms"]))
    cmd.delete("kmatch_input")
    if any(not cmd.count_atoms(name) for name in objects):
        raise RuntimeError("Each exported object must contain coordinates.")
    cmd.hide("everything", "all")
    colors = {
        "kmatch_query": [0.545, 0.768, 0.616],
        "kmatch_hit": [0.714, 0.741, 0.757],
        "kmatch_query_site": [0.153, 0.467, 0.298],
        "kmatch_hit_site": [0.337, 0.408, 0.459],
        "kmatch_ligand": [0.871, 0.635, 0.231],
        "kmatch_nitrogen": [0.19, 0.39, 0.81],
        "kmatch_oxygen": [0.85, 0.25, 0.24],
        "kmatch_sulfur": [0.83, 0.70, 0.18],
        "kmatch_ion": [0.70, 0.36, 0.65],
    }
    for name, rgb in colors.items():
        cmd.set_color(name, rgb)
    cmd.color("kmatch_query", query)
    cmd.color("kmatch_hit", target)
    cmd.set("cartoon_color", "kmatch_query", query)
    cmd.set("cartoon_color", "kmatch_hit", target)
    cmd.cartoon("automatic", query + " or " + target)
    if full:
        cmd.dss(query + " or " + target)
    else:
        # Explicit local-only fallback cannot reconstruct secondary structure/side chains.
        cmd.set("cartoon_trace_atoms", 1, target)
        cmd.cartoon("tube", target)
    cmd.show("cartoon", query + " or " + target)
    cmd.set("valence", 0)
    cmd.set("stick_radius", 0.18)
    cmd.set("stick_transparency", 0)

    def residue_selection(name, residues):
        wanted = {(str(r["number"]) + r["insertion_code"], r["name"]) for r in residues}
        indices = [
            str(a.index)
            for a in cmd.get_model(name).atom
            if (a.resi, a.resn) in wanted
            and a.name not in {"N", "C", "O", "OXT"}
            and a.symbol not in {"H", "D"}
        ]
        return name + " and index " + "+".join(indices) if indices else "none"

    def element_colors(selection):
        for element, color in (
            ("N", "kmatch_nitrogen"),
            ("O", "kmatch_oxygen"),
            ("S", "kmatch_sulfur"),
            ("P", "kmatch_ion"),
        ):
            cmd.color(color, f"({selection}) and elem {element}")

    site_selections = {}
    if full and full["hit_contact_residues"]:
        site_selections["hit_ligand_contacts"] = residue_selection(
            target, full["hit_contact_residues"]
        )
        site_selections["query_aligned_contacts"] = residue_selection(
            query, [p["query_residue"] for p in full["aligned_query_contacts"]]
        )
        for selection_name, selection in site_selections.items():
            cmd.select(selection_name, selection, enable=0)
            cmd.show("sticks", selection_name)
            cmd.color(
                "kmatch_hit_site" if selection_name.startswith("hit_") else "kmatch_query_site",
                selection_name,
            )
            element_colors(selection_name)
        cmd.set("cartoon_transparency", 0.60, query)
        cmd.set("cartoon_transparency", 0.60, target)
    for name, count in ligand_objects:
        cmd.color("kmatch_ligand", name)
        element_colors(name)
        cmd.show("sticks", name)
        if count == 1:
            # A free metal has no bond to draw as a stick; show a compact sphere instead.
            cmd.color("kmatch_ion", name)
            cmd.show("spheres", name)
            cmd.set("sphere_scale", 0.35, name)
    cmd.hide("everything", "hydro")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("orthoscopic", 1)
    cmd.orient("all")
    cmd.zoom("all", buffer=3)
    pymol.session.kmatch_alignment = provenance
    expected = {name: coordinates(name) for name in objects}
    cmd.save(request["output"], format="pse", quiet=1)
    cmd.reinitialize()
    del pymol.session.kmatch_alignment
    cmd.load(request["output"], quiet=1)
    if set(cmd.get_names("objects")) != set(objects):
        raise RuntimeError("PyMOL session did not retain the named query, hit and ligand objects.")
    for name, before in expected.items():
        after = coordinates(name)
        if len(before) != len(after) or any(
            not math.isclose(a, b, abs_tol=0.0001, rel_tol=0)
            for index, p in enumerate(before)
            for a, b in ((p[axis], after[index][axis]) for axis in range(3))
        ):
            raise RuntimeError("PyMOL changed coordinates during session round-trip.")
    if getattr(pymol.session, "kmatch_alignment", None) != provenance:
        raise RuntimeError("PyMOL session did not retain alignment provenance.")
    Path(request["verification"]).write_text(
        json.dumps(
            {
                "roundtrip_verified": True,
                "pymol_version": cmd.get_version()[0],
                "objects": {name: len(value) for name, value in expected.items()},
                "coordinates_unchanged_on_reload": True,
                "provenance_retained": True,
                "presentation": {
                    "colors": colors,
                    "ribbon_transparency": float(cmd.get("cartoon_transparency", query)),
                    "query_stick_atoms": cmd.count_atoms(query + " and rep sticks"),
                    "hit_stick_atoms": cmd.count_atoms(target + " and rep sticks"),
                    "ligand_objects": [name for name, _ in ligand_objects],
                    "valence": cmd.get_setting_int("valence"),
                    "visible_hydrogen_atoms": cmd.count_atoms(
                        "hydro and (rep sticks or rep spheres or rep lines)"
                    ),
                },
            },
            indent=2,
        )
        + "\n"
    )


try:
    main()
except Exception as exc:
    print("Kmatch PyMOL export failed: " + str(exc), flush=True)
    cmd.quit(1)
else:
    cmd.quit(0)
