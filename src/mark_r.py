
import io
import os
import random
import re
from tkinter import Image
from rdkit import Chem
from rdkit.Chem import rdmolops
from PIL import Image
from rdkit.Chem import Draw, AllChem
from rdkit import DataStructs
from rdkit.Chem.Draw import MolDraw2DCairo
from copy import deepcopy
from src.fragment import fragment_and_draw, sanitize_filename
from src.extend_frag import fragment_data  
from src.extend import get_sparse_fingerprint, parse_fingerprint, calculate_tanimoto

from src import global_counters


def _clean_replacement_symbols(replaced_atoms):
    cleaned = []
    seen = set()

    for atom in replaced_atoms:
        if atom is None:
            continue

        symbol = str(atom).strip()
        if not symbol or symbol in {"*", "[*]"}:
            continue

        if symbol not in seen:
            cleaned.append(symbol)
            seen.add(symbol)

    return cleaned


def _get_original_symbols_for_mcs_atom(molecule_list, mcs_mol, mcs_atom_idx):
    symbols = []
    seen = set()

    for mol in molecule_list:
        if mol is None:
            continue

        match = mol.GetSubstructMatch(mcs_mol)
        if not match or mcs_atom_idx >= len(match):
            continue

        original_atom_idx = match[mcs_atom_idx]
        atom = mol.GetAtomWithIdx(original_atom_idx)
        symbol = atom.GetSymbol()

        if symbol and symbol not in {"*", "[*]"} and symbol not in seen:
            symbols.append(symbol)
            seen.add(symbol)

    if not symbols and mcs_atom_idx < mcs_mol.GetNumAtoms():
        atom = mcs_mol.GetAtomWithIdx(mcs_atom_idx)
        symbol = atom.GetSymbol()
        if symbol and symbol not in {"*", "[*]"}:
            symbols.append(symbol)

    return symbols


def _has_non_carbon_replacement(symbols):
    return any(symbol not in {"C", "H"} for symbol in symbols)


def mark_mcs_with_r(molecule_list, mcs_smiles, output_folder, differences):
    mcs_mol = Chem.MolFromSmarts(mcs_smiles)
    if not mcs_mol:
        return None, None, None, None, None

    original_mcs_atom_count = mcs_mol.GetNumAtoms()

    r_group_counts = {idx: 0 for idx in range(mcs_mol.GetNumAtoms())}
    mcs_with_r = Chem.RWMol(mcs_mol)
    r_group_mapping = {}
    x_group_mapping = {}
    z_group_mapping = {}
    bond_indices = []
    atom_indices = set()
    non_carbon_atom_replacements = {}
    z_group_replacements = {}

    for mol in molecule_list:
        match = mol.GetSubstructMatch(mcs_mol)
        r_group_temp_counts = {idx: 0 for idx in range(mcs_mol.GetNumAtoms())}

        for atom_idx in range(mol.GetNumAtoms()):
            if atom_idx not in match:
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx in match:
                        r_atom_idx = match.index(neighbor_idx)
                        r_group_temp_counts[r_atom_idx] += 1

        for idx in r_group_temp_counts:
            r_group_counts[idx] = max(r_group_counts[idx], r_group_temp_counts[idx])

    for r_atom_idx, count in r_group_counts.items():
        for _ in range(count):
            if mcs_with_r.GetAtomWithIdx(r_atom_idx).GetDegree() == 1:
                if r_atom_idx not in r_group_mapping:
                    r_group_mapping[r_atom_idx] = f"R{str(global_counters.global_r_group_counter)}"
                    global_counters.global_r_group_counter += 1

                r_group_atom = mcs_with_r.GetAtomWithIdx(r_atom_idx)
                r_group_atom.SetProp("atomLabel", r_group_mapping[r_atom_idx])

                for neighbor in r_group_atom.GetNeighbors():
                    atom_indices.add(neighbor.GetIdx())
            else:
                r_group_atom = Chem.Atom(0)
                r_group_idx = mcs_with_r.AddAtom(r_group_atom)
                mcs_with_r.AddBond(r_atom_idx, r_group_idx, Chem.BondType.SINGLE)

                r_group_mapping[r_group_idx] = f"R{str(global_counters.global_r_group_counter)}"
                mcs_with_r.GetAtomWithIdx(r_group_idx).SetProp("atomLabel", r_group_mapping[r_group_idx])

                atom_indices.add(r_atom_idx)
                global_counters.global_r_group_counter += 1

    rings = rdmolops.GetSymmSSSR(mcs_mol)

    for ring in rings:
        for atom_idx in ring:
            atom_idx = int(atom_idx)

            if atom_idx >= original_mcs_atom_count:
                continue

            if atom_idx in r_group_mapping:
                continue

            if atom_idx in x_group_mapping:
                continue

            replacement_symbols = _get_original_symbols_for_mcs_atom(
                molecule_list,
                mcs_mol,
                atom_idx
            )

            if not replacement_symbols:
                continue

            if not _has_non_carbon_replacement(replacement_symbols):
                continue

            source_atom = mcs_mol.GetAtomWithIdx(atom_idx)
            source_symbol = source_atom.GetSymbol()

            if source_symbol in {"C", "H"} and set(replacement_symbols).issubset({"C", "H"}):
                continue

            x_group_mapping[atom_idx] = f"X{str(global_counters.global_x_group_counter)}"
            global_counters.global_x_group_counter += 1

            mcs_with_r.ReplaceAtom(atom_idx, Chem.Atom(0))
            x_group_atom = mcs_with_r.GetAtomWithIdx(atom_idx)
            x_group_atom.SetProp("atomLabel", x_group_mapping[atom_idx])

            non_carbon_atom_replacements[x_group_mapping[atom_idx]] = replacement_symbols

    for atom_idx in range(original_mcs_atom_count):
        if atom_idx in r_group_mapping:
            continue

        if atom_idx in x_group_mapping:
            continue

        atom = mcs_with_r.GetAtomWithIdx(atom_idx)

        if atom.HasProp("atomLabel"):
            label = atom.GetProp("atomLabel")
            if label.startswith(("R", "X")):
                continue

        if atom.GetDegree() != 1 or len(atom.GetNeighbors()) != 1:
            continue

        replacement_symbols = _get_original_symbols_for_mcs_atom(
            molecule_list,
            mcs_mol,
            atom_idx
        )

        if not replacement_symbols:
            continue

        if not _has_non_carbon_replacement(replacement_symbols):
            continue

        z_group_mapping[atom_idx] = f"Z{str(global_counters.global_z_group_counter)}"
        global_counters.global_z_group_counter += 1

        atom.SetProp("atomLabel", z_group_mapping[atom_idx])
        z_group_replacements[z_group_mapping[atom_idx]] = replacement_symbols

    for bond in mcs_with_r.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if (begin_idx in atom_indices or end_idx in atom_indices) and not (
            begin_idx in r_group_mapping or end_idx in r_group_mapping
        ):
            bond_indices.append((begin_idx, end_idx))

    atom_r_mapping = {}

    for atom_idx in atom_indices:
        connected_r_groups = []

        for neighbor in mcs_with_r.GetAtomWithIdx(atom_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in r_group_mapping:
                connected_r_groups.append(r_group_mapping[neighbor_idx])

        atom_r_mapping[atom_idx] = connected_r_groups

    original_atom_mappings = map_indices_to_original(
        molecule_list,
        mcs_mol,
        atom_r_mapping
    )

    original_bond_mappings = map_bonds_to_original(
        molecule_list,
        mcs_mol,
        bond_indices
    )

    for atom_idx in range(mcs_with_r.GetNumAtoms()):
        atom = mcs_with_r.GetAtomWithIdx(atom_idx)

        if atom.GetSymbol() == "N" and atom.GetDegree() == 2:
            single_bonds = sum(
                1 for bond in atom.GetBonds()
                if bond.GetBondType() == Chem.BondType.SINGLE
            )

            if single_bonds == 2:
                atom.SetProp("atomLabel", "NH")

    for x_group, replaced_atoms in non_carbon_atom_replacements.items():
        create_x_group_folder(output_folder, x_group, replaced_atoms)

    for z_group, replaced_atoms in z_group_replacements.items():
        create_z_group_folder(output_folder, z_group, replaced_atoms)

    return (
        mcs_with_r,
        list(atom_indices),
        bond_indices,
        original_atom_mappings,
        original_bond_mappings
    )


def create_z_group_folder(output_folder, z_group, replaced_atoms):
    replaced_atoms = _clean_replacement_symbols(replaced_atoms)

    if not replaced_atoms:
        return

    z_group_folder = os.path.join(output_folder, z_group)
    os.makedirs(z_group_folder, exist_ok=True)

    z_group_file = os.path.join(z_group_folder, "replacements.txt")

    with open(z_group_file, "w") as f:
        for atom in replaced_atoms:
            f.write(f"{atom}\n")


def create_x_group_folder(output_folder, x_group, replaced_atoms):
    replaced_atoms = _clean_replacement_symbols(replaced_atoms)

    if not replaced_atoms:
        return

    x_group_folder = os.path.join(output_folder, x_group)
    os.makedirs(x_group_folder, exist_ok=True)

    x_group_file = os.path.join(x_group_folder, "replacements.txt")

    with open(x_group_file, "w") as f:
        for atom in replaced_atoms:
            f.write(f"{atom}\n")
            

def map_indices_to_original(molecule_list, mcs_mol, atom_r_mapping):
    original_mappings = []
    for mol in molecule_list:
        match = mol.GetSubstructMatch(mcs_mol)
        mol_mapping = {}
        for mcs_atom_idx, r_groups in atom_r_mapping.items():
            if mcs_atom_idx < len(match):
                original_atom_idx = match[mcs_atom_idx]
                mol_mapping[original_atom_idx] = r_groups
        original_mappings.append(mol_mapping)
    return original_mappings


def map_bonds_to_original(molecule_list, mcs_mol, bond_indices):
    original_bond_mappings = []
    for mol in molecule_list:
        match = mol.GetSubstructMatch(mcs_mol)
        bond_mapping = []
        if match: 
            for begin_idx, end_idx in bond_indices:
                if begin_idx < len(match) and end_idx < len(match):
                    original_begin_idx = match[begin_idx]
                    original_end_idx = match[end_idx]
                    original_bond = mol.GetBondBetweenAtoms(original_begin_idx, original_end_idx)
                    if original_bond:
                        original_bond_idx = original_bond.GetIdx()
                        bond_mapping.append((original_begin_idx, original_end_idx, original_bond_idx))
        original_bond_mappings.append(bond_mapping)
    return original_bond_mappings


def mark_nest_mcs_with_r(mcs_mapping, output_folder):
    initial_r_group_counter = global_counters.global_r_group_counter
    max_r_in_mcs = initial_r_group_counter
    marked_r_groups = []
    valid_molecules = []
    for i, (mcs_smarts, mcs_data) in enumerate(mcs_mapping.items(), start=1):
        global_counters.global_r_group_counter = initial_r_group_counter
        mcs_smiles = mcs_data['mcs_structure']
        fragments = mcs_data['fragments']
        mcs_mol = Chem.MolFromSmiles(mcs_smiles)

        if not mcs_mol:
            continue

        molecule_index = []  
        molecule_list = []
        for frag in fragments:
            frag_smiles = frag['smiles']  
            frag_file_name = frag['file_name'] 
            frag_mol = Chem.MolFromSmiles(frag_smiles)
            if frag_mol:
                molecule_list.append(frag_mol)
                fragment_index = int(frag_file_name.split('_')[1])  
                molecule_index.append((fragment_index, frag_mol)) 

        r_group_counts = {idx: 0 for idx in range(mcs_mol.GetNumAtoms())}
        r_group_mapping = {}
        bond_indices = []
        atom_indices = set()

        for mol in molecule_list:
            match = mol.GetSubstructMatch(mcs_mol)
            r_group_temp_counts = {idx: 0 for idx in range(mcs_mol.GetNumAtoms())}
            for atom_idx in range(mol.GetNumAtoms()):
                if atom_idx not in match:
                    for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx in match:
                            r_atom_idx = match.index(neighbor_idx)
                            r_group_temp_counts[r_atom_idx] += 1
            for idx in r_group_temp_counts:
                r_group_counts[idx] = max(r_group_counts[idx], r_group_temp_counts[idx])

        mcs_with_r = Chem.RWMol(mcs_mol)
        atom_replacements = []
        for atom in mcs_with_r.GetAtoms():
            if atom.GetSymbol() == '*':
                atom_replacements.append(atom.GetIdx())
                atom.SetAtomicNum(1)

        try:
            smiles_before_kekulize = Chem.MolToSmiles(mcs_with_r, isomericSmiles=True)
            Chem.Kekulize(mcs_with_r, clearAromaticFlags=True)
        except Chem.rdchem.KekulizeException as e:
            print(f"Error during Kekulization of molecule: {e}")

        for idx in atom_replacements:
            mcs_with_r.GetAtomWithIdx(idx).SetAtomicNum(0)

        for r_atom_idx, count in r_group_counts.items():
            for _ in range(count):
                r_group_atom = mcs_with_r.GetAtomWithIdx(r_atom_idx)

                if r_group_atom.GetDegree() == 1:  
                    if r_atom_idx not in r_group_mapping:
                        r_group_mapping[r_atom_idx] = f'R{(str(global_counters.global_r_group_counter))}'
                        marked_r_groups.append(f'R{(str(global_counters.global_r_group_counter))}')
                        global_counters.global_r_group_counter += 1
                    r_group_atom.SetProp('atomLabel', r_group_mapping[r_atom_idx])
                    for neighbor in r_group_atom.GetNeighbors():
                        atom_indices.add(neighbor.GetIdx())
                else:
                    new_r_group_atom = Chem.Atom(0) 
                    new_r_group_idx = mcs_with_r.AddAtom(new_r_group_atom)
                    mcs_with_r.AddBond(r_atom_idx, new_r_group_idx, Chem.BondType.SINGLE)
                    r_group_mapping[new_r_group_idx] = f'R{(str(global_counters.global_r_group_counter))}'
                    marked_r_groups.append(f'R{(str(global_counters.global_r_group_counter))}')
                    mcs_with_r.GetAtomWithIdx(new_r_group_idx).SetProp('atomLabel', r_group_mapping[new_r_group_idx])
                    atom_indices.add(r_atom_idx)
                    global_counters.global_r_group_counter += 1
                    max_r_in_mcs = max(max_r_in_mcs, global_counters.global_r_group_counter)

                    atom_info = {
                        "Index": r_atom_idx,
                        "Symbol": r_group_atom.GetSymbol(),
                        "Degree": r_group_atom.GetDegree(),
                    }

                    if atom_info['Symbol'] in ['N', 'S', 'P', 'O']:
                        if r_group_atom.GetDegree() > 0:
                            new_atom = Chem.Atom(r_group_atom.GetAtomicNum())
                            mcs_with_r.ReplaceAtom(r_atom_idx, new_atom)

        for bond in mcs_with_r.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if (begin_idx in atom_indices or end_idx in atom_indices) and not (begin_idx in r_group_mapping or end_idx in r_group_mapping):
                bond_indices.append((begin_idx, end_idx))

        atom_r_mapping = {}
        for atom_idx in atom_indices:
            connected_r_groups = []
            for neighbor in mcs_with_r.GetAtomWithIdx(atom_idx).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in r_group_mapping:
                    connected_r_groups.append(r_group_mapping[neighbor_idx])
            atom_r_mapping[atom_idx] = connected_r_groups

        original_atom_mappings = map_indices_to_original(molecule_list, mcs_mol, atom_r_mapping)
        original_bond_mappings = map_bonds_to_original(molecule_list, mcs_mol, bond_indices)

        if 'R' in mcs_smarts:
            continue
        
        # Here is the fix for properly generating and saving images
        asterisk_indices = [atom.GetIdx() for atom in mcs_with_r.GetAtoms() if atom.GetSymbol() == '*']
        r_indices = set(r_group_mapping.keys()) 
        true_asterisk_indices = [idx for idx in asterisk_indices if idx not in r_indices]  

        processed_asterisk_removal = False  
        
        if len(true_asterisk_indices) == 1:
            valid_molecules.append(mcs_with_r)
            
        if len(true_asterisk_indices) >= 2:
            try:
                # Phase 1: Record the non-asterisk end atom information of deleted double bonds
                retained_atom_mapnums = {}  # {original index: original mapping number}
                double_bond_pairs = []      # Record the non-asterisk end indices of cut double bonds

                for asterisk_idx in true_asterisk_indices:
                    atom = mcs_with_r.GetAtomWithIdx(asterisk_idx)
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other_idx = bond.GetOtherAtomIdx(asterisk_idx)
                            if other_idx not in true_asterisk_indices:
                                other_atom = mcs_with_r.GetAtomWithIdx(other_idx)
                                original_mapnum = other_atom.GetAtomMapNum()
                                retained_atom_mapnums[other_idx] = original_mapnum
                                other_atom.SetAtomMapNum(other_idx)  # Temporary marking
                                double_bond_pairs.append(other_idx)

                # Phase 2: Perform the original deletion operation (unchanged)
                for idx in sorted(true_asterisk_indices, reverse=True):
                    bonds = list(mcs_with_r.GetAtomWithIdx(idx).GetBonds())
                    for bond in bonds:
                        mcs_with_r.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    mcs_with_r.RemoveAtom(idx)
                processed_asterisk_removal = True

                # Phase 3: Double bond repair logic
                if len(double_bond_pairs) >= 2:
                    retained_atoms = []
                    for atom in mcs_with_r.GetAtoms():
                        if atom.GetAtomMapNum() in double_bond_pairs:
                            original_idx = atom.GetAtomMapNum()
                            atom.SetAtomMapNum(retained_atom_mapnums[original_idx])  # Restore original mapping number
                            retained_atoms.append(atom.GetIdx())
                    
                    if len(retained_atoms) == 2:
                        a, b = retained_atoms
                        bond = mcs_with_r.GetBondBetweenAtoms(a, b)
                        if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                            backup_mol = Chem.Mol(mcs_with_r.ToBinary())
                            try:
                                bond.SetBondType(Chem.BondType.DOUBLE)
                                Chem.SanitizeMol(mcs_with_r)
                            except:
                                mcs_with_r = backup_mol

                for atom in mcs_with_r.GetAtoms():
                    atom.SetAtomMapNum(0)  # Reset mapping number to zero

            except Exception as e:
                processed_asterisk_removal = False
            
        AllChem.Compute2DCoords(mcs_with_r)
        drawer = Draw.MolDraw2DCairo(300, 300)
        options = drawer.drawOptions()
        options.useBWAtomPalette()
        options.colorAtoms = False
        options.highlightColor = None

        drawer.DrawMolecule(mcs_with_r)
        drawer.FinishDrawing()

        png_data = drawer.GetDrawingText()
        img = Image.open(io.BytesIO(png_data))
        img = img.convert("RGB")

        # Saving the image for the structure
        if processed_asterisk_removal:
            img_path = os.path.join(output_folder, f'0_mcs_with_r_{i}_{sanitize_filename(mcs_smarts)}.png')
        else:
            img_path = os.path.join(output_folder, f'mcs_with_r_{i}_{sanitize_filename(mcs_smarts)}.png')
        img.save(img_path)

        mcs_mapping[mcs_smarts]['r_group_mapping'] = r_group_mapping
        fragment_and_draw(molecule_list, molecule_index, original_atom_mappings, original_bond_mappings, output_folder)

        import re

        try:
            mol_for_smiles = Chem.Mol(mcs_with_r)

            for atom in mol_for_smiles.GetAtoms():
                if atom.HasProp('atomLabel'):
                    label = atom.GetProp('atomLabel')
                    if isinstance(label, str) and label.startswith('R'):
                        digits = ''.join(ch for ch in label if ch.isdigit())
                        if digits:
                            try:
                                atom.SetAtomMapNum(int(digits))
                            except Exception:
                                pass

            raw_smiles = Chem.MolToSmiles(mol_for_smiles, canonical=False)

            smiles_with_r = raw_smiles
            replaced = set()
            for atom in mol_for_smiles.GetAtoms():
                if atom.HasProp('atomLabel'):
                    label = atom.GetProp('atomLabel')
                    if isinstance(label, str) and label.startswith('R'):
                        num = ''.join(ch for ch in label if ch.isdigit())
                        if not num or num in replaced:
                            continue
                        pattern = re.compile(r'\[[^\]]*:' + re.escape(num) + r'\]')
                        smiles_with_r, nsubs = pattern.subn(label, smiles_with_r, count=1)
                        if nsubs == 0:
                            pattern_alt = re.compile(r'\[[^\]]*' + re.escape(num) + r'\]')
                            smiles_with_r, _ = pattern_alt.subn(label, smiles_with_r, count=1)
                        replaced.add(num)

            txt_path = os.path.join(output_folder, "mcs_with_r_smiles.txt")
            with open(txt_path, "a", encoding="utf-8") as f:
                f.write(f"{i}. {smiles_with_r}\n")

        except Exception as e:
            print(f"[Warning] Failed to write R-labeled SMILES for {mcs_smarts}: {e}")

        
    generated_smiles, original_smiles = extend_mcs_with_r(valid_molecules, output_folder)
    extend_mcs_with_ring(original_smiles, output_folder)
    
    global_counters.global_r_group_counter = max_r_in_mcs
    return marked_r_groups


def is_equivalent_to_existing(smiles, existing_smiles_set):
    mol1 = Chem.MolFromSmiles(smiles)
    if mol1 is None:
        return False
    smi1 = Chem.MolToSmiles(mol1, canonical=True)
    
    for ex_smiles in existing_smiles_set:
        mol2 = Chem.MolFromSmiles(ex_smiles)
        if mol2 is None:
            continue
        smi2 = Chem.MolToSmiles(mol2, canonical=True)
        if smi1 == smi2:
            return True
    return False
    
    
def generate_star_extensions(mol):
    extended_mols = []

    # count complex bicyclic systems
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    complex_ring_count = 0
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            shared = set(atom_rings[i]) & set(atom_rings[j])
            if len(shared) > 1:
                complex_ring_count += 1
                break

    allow_complex = complex_ring_count < 3

    # find star
    star_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == '*']
    if not star_atoms:
        return extended_mols

    star_idx = star_atoms[0]
    star_bonds = mol.GetAtomWithIdx(star_idx).GetBonds()
    if len(star_bonds) != 1:
        return extended_mols

    anchor_idx = star_bonds[0].GetOtherAtomIdx(star_idx)

    # find ring system containing anchor
    ring_atoms = set()
    for ring in atom_rings:
        if anchor_idx in ring:
            ring_atoms.update(ring)

    if not ring_atoms:
        return extended_mols

    for target_idx in ring_atoms:
        if target_idx == anchor_idx:
            continue

        target_atom = mol.GetAtomWithIdx(target_idx)

        if target_atom.GetSymbol() != 'C':
            continue

        if any(
            mol.GetAtomWithIdx(b.GetOtherAtomIdx(target_idx)).GetSymbol() in ['*', 'R']
            for b in target_atom.GetBonds()
        ):
            continue

        # if this is a complex bicyclic position and not allowed, skip
        if not allow_complex:
            shared_count = 0
            for ring in atom_rings:
                if target_idx in ring and anchor_idx in ring:
                    shared_count += 1
            if shared_count > 1:
                continue

        mol_copy = deepcopy(mol)

        mol_copy.RemoveBond(star_idx, anchor_idx)
        mol_copy.AddBond(star_idx, target_idx, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(mol_copy)
        except Exception:
            continue

        extended_mols.append((mol_copy, target_idx))

    return extended_mols


def apply_fragment_extension_with_r(mol, folder, num_input, r_numbers=None):
    if mol is None:
        return []

    if is_complex_bicyclic(mol):
        return []

    if num_input == 1:
        top_n = 6
    elif num_input == 2:
        top_n = 4
    elif num_input == 3:
        top_n = 3
    elif num_input == 4:
        top_n = 2
    else:
        top_n = 1

    if r_numbers is None:
        r_numbers = []

    input_star_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() == '*']
    if not input_star_atoms:
        return []

    input_star_positions = [a.GetIdx() for a in input_star_atoms]
    cleaned_smiles = Chem.MolToSmiles(mol)
    input_starts_with_star_eq = cleaned_smiles.startswith("*=")

    try:
        target_fp_str = get_sparse_fingerprint(cleaned_smiles)
        target_fp = parse_fingerprint(target_fp_str)
    except ValueError:
        return []

    def valid_fragment(frag_mol):
        if frag_mol is None:
            return False
        ri = frag_mol.GetRingInfo()
        if ri.NumRings() == 0:
            return False
        for ring in ri.AtomRings():
            for idx in ring:
                atom = frag_mol.GetAtomWithIdx(idx)
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetSymbol() != '*':
                        return False
        stars = [a for a in frag_mol.GetAtoms() if a.GetSymbol() == '*']
        if not stars:
            return False
        for a in stars:
            if a.IsInRing():
                return False
            connected_to_ring = False
            for nb in a.GetNeighbors():
                bond = frag_mol.GetBondBetweenAtoms(a.GetIdx(), nb.GetIdx())
                if nb.IsInRing():
                    if input_starts_with_star_eq and bond.GetBondType() == Chem.BondType.DOUBLE:
                        connected_to_ring = True
                    elif not input_starts_with_star_eq and bond.GetBondType() == Chem.BondType.SINGLE:
                        connected_to_ring = True
            if not connected_to_ring:
                return False
        return True

    candidates = []
    for frag_smiles, frag_fp_str in fragment_data.items():
        if not isinstance(frag_smiles, str):
            continue
        if "@" in frag_smiles or "/" in frag_smiles or "\\" in frag_smiles:
            continue
        brackets = re.findall(r"\[.*?\]", frag_smiles)
        if any(b not in ["[nH]", "[SH]"] for b in brackets):
            continue
        frag_starts_with_star_eq = frag_smiles.startswith("*=")
        if input_starts_with_star_eq != frag_starts_with_star_eq:
            continue

        try:
            ext_mol = Chem.MolFromSmiles(frag_smiles)
        except Exception:
            continue
        if ext_mol is None or not valid_fragment(ext_mol):
            continue
        try:
            frag_fp = parse_fingerprint(frag_fp_str)
        except Exception:
            continue
        sim = calculate_tanimoto(target_fp, frag_fp)
        candidates.append((sim, ext_mol, frag_smiles))

    if not candidates:
        print(f"[DEBUG] Input: {cleaned_smiles}, no valid fragment candidates.")
        return []

    candidates.sort(reverse=True, key=lambda x: x[0])
    top_candidates = candidates[:top_n]

    results = []
    for sim, frag_mol, frag_smi in top_candidates:
        frag_rw = Chem.RWMol(frag_mol)

        frag_star_atoms = [a for a in frag_rw.GetAtoms() if a.GetSymbol() == '*']
        frag_star_positions = [a.GetIdx() for a in frag_star_atoms]

        print(f"[DEBUG] Original * positions in fragment: {frag_star_positions}")

        num_to_add = len(input_star_positions) - len(frag_star_positions)
        if num_to_add > len(r_numbers):
            print(f"[DEBUG] Not enough R numbers for new *")
            continue

        used_r_numbers = r_numbers[:]

        for input_idx, input_star_idx in enumerate(input_star_positions):
            if input_idx < len(frag_star_positions):
                continue

            attach_idx = min(input_star_idx, frag_rw.GetNumAtoms() - 2)
            atom = frag_rw.GetAtomWithIdx(attach_idx)

            if atom.GetSymbol() == '*' or any(nb.GetSymbol() == '*' for nb in atom.GetNeighbors()):
                continue

            new_atom = Chem.Atom("*")
            new_idx = frag_rw.AddAtom(new_atom)
            frag_rw.AddBond(attach_idx, new_idx, Chem.BondType.SINGLE)

            r_num = used_r_numbers.pop(0)
            frag_rw.GetAtomWithIdx(new_idx).SetProp("atomLabel", f"R{r_num}")
            frag_rw.GetAtomWithIdx(new_idx).SetProp("added_star", "1")

            print(f"[DEBUG] Adding * at frag atom {attach_idx}, marked as R{r_num}")

        frag_final = frag_rw.GetMol()
        if frag_final is None:
            continue

        try:
            Chem.SanitizeMol(frag_final)
        except Exception:
            continue

        bad = False
        for atom in frag_final.GetAtoms():
            if atom.GetSymbol() == '*':
                if any(nb.GetSymbol() == '*' for nb in atom.GetNeighbors()):
                    bad = True
                    break
        if bad:
            continue

        ext_smiles_final = Chem.MolToSmiles(frag_final)
        print(f"[DEBUG] Input: {cleaned_smiles}, Fragment candidate: {frag_smi}, Selected: {ext_smiles_final}")
        results.append((frag_final, 0))

    return results


def extend_mcs_with_r(valid_molecules, output_folder):
    generated_smiles = set()
    original_smiles = set()

    for i, mcs_with_r in enumerate(valid_molecules):
        smiles = Chem.MolToSmiles(mcs_with_r)
        generated_smiles.add(smiles)
        original_smiles.add(smiles)

    for i, mcs_with_r in enumerate(valid_molecules):
        if mcs_with_r is None:
            continue

        mol = Chem.RWMol(mcs_with_r)
        if mol is None:
            continue

        r_numbers = []
        for atom in mol.GetAtoms():
            if atom.HasProp("atomLabel"):
                label = atom.GetProp("atomLabel")
                if label.startswith("R"):
                    num = ''.join(ch for ch in label if ch.isdigit())
                    if num:
                        r_numbers.append(int(num))
        r_numbers = sorted(list(set(r_numbers)))

        star_results = generate_star_extensions(mol)
        frag_results = apply_fragment_extension_with_r(mol, folder=output_folder, num_input=len(valid_molecules), r_numbers=r_numbers)
        extended_results = star_results + frag_results

        if not extended_results:
            continue

        for idx, (mol_copy, pos) in enumerate(extended_results):
            try:
                Chem.SanitizeMol(mol_copy)
            except Exception:
                continue

            extend_smiles = Chem.MolToSmiles(mol_copy)
            if extend_smiles in generated_smiles:
                continue
            if is_equivalent_to_existing(extend_smiles, generated_smiles):
                continue
            generated_smiles.add(extend_smiles)

            AllChem.Compute2DCoords(mol_copy)
            drawer = Draw.MolDraw2DCairo(300, 300)
            options = drawer.drawOptions()
            options.useBWAtomPalette()
            options.colorAtoms = False
            options.highlightColor = None
            drawer.DrawMolecule(mol_copy)
            drawer.FinishDrawing()

            png_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(png_data)).convert("RGB")
            img_path = os.path.join(output_folder, f'mcs_with_r_{i + 1}_extended_{idx}_{pos}.png')
            img.save(img_path)

            try:
                mol_for_smiles = Chem.Mol(mol_copy)
                raw_smiles = Chem.MolToSmiles(mol_for_smiles)
                smiles_with_r = raw_smiles
                replaced = set()
                for atom in mol_for_smiles.GetAtoms():
                    if atom.HasProp("atomLabel"):
                        label = atom.GetProp("atomLabel")
                        if label.startswith("R"):
                            num = ''.join(ch for ch in label if ch.isdigit())
                            if not num or num in replaced:
                                continue
                            pattern = re.compile(r'\[[^\]]*:' + re.escape(num) + r'\]')
                            smiles_with_r, nsubs = pattern.subn(label, smiles_with_r, count=1)
                            if nsubs == 0:
                                pattern_alt = re.compile(r'\[[^\]]*' + re.escape(num) + r'\]')
                                smiles_with_r, _ = pattern_alt.subn(label, smiles_with_r, count=1)
                            replaced.add(num)

                txt_path = os.path.join(output_folder, "mcs_with_r_smiles.txt")
                with open(txt_path, "a", encoding="utf-8") as f:
                    f.write(f"{i + 1}_extended_{idx}_{pos}. {smiles_with_r}\n")

            except Exception:
                pass

    return generated_smiles, original_smiles


def is_benzene_ring(mol, ring_atoms):
    """Determine if the ring is a benzene ring (aromatic six-membered carbon ring)"""
    if len(ring_atoms) != 6:
        return False
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'C' or not atom.GetIsAromatic():
            return False
    bonds = mol.GetBonds()
    for bond in bonds:
        if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
            if bond.GetBondType() != Chem.BondType.AROMATIC:
                return False
    return True


def is_cyclohexane_ring(mol, ring_atoms):
    """Determine if the ring is a cyclohexane (non-aromatic six-membered single-bond carbon ring)"""
    if len(ring_atoms) != 6:
        return False
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'C':
            return False
    bonds = mol.GetBonds()
    for bond in bonds:
        if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False
    return True


def adjust_nitrogen_hydrogens(atom):
    """Adjust the number of hydrogens on a nitrogen atom to match its valence"""
    if atom.GetSymbol() != 'N':
        return
    current_h = atom.GetNumExplicitHs()
    if current_h > 0:
        atom.SetNumExplicitHs(current_h - 1)
    else:
        implicit_h = atom.GetNumImplicitHs()
        if implicit_h > 0:
            atom.SetNumExplicitHs(implicit_h - 1)


def extend_mcs_with_ring(generated_smiles, output_folder):
    if isinstance(generated_smiles, set):
        generated_smiles = " ".join(generated_smiles)
    smiles_list = generated_smiles.split()
    
    for i, smile in enumerate(smiles_list):
        if smile.count('*') != 1:
            continue
        
        mol = Chem.MolFromSmiles(smile)
        if not mol:
            continue
        
        # if is_complex_bicyclic(mol):
        #     # print(f"Skipping complex bicyclic structure: {smile}")
        #     continue  
              
        try:
            Chem.Kekulize(mol)
        except Chem.KekulizeException:
            print("Kekulization failed, trying to continue")
        
        Chem.SanitizeMol(mol)
        
        ring_info = list(rdmolops.GetSymmSSSR(mol))
        reactive_atoms = set()
        
        # Step 1: Identify reactive atoms
        for ring in ring_info:
            ring_atoms = list(ring)
            if is_benzene_ring(mol, ring_atoms) or is_cyclohexane_ring(mol, ring_atoms):
                # print(f"Skipping excluded ring: {ring_atoms}")
                continue
            
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                atom.UpdatePropertyCache()
                
                if atom.GetSymbol() not in ['C', 'N']:
                    continue
                
                used_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
                
                if atom.GetSymbol() == 'C':
                    remaining_valence = 4 - used_valence
                else:
                    explicit_degree = atom.GetDegree()
                    remaining_valence = 3 - explicit_degree
                    if remaining_valence < 1 and atom.GetNumImplicitHs() > 0:
                        remaining_valence += 1
                
                if remaining_valence >= 1:
                    reactive_atoms.add(atom_idx)
        
        # Step 2: Find adjacent pairs
        adj_pairs = []
        for atom_idx in reactive_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in reactive_atoms and neighbor_idx > atom_idx:
                    for ring in ring_info:
                        if atom_idx in ring and neighbor_idx in ring:
                            ring_size = len(ring)
                            adj_pairs.append( (atom_idx, neighbor_idx, ring_size) )
                            break
        
        if not adj_pairs:
            # print(f"No adjacent reactive pairs in {smile}")
            continue
        
        a1, a2, ring_size = adj_pairs[0]
        
        emol = Chem.EditableMol(mol)
        new_atoms = []
        
        if ring_size == 5:
            # Add a six-membered ring
            for _ in range(4):
                new_atoms.append(emol.AddAtom(Chem.Atom('C')))
            emol.AddBond(a1, new_atoms[0], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[0], new_atoms[1], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[1], new_atoms[2], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[2], new_atoms[3], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[3], a2, Chem.BondType.SINGLE)
        elif ring_size == 6:
            # Add a five-membered ring
            for _ in range(3):
                new_atoms.append(emol.AddAtom(Chem.Atom('C')))
            emol.AddBond(a1, new_atoms[0], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[0], new_atoms[1], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[1], new_atoms[2], Chem.BondType.SINGLE)
            emol.AddBond(new_atoms[2], a2, Chem.BondType.SINGLE)
        else:
            continue
        
        new_mol = emol.GetMol()
        
        # Check again after extension: did we create a complex fused/bridged ring?
        if is_complex_bicyclic(new_mol):
            # print(f"Created bridged ring after expansion, skipping: {smile}")
            continue   
             
        # Adjust hydrogens on nitrogen atoms
        for atom_idx in [a1, a2]:
            atom = new_mol.GetAtomWithIdx(atom_idx)
            adjust_nitrogen_hydrogens(atom)
        
        try:
            Chem.SanitizeMol(new_mol)
            AllChem.Compute2DCoords(new_mol)
        except Exception as e:
            # print(f"Sanitization error: {e}")
            continue

        # Create an image of the molecule using your provided method
        drawer = Draw.MolDraw2DCairo(300, 300)
        options = drawer.drawOptions()
        options.useBWAtomPalette()  # Use black and white atom palette
        options.colorAtoms = False  # Do not color atoms
        options.highlightColor = None  # No highlight color

        drawer.DrawMolecule(new_mol)
        drawer.FinishDrawing()

        # Save the image to the specified output folder
        png_data = drawer.GetDrawingText()
        img = Image.open(io.BytesIO(png_data)).convert("RGB")
        img_path = os.path.join(output_folder, f'mcs_with_r_ring_{i + 1}_extended_{a1}_{a2}.png')  # i starts from 1 
        img.save(img_path)


        try:
            mol_for_smiles = Chem.Mol(new_mol)

            for atom in mol_for_smiles.GetAtoms():
                if atom.HasProp("atomLabel"):
                    label = atom.GetProp("atomLabel")
                    if label.startswith("R"):
                        num = ''.join(ch for ch in label if ch.isdigit())
                        if num:
                            try:
                                atom.SetAtomMapNum(int(num))
                            except Exception:
                                pass

            raw_smiles = Chem.MolToSmiles(mol_for_smiles)

            smiles_with_r = raw_smiles
            replaced = set()
            for atom in mol_for_smiles.GetAtoms():
                if atom.HasProp("atomLabel"):
                    label = atom.GetProp("atomLabel")
                    if label.startswith("R"):
                        num = ''.join(ch for ch in label if ch.isdigit())
                        if not num or num in replaced:
                            continue
                        pattern = re.compile(r'\[[^\]]*:' + re.escape(num) + r'\]')
                        smiles_with_r, nsubs = pattern.subn(label, smiles_with_r, count=1)
                        if nsubs == 0:
                            pattern_alt = re.compile(r'\[[^\]]*' + re.escape(num) + r'\]')
                            smiles_with_r, _ = pattern_alt.subn(label, smiles_with_r, count=1)
                        replaced.add(num)

            txt_path = os.path.join(output_folder, "mcs_with_r_smiles.txt")
            with open(txt_path, "a", encoding="utf-8") as f:
                f.write(f"{i + 1}_ring_{a1}_{a2}. {smiles_with_r}\n")

        except Exception as e:
            print(f"[Warning] Failed to write ring-extended SMILES for molecule {i+1}_ring_{a1}_{a2}: {e}")
        # ===========================================================


    return


def is_complex_bicyclic(mol):
    """
    Check if the molecule contains overlapping rings (bridged or fused),
    suggesting it's a complex non-planar ring system.
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Check for atom overlap across multiple rings (shared atoms)
    for i in range(len(atom_rings)):
        for j in range(i + 1, len(atom_rings)):
            shared = set(atom_rings[i]) & set(atom_rings[j])
            if len(shared) > 1:  # Shared more than 1 atom = fused or bridged ring
                return True
    return False

