from src.extend_frag import fragment_data  
from rdkit import Chem
from rdkit.Chem import AllChem
import re
from src.description import description_data 



def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        standardized_smiles = Chem.MolToSmiles(mol, canonical=True)
        return standardized_smiles
    return None


def clean_smiles(smiles):
    bond_marker = re.findall(r'\[\d*\*\]', smiles)
    
    if len(bond_marker) == 1:
        cleaned_smiles = re.sub(r'\[\d*\*\]', '*', smiles)
        return cleaned_smiles
    
    elif len(bond_marker) >= 2:
        cleaned_smiles = re.sub(r'\[\d*\*\]', '[*]', smiles)
        cleaned_smiles = re.sub(r'\(\[\*\]\)', '', cleaned_smiles)  
        cleaned_smiles = re.sub(r'\=\[\*\]', '', cleaned_smiles)
        cleaned_smiles = re.sub(r'\#\[\*\]', '', cleaned_smiles)
        cleaned_smiles = re.sub(r'\[\*\]\#', '', cleaned_smiles)
        cleaned_smiles = re.sub(r'\[\*\]\=', '', cleaned_smiles)
        cleaned_smiles = re.sub(r'\[\*\]', '', cleaned_smiles)
        cleaned_smiles = '*' + cleaned_smiles
        return cleaned_smiles

    return smiles


def get_sparse_fingerprint(smiles, radius=2):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    fp = AllChem.GetMorganFingerprint(mol, radius)
    fp_dict = fp.GetNonzeroElements()
    fp_str = ",".join(f"{bit}:{count}" for bit, count in fp_dict.items())
    return fp_str


def parse_fingerprint(fp_str):
    fp_dict = {}
    if not fp_str:
        return fp_dict
    for pair in fp_str.split(','):
        if not pair:
            continue
        bit, count = pair.split(':')
        fp_dict[int(bit)] = int(count)
    return fp_dict


def calculate_tanimoto(target_fp, query_fp):
    common_bits = set(target_fp.keys()) & set(query_fp.keys())
    intersection = sum(min(target_fp[bit], query_fp[bit]) for bit in common_bits)
    sum_target = sum(target_fp.values())
    sum_query = sum(query_fp.values())
    denominator = sum_target + sum_query - intersection
    
    if denominator == 0:
        return 0.0
    return intersection / denominator


def fragment_extension(cleaned_smiles, folder):
    if cleaned_smiles == "[H]":
        return []

    try:
        target_fp_str = get_sparse_fingerprint(cleaned_smiles)
        target_fp = parse_fingerprint(target_fp_str)
    except ValueError as e:
        print(f"Error processing SMILES: {e}")
        return []

    similarities = []
    for frag_smiles, frag_fp_str in fragment_data.items():
        if "@" in frag_smiles or "/" in frag_smiles or "\\" in frag_smiles or "[CH]" in frag_smiles:
            continue
        frag_fp = parse_fingerprint(frag_fp_str)
        similarity = calculate_tanimoto(target_fp, frag_fp)
        similarities.append((similarity, frag_smiles))

    similarities.sort(reverse=True, key=lambda x: x[0])
    top_results = similarities[:20]

    results = []
    matched_fragments = set(smiles for _, smiles in top_results)

    def find_description(smiles):
        for entry in description_data:
            smiles_list = []
            for sm in entry.get("SMILE", "").split(";"):
                cleaned_smile = standardize_smiles(clean_smiles(sm))
                if cleaned_smile is not None:
                    smiles_list.append(cleaned_smile)
                else:
                    smiles_list.append(sm)

            if smiles in smiles_list:
                return entry.get("Description")
        return None

    for _, smiles in top_results:
        description = find_description(smiles)
        if description:
            results.append(f"Molecule_0_extend.smiles: {description}")

    special_cases = {"*C", "*CC", "*C(C)C", "*CCC"}
    if cleaned_smiles in special_cases and "*C1CC1" not in matched_fragments:
        description = find_description("*C1CC1")
        if description:
            results.append(f"Molecule_0_extend.smiles: {description}")

    return results


def mol_extension(cleaned_smiles, folder):
    if cleaned_smiles == "[H]":
        return []

    extension_map = {
        "*C#N": "*C(=O)OC(C)C",
        "*C(=O)OC": "*C#N",
        "*N#C": "*OC(=O)C",
    }

    smiles_to_process = [(cleaned_smiles, "original")]
    if cleaned_smiles in extension_map:
        smiles_to_process.append((extension_map[cleaned_smiles], "extended"))

    results = []
    seen_smiles = set()

    # ---------- Part 1: fingerprint-based extension ----------
    for query_smiles, tag in smiles_to_process:
        try:
            target_fp_str = get_sparse_fingerprint(query_smiles)
            target_fp = parse_fingerprint(target_fp_str)
        except ValueError:
            continue

        similarities = []
        for frag_smiles, frag_fp_str in fragment_data.items():
            frag_fp = parse_fingerprint(frag_fp_str)
            similarity = calculate_tanimoto(target_fp, frag_fp)
            similarities.append((similarity, frag_smiles))

        similarities.sort(reverse=True, key=lambda x: x[0])
        top_results = similarities[:20]

        for _, smiles in top_results:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            if mol.GetRingInfo().AtomRings():
                continue

            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            if canonical_smiles not in seen_smiles:
                seen_smiles.add(canonical_smiles)
                results.append(canonical_smiles)

    # ---------- Part 2: heuristic SMILES expansion ----------
    heuristic_candidates = set()

    if cleaned_smiles.startswith("*"):
        base = cleaned_smiles[1:]

        # basic chain growth
        heuristic_candidates.add("*C" + base)
        heuristic_candidates.add("*CC" + base)

        # chain shrink
        if base.startswith("CC"):
            heuristic_candidates.add("*C" + base[2:])

        # carbonyl substitution
        for i in range(len(base)):
            if base[i] == "C":
                heuristic_candidates.add("*" + base[:i] + "C(=O)" + base[i+1:])

        # C -> CC expansion
        for i in range(len(base)):
            if base[i] == "C":
                heuristic_candidates.add("*" + base[:i] + "CC" + base[i+1:])

    # validate heuristic candidates with RDKit
    for smi in heuristic_candidates:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue

        invalid = False

        for atom in mol.GetAtoms():
            smarts = atom.GetSmarts()

            if smarts.startswith("[") and atom.GetSymbol() == "C":
                invalid = True
                break

            if atom.GetNumRadicalElectrons() != 0:
                invalid = True
                break

            symbol = atom.GetSymbol()
            valence = atom.GetTotalValence()

            if symbol == "C" and valence > 4:
                invalid = True
                break
            if symbol == "O" and valence > 2:
                invalid = True
                break
            if symbol == "N" and valence > 3:
                invalid = True
                break

        if invalid:
            continue

        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "*":
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() in ("O", "N"):
                        invalid = True
                        break
            if invalid:
                break

        if invalid:
            continue

        canonical = Chem.MolToSmiles(mol, canonical=True)
        if canonical not in seen_smiles:
            seen_smiles.add(canonical)
            results.append(canonical)

    return results
