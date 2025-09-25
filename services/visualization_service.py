import numpy as np
from rdkit import Chem

class VisualizationService:
    def calculate_bond_angles(self, mol):
        """Calculate bond angles of the molecule"""
        bond_angles = []
        conf = mol.GetConformer()
        
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()
            
            neighbors1 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(atom1).GetNeighbors() if atom.GetIdx() != atom2]
            neighbors2 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(atom2).GetNeighbors() if atom.GetIdx() != atom1]
            
            for n1 in neighbors1:
                for n2 in neighbors2:
                    pos1 = np.array(conf.GetAtomPosition(n1))
                    pos2 = np.array(conf.GetAtomPosition(atom1))
                    pos3 = np.array(conf.GetAtomPosition(atom2))
                    pos4 = np.array(conf.GetAtomPosition(n2))
                    
                    vec1 = pos1 - pos2
                    vec2 = pos3 - pos2
                    angle = np.degrees(np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))))
                    bond_angles.append((n1, atom1, atom2, angle))
                    
        return bond_angles
        
    def get_atom_colors(self):
        """Return atom color scheme"""
        return {
            "C": "gray",
            "H": "white",
            "O": "red",
            "N": "blue",
            "Cl": "green",
            "Br": "brown",
            "F": "cyan",
            "P": "orange",
            "S": "yellow",
        }
        
    def get_bond_colors(self):
        """Return bond color scheme"""
        return {
            Chem.rdchem.BondType.SINGLE: "black",
            Chem.rdchem.BondType.DOUBLE: "blue",
            Chem.rdchem.BondType.TRIPLE: "purple",
        }