# models/molecule.py
from extensions import db  # Import db from extensions
from datetime import datetime

class Molecule(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False)
    smiles = db.Column(db.String(500), nullable=False)
    formula = db.Column(db.String(100))
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    file_path = db.Column(db.String(200))
    
    def __repr__(self):
        return f'<Molecule {self.name}: {self.smiles}>'