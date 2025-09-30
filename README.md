# 🧬 ChemViz - Interactive 3D Molecular Visualizer

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![Flask](https://img.shields.io/badge/Flask-2.3.3-green.svg)](https://flask.palletsprojects.com)
[![RDKit](https://img.shields.io/badge/RDKit-2023.3.3-orange.svg)](https://rdkit.org)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Transform chemical structures into stunning 3D visualizations with ChemViz - a modern web application that brings molecules to life through interactive 3D rendering, powered by RDKit and 3Dmol.js.

![ChemViz Demo](https://via.placeholder.com/800x400/1a1a1a/7c3aed?text=ChemViz+3D+Molecular+Visualization)

## ✨ Features

### 🎯 Core Functionality
- **SMILES to 3D**: Convert chemical SMILES strings into interactive 3D molecular models
- **Real-time Visualization**: Instant 3D rendering with smooth animations and controls
- **Multiple Display Modes**: Switch between ball-and-stick and space-filling representations
- **Interactive Controls**: Rotate, zoom, spin, and manipulate molecules in real-time

### 🎨 Modern UI/UX
- **Dark Theme**: Professional dark interface optimized for molecular visualization
- **Responsive Design**: Works seamlessly on desktop, tablet, and mobile devices
- **Intuitive Controls**: User-friendly sidebar with organized control panels
- **Visual Feedback**: Real-time property display and molecular information

### 📊 Chemical Intelligence
- **Molecular Properties**: Automatic calculation of formula, molar mass, and other properties
- **Atom Legend**: Color-coded atomic visualization with standard CPK colors
- **Bond Analysis**: Advanced bond angle calculations and structural insights
- **Data Persistence**: SQLite database for storing and managing molecular data

## 🚀 Why Choose ChemViz?

### 🔬 **For Researchers & Educators**
- Visualize complex molecular structures instantly
- Create engaging educational content
- Analyze chemical properties at a glance
- Export high-quality molecular images

### 💻 **Modern Tech Stack**
- **Backend**: Flask (Python) - Robust and scalable web framework
- **Database**: SQLAlchemy with SQLite - Lightweight yet powerful data storage
- **Chemistry**: RDKit - Industry-standard computational chemistry toolkit
- **3D Rendering**: 3Dmol.js - High-performance WebGL molecular visualization
- **Frontend**: Modern HTML5/CSS3/JavaScript with responsive design

### 🌟 **User Experience**
- **Zero Installation**: Web-based interface accessible from any browser
- **Fast Performance**: Optimized rendering and caching for smooth interactions
- **Professional Interface**: Clean, modern design focused on usability
- **Keyboard Shortcuts**: Power-user features for efficient navigation

## 🛠️ Installation & Setup

### Prerequisites
- **Python 3.8+** (Python 3.10+ recommended)
- **Conda** (Anaconda/Miniconda) - *Why Conda? See below!*

### 🐍 Why Use Conda?

Conda is **essential** for this project because:

1. **RDKit Dependency Management**: RDKit has complex C++ dependencies that are notoriously difficult to install with pip alone
2. **Scientific Computing**: Conda excels at managing scientific Python packages with native dependencies
3. **Environment Isolation**: Prevents conflicts between different Python projects
4. **Cross-Platform Compatibility**: Works consistently across Windows, macOS, and Linux
5. **Binary Package Distribution**: Pre-compiled packages mean faster, more reliable installations

### Step-by-Step Installation

#### 1️⃣ **Clone the Repository**
```bash
git clone https://github.com/rhr18818/chem_3d_model.git
cd chem_3d_model
```

#### 2️⃣ **Create Conda Environment**
```bash
# Create a new conda environment with Python 3.10
conda create -n chemviz python=3.10 -y

# Activate the environment
conda activate chemviz
```

#### 3️⃣ **Install RDKit via Conda**
```bash
# Install RDKit from conda-forge (this is why we need conda!)
conda install -c conda-forge rdkit -y
```

#### 4️⃣ **Install Python Dependencies**
```bash
# Install remaining packages via pip
pip install -r requirements.txt
```

#### 5️⃣ **Initialize Database**
```bash
# Create the database tables
flask init-db
```

#### 6️⃣ **Launch the Application**
```bash
# Start the development server
python app.py
```

🎉 **Success!** Open your browser and navigate to `http://localhost:5001`

## 📖 Usage Guide

### 🔤 **Input Methods**

#### **SMILES Strings** (Recommended)
SMILES (Simplified Molecular Input Line Entry System) is a chemical notation that describes molecular structure:

```
Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Benzene: c1ccccc1
Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
Morphine: CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O
```

#### **Common Examples to Try**
- `CCO` - Ethanol (alcohol)
- `C6H12O6` - Won't work! Use SMILES: `C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O` (Glucose)
- `c1ccccc1` - Benzene (aromatic ring)
- `CC(C)C` - Isobutane (branched alkane)

### 🎮 **Interactive Controls**

#### **3D Viewer Controls**
- **Mouse Drag**: Rotate molecule in 3D space
- **Mouse Wheel**: Zoom in/out
- **Double Click**: Center and fit molecule to view

#### **Display Modes**
- **Ball-and-Stick**: Classic representation showing atoms as spheres and bonds as cylinders
- **Space-Filling**: Van der Waals representation showing molecular volume

#### **Action Buttons**
- **🔄 Reset**: Return to default view and zoom level
- **▶️ Spin**: Toggle automatic rotation animation
- **💾 Save**: Save molecule to your collection (future feature)
- **📥 Download**: Export current view as PNG image

#### **Keyboard Shortcuts**
- `Ctrl + K`: Focus search box
- `Ctrl + 1`: Toggle left sidebar (future feature)
- `Ctrl + 2`: Toggle right sidebar (future feature)

### 📊 **Understanding the Interface**

#### **Left Sidebar - Controls**
- Search functionality
- Display mode toggles
- Interactive controls
- Action buttons

#### **Center Panel - 3D Viewer**
- Main molecular visualization area
- Real-time 3D rendering
- Interactive manipulation space

#### **Right Sidebar - Information**
- Atom legend with CPK colors
- Molecular properties
- Chemical data display

## 🗂️ Project Structure

```
chem_3d_model/
├── 📁 app.py                 # Main Flask application
├── 📁 config.py             # Configuration settings
├── 📁 extensions.py         # Database extensions
├── 📁 requirements.txt      # Python dependencies
├── 📁 data/
│   └── molecules/           # Stored molecule files (.pkl)
├── 📁 instance/
│   └── molecules.db         # SQLite database
├── 📁 models/
│   └── molecule.py          # Database models
├── 📁 services/
│   ├── molecule_service.py  # Business logic for molecules
│   └── visualization_service.py  # 3D rendering services
├── 📁 static/
│   ├── css/                 # Stylesheets
│   ├── js/                  # JavaScript files
│   └── images/              # Static images
├── 📁 templates/
│   ├── base.html           # Base template
│   ├── index.html          # Home page
│   ├── gallery.html        # Molecule gallery
│   └── visualize.html      # 3D visualization page
└── 📁 tests/
    └── test_app.py         # Unit tests
```

## 🧪 Technical Details

### **Core Workflow**
```python
# Your core workflow:
SMILES → RDKit (Python) → 3D Coordinates → Database → 3Dmol.js
# This is naturally Python-centric
```

The application follows a streamlined data pipeline:
1. **Input**: User provides SMILES string through web interface
2. **Processing**: RDKit converts SMILES to 3D molecular structure
3. **Storage**: Molecule data and coordinates saved to SQLite database
4. **Visualization**: 3Dmol.js renders interactive 3D representation
5. **Interaction**: Real-time manipulation and analysis in the browser

### **Backend Architecture**
- **Flask Framework**: Lightweight, flexible web framework
- **SQLAlchemy ORM**: Database abstraction and management
- **RDKit Integration**: Chemical informatics and molecular processing
- **RESTful API**: Clean API endpoints for data exchange

### **Frontend Stack**
- **3Dmol.js**: WebGL-based 3D molecular visualization
- **Bootstrap 5**: Responsive CSS framework
- **Font Awesome**: Professional icon library
- **Custom CSS**: Modern dark theme with glassmorphism effects

### **Database Schema**
```sql
Molecule:
├── id (Primary Key)
├── name (String)
├── smiles (String)
├── formula (String)
├── molar_mass (Float)
├── file_path (String)
├── created_at (DateTime)
└── properties (JSON)
```

## 🚀 Advanced Usage

### **API Endpoints**
```http
GET  /api/molecule/<id>/data     # Get molecule 3D data
POST /api/molecule/<id>/save     # Save to collection
GET  /api/molecule/<id>/download # Download as SDF
```

### **Development Mode**
```bash
# Enable Flask debug mode
export FLASK_ENV=development
export FLASK_DEBUG=1
python app.py
```

### **Database Management**
```bash
# Clear all data
flask clear-db

# Reinitialize database
flask init-db
```

## 🤝 Contributing

We welcome contributions! Here's how to get started:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### **Development Guidelines**
- Follow PEP 8 style guidelines
- Add tests for new features
- Update documentation as needed
- Ensure cross-platform compatibility

## 🐛 Troubleshooting

### **Common Issues**

#### **Port Already in Use**
```bash
# Kill process using port 5001
lsof -i :5001
kill <PID>
```

#### **RDKit Installation Issues**
```bash
# Reinstall RDKit via conda
conda remove rdkit
conda install -c conda-forge rdkit
```

#### **Database Errors**
```bash
# Reset database
rm instance/molecules.db
flask init-db
```

#### **3D Visualization Not Loading**
- Check browser console for JavaScript errors
- Ensure WebGL is supported in your browser
- Try a different browser (Chrome, Firefox recommended)

## 📋 Requirements

### **System Requirements**
- **OS**: Windows 10+, macOS 10.14+, Linux (Ubuntu 18.04+)
- **RAM**: 4GB minimum, 8GB recommended
- **Browser**: Chrome 80+, Firefox 75+, Safari 13+, Edge 80+

### **Python Dependencies**
See `requirements.txt` for complete list:
- Flask 2.3.3
- Flask-SQLAlchemy 3.0.5
- RDKit 2023.3.3
- NumPy 1.24.3
- Pillow 10.0.0

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🌟 Acknowledgments

- **RDKit Team** - For the incredible computational chemistry toolkit
- **3Dmol.js Developers** - For the powerful 3D visualization library
- **Flask Community** - For the excellent web framework
- **Chemical Community** - For continuous feedback and support

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/rhr18818/chem_3d_model/issues)
- **Documentation**: [Project Wiki](https://github.com/rhr18818/chem_3d_model/wiki)
- **Email**: [rhr18818@example.com](mailto:rhr18818@example.com)

---

<div align="center">

**Made with ❤️ for the Chemistry Community**

[⭐ Star this project](https://github.com/rhr18818/chem_3d_model) | [🐛 Report Bug](https://github.com/rhr18818/chem_3d_model/issues) | [💡 Request Feature](https://github.com/rhr18818/chem_3d_model/issues)

</div>
