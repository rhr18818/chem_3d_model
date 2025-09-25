// 3Dmol viewer integration

class MoleculeViewer {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.options = {
            width: options.width || 800,
            height: options.height || 600,
            backgroundColor: options.backgroundColor || 'white',
            ...options
        };
        
        this.viewer = null;
        this.isSpinning = false;
        this.init();
    }
    
    init() {
        this.viewer = $3Dmol.createViewer(this.containerId, {
            width: this.options.width,
            height: this.options.height,
            backgroundColor: this.options.backgroundColor
        });
    }
    
    loadMolecule(data, format = 'sdf') {
        this.viewer.clear();
        this.viewer.addModel(data, format);
        this.viewer.zoomTo();
        this.viewer.render();
        return this;
    }
    
    applyStyle(styleOptions = {}) {
        const defaultStyle = {
            stick: {
                radius: 0.15,
                colorscheme: 'Jmol'
            },
            sphere: {
                radius: 0.3
            }
        };
        
        const style = {...defaultStyle, ...styleOptions};
        this.viewer.setStyle({}, style);
        this.viewer.render();
        return this;
    }
    
    addLabels(labels = []) {
        labels.forEach(label => {
            this.viewer.addLabel(label.text, {
                position: label.position,
                backgroundColor: label.backgroundColor || 'white',
                fontColor: label.fontColor || 'black',
                fontSize: label.fontSize || 14
            });
        });
        this.viewer.render();
        return this;
    }
    
    toggleSpin() {
        this.isSpinning = !this.isSpinning;
        this.viewer.spin(this.isSpinning);
        return this.isSpinning;
    }
    
    resetView() {
        this.viewer.zoomTo();
        this.viewer.render();
        return this;
    }
    
    downloadImage(filename = 'molecule.png') {
        this.viewer.pngURI().then(pngURI => {
            const link = document.createElement('a');
            link.href = pngURI;
            link.download = filename;
            link.click();
        });
        return this;
    }
    
    setBackgroundColor(color) {
        this.viewer.setBackgroundColor(color);
        this.viewer.render();
        return this;
    }
}

// Export for use in other modules
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MoleculeViewer;
} else if (typeof window !== 'undefined') {
    window.MoleculeViewer = MoleculeViewer;
}