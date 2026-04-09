import nbformat
import os

for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if not d.startswith('.')] # Skip hidden folders
    for file in files:
        if file.endswith(".ipynb"):
            path = os.path.join(root, file)
            print(f"Cleaning: {path}")
            with open(path, 'r', encoding='utf-8') as f:
                nb = nbformat.read(f, as_version=4)
            
            # Manually clear outputs and execution counts
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    cell.outputs = []
                    cell.execution_count = None
            
            with open(path, 'w', encoding='utf-8') as f:
                nbformat.write(nb, f)