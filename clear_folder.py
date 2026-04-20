import nbformat
import os

# Define the extensions we want to delete
EXTENSIONS_TO_DELETE = {".png", ".pdf",".df"}

for root, dirs, files in os.walk("."):
    # Skip hidden folders like .ipynb_checkpoints or .git
    dirs[:] = [d for d in dirs if not d.startswith('.')] 
    
    for file in files:
        path = os.path.join(root, file)
        
        # --- Task 1: Delete PNG and PDF files ---
        if any(file.lower().endswith(ext) for ext in EXTENSIONS_TO_DELETE):
            try:
                os.remove(path)
                print(f"Deleted: {path}")
            except OSError as e:
                print(f"Error deleting {path}: {e}")
            continue # Move to next file

        # --- Task 2: Clear Notebook Outputs ---
        if file.endswith(".ipynb"):
            print(f"Cleaning Notebook: {path}")
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    nb = nbformat.read(f, as_version=4)
                
                for cell in nb.cells:
                    if cell.cell_type == 'code':
                        cell.outputs = []
                        cell.execution_count = None
                
                with open(path, 'w', encoding='utf-8') as f:
                    nbformat.write(nb, f)
            except Exception as e:
                print(f"Failed to process notebook {path}: {e}")