import os, re
import shutil

def clean(path):
    os.chdir(path)
    delete = [".aux", ".fdb_latexmk", ".fls", ".log", ".maf", 
                ".out", ".synctex.gz", ".toc", ".mtc", '.bbl', '.blg', '.nav', '.snm', '.hd']
    for file in os.listdir():
        if any([file.endswith(i) for i in delete]):
            os.remove(file)
        elif re.match('.+\.mtc\d+', file):
            os.remove(file)
        elif os.path.isdir(file) and not file.startswith("."):
            if file.startswith("_minted"):
                shutil.rmtree(file)
            else:
                clean(path+"/"+file)
                os.chdir(path)

if __name__ == "__main__":
    path = os.path.dirname(os.path.abspath(__file__))
    clean(path)
