import gau2grid as gg
from pathlib import Path
import os
import zipfile
import tempfile
import shutil

am_list = [6, 8]

for am in am_list:

    folder = f"gau2grid-am{am}-{gg.__version__}"
    zip_filename = folder + '.zip'
    zipf = zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED)

    path = Path(folder)
    path.mkdir(parents=True)
    gg.c_gen.generate_c_gau2grid(am, path=path.resolve())

    for filename in path.iterdir():
        zipf.write(filename)

    shutil.rmtree(path.resolve())

#    with tempfile.TemporaryDirectory() as tmp:
#        os.chdir(tmp)
#
#        folder = f"gau2grid-am{am}-{gg.__version__}"
#        zip_filename = folder + '.zip'
#        zip_path = os.path.join(tmp, zip_filename)
#        zipf = zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED)
#
#        path = Path(tmp) / folder
#        path.mkdir(parents=True)
#        gg.c_gen.generate_c_gau2grid(am, path=path.resolve())
#
#        for filename in path.iterdir():
#            arcname = os.path.join(*str(filename).split(os.path.sep)[-2:])
#            print(filename, arcname)
#            zipf.write(filename, arcname=arcname)
#        
