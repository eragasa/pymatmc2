import pytest
import os
import shutil

from pymatmc2.utils import validate_vasp_simulation_input_files_exist

@pytest.fixture
def vasp_simulation_dir(tmpdir, request):
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    tmpdir.mkdir('vasp_simulation')
    if os.path.isdir(test_dir):
        for f in ['POSCAR', 'INCAR', 'POTCAR', 'KPOINTS']:
            shutil.copy(
                src = os.path.join(test_dir, f),
                dst = os.path.join(str(tmpdir), 'vasp_simulation', f)
            )
    
    return tmpdir

def test__files_do_exist(vasp_simulation_dir, request):
    path = os.path.join(str(vasp_simulation_dir), 'vasp_simulation')
    assert validate_vasp_simulation_input_files_exist(path=path)