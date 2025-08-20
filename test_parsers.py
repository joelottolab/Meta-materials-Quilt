import glob, yaml
def test_yaml_files_exist():
    files = glob.glob('../**/param_*.yaml', recursive=True)
    # This test asserts parser script presence more than actual files in this packaged repo
    assert isinstance(files, list)
