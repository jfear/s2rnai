import os
from pathlib import Path
import yaml


def read_config(fname, keepers=None):
    """Reads a YAML file.

    If a list of keepers is provided, will look through the YAML and only
    return those keys.

    """
    with open(fname, "r") as fh:
        c = yaml.full_load(fh)

    if keepers is None:
        return c

    if isinstance(keepers, str):
        return c.get(keepers, None)

    config = {}
    for k in keepers:
        v = c.get(k, None)
        if v is not None:
            config[k] = v

    return config


# Useful directories
PROJECT_DIR = Path(__file__).absolute().parents[2].as_posix()
CONFIG_DIR = Path(PROJECT_DIR, "config").as_posix()
CACHE_DIR = Path("~/.cache").expanduser().as_posix()

# Make missing dirs
Path(CACHE_DIR).mkdir(exist_ok=True, parents=True)
Path(CONFIG_DIR).mkdir(exist_ok=True, parents=True)
