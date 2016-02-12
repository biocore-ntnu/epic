
from joblib import Memory
from os.path import expanduser

MEMORY = Memory(cachedir=expanduser("~/.epic"), verbose=0)
