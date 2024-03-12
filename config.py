import os


RAMA_PREFERENCES = {
    "PRO": {
        "file": os.path.join('data', 'pref_proline.data'),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRE-PRO": {
        "file": os.path.join('data', 'pref_preproline.data'),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "GLY": {
        "file": os.path.join('data', 'pref_glycine.data'),
        "bounds": [0, 0.002, 0.02, 1],
    }
}
