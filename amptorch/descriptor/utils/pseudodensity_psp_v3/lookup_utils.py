import os
import json


def generate_lookup_dictionary():
    lookup = {}
    for filename in os.listdir("."):
        if filename.endswith(".g"):
            element = filename.split("_")[0]
            lookup[element] = filename

    with open("lookup.json", "w") as f:
        json.dump(lookup, f)


generate_lookup_dictionary()
