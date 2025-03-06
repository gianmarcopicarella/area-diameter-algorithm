import os
import shutil
import json


def prepare_path(path):
    if os.path.isdir(path): shutil.rmtree(path)
    os.makedirs(path)


def to_kb(bytes):
    return bytes / 1000


def to_seconds(ms):
    return ms / 1000


def read_json(path, mode='r'):
    with open(path, mode) as f: return json.load(f)


def write_json(path, data):
    with open(path, 'w') as f: json.dump(data, f)
