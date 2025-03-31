import os
import pickle
from pathlib import Path
from logging import getLogger, handlers, Formatter, DEBUG
from scipy.stats import norm
from sklearn.preprocessing import robust_scale

def pickle_dump(obj, path):
    with open(path, mode='wb') as f:
        pickle.dump(obj,f)

def pickle_load(path):
    with open(path, mode='rb') as f:
        data = pickle.load(f)
        return data

def set_logger():
    root_logger = getLogger()
    root_logger.setLevel(DEBUG)
    rotating_handler = handlers.RotatingFileHandler(
        r'../../log/app.log',
        mode="a",
        maxBytes=100 * 1024,
        backupCount=3,
        encoding="utf-8"
    )
    format = Formatter('%(asctime)s : %(levelname)s : %(filename)s - %(message)s')
    rotating_handler.setFormatter(format)
    root_logger.addHandler(rotating_handler)

def robust_z(x):
    coefficient = norm.ppf(0.75) - norm.ppf(0.25)
    robust_z_score = robust_scale(x)*coefficient
    return robust_z_score

def file_checker(path, overwrite = False):
    path_obj = Path(path)
    if path_obj.exists():
        if overwrite:
            print(f"{path} is already exists, but OVERWRITE!")
            return True
        else:
            print(f"{path} is already exists!")
            return False
    else:
        return True