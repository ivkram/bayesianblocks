import os
import glob

def find_gamma(file, radio_name):
    file.seek(0)
    for s in file:
        s = s.split()
        if s[0] == radio_name:
            return s[3][:12]
    return False

def get_J2000(file, radio_name):
    file.seek(0)
    for s in file:
        s = s.split()
        if s[0] == radio_name:
            return s[1]
