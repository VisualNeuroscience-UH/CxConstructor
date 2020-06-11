import zlib
import pickle

def getData(filename):
    fi = open(filename, 'rb')
    data_pickle = zlib.decompress(fi.read())
    data = pickle.loads(data_pickle)
    return data
