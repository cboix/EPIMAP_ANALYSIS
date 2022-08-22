#!/usr/bin/python
import fire
import encode_utils as eu
from encode_utils.connection import Connection

class encode_handler(object):
    def __init__(self, fileid, filepath, mode='prod'):
        self.fileid = fileid
        self.filepath = filepath
        self.mode = mode

    def upload(self):
        conn = Connection(self.mode)
        conn.upload_file(file_id=self.fileid,
                         file_path=self.filepath)

if __name__ == "__main__":
    fire.Fire(encode_handler)
