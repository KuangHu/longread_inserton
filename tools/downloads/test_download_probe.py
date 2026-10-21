import argparse
import gzip
import hashlib
from io import BytesIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import download_probe as d


class Response(BytesIO):
    status = 206
    headers = {"Content-Range": "bytes 0-2/3", "ETag": "test"}


class IntegrityTests(unittest.TestCase):
    def test_range_bytes(self):
        with patch.object(d,"request",return_value=Response(b"abc")):
            data, info=d.fetch_range("https://example.test/data",0,2)
            self.assertEqual(data,b"abc")
            self.assertEqual(info["bytes"],3)

    def test_ignored_range_rejected(self):
        response=Response(b"abc")
        response.status=200
        with patch.object(d,"request",return_value=response):
            with self.assertRaisesRegex(ValueError,"honor Range"):
                d.fetch_range("https://example.test/data",0,2)

    def test_incomplete_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/"bad.gz"
            with gzip.open(p,"wt") as f:f.write("@read\nACGT\n+\n")
            with self.assertRaises(EOFError):list(d.read_fastq(p))

    def test_sequence_quality_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/"bad.gz"
            with gzip.open(p,"wt") as f:f.write("@read\nACGT\n+\nIII\n")
            with self.assertRaises(ValueError):list(d.read_fastq(p))

    def test_full_checksum_before_publication(self):
        payload=gzip.compress(b"@read\nACGT\n+\nIIII\n")
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)
            resolved=dict(url="https://example.test/full.fastq.gz",source_format="FASTQ.GZ",
                          source_size=len(payload),archive_md5="0"*32)
            args=argparse.Namespace(out=out,max_full_gb=1)
            with patch.object(d,"resolve",return_value=resolved),patch.object(d,"request",return_value=Response(payload)):
                with self.assertRaisesRegex(ValueError,"MD5"):
                    d.full_download({},args)
            self.assertFalse((out/"full.fastq.gz").exists())
            resolved["archive_md5"]=hashlib.md5(payload).hexdigest()
            with patch.object(d,"resolve",return_value=resolved),patch.object(d,"request",return_value=Response(payload)):
                d.full_download({},args)
            self.assertEqual(list(d.read_fastq(out/"full.fastq.gz")),[("@read","ACGT","IIII")])


if __name__=="__main__":unittest.main()

