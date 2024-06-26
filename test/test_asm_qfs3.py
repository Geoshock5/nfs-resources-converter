import io
import os
import unittest

from resources.eac.compressions.qfs3 import Qfs3Compression


class TestAsmQFS3Algorythm(unittest.TestCase):

    def test_0_al1(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/AL1.QFS'
        with open(file_name, 'rb') as file:
            uncompressed = parser.uncompress(file, os.path.getsize(file_name))
            print('CHECKING OUTPUT....')
            with open('test/samples/AL1.FSH', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])

    def test_1_vertbst(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/VERTBST.QFS'
        with open(file_name, 'rb') as file:
            uncompressed = parser.uncompress(file, os.path.getsize(file_name))
            print('CHECKING OUTPUT....')
            with open('test/samples/VERTBST.FSH', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])

    def test_2_ldiabl_pbs(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/LDIABL.PBS'
        with open(file_name, 'rb') as file:
            uncompressed = parser.uncompress(file, os.path.getsize(file_name))
            print('CHECKING OUTPUT....')
            with open('test/samples/LDIABL.PBS.BIN', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])

    def test_3_gtitle(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/GTITLE.QFS'
        with open(file_name, 'rb') as file:
            uncompressed = parser.uncompress(file, os.path.getsize(file_name))
            print('CHECKING OUTPUT....')
            with open('test/samples/GTITLE.FSH', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])

    def test_4_gvertbst(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/GVERTBST.QFS'
        with open(file_name, 'rb') as file:
            uncompressed = parser.uncompress(file, os.path.getsize(file_name))
            print('CHECKING OUTPUT....')
            with open('test/samples/GVERTBST.FSH', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])

    def test_5_ldiabl_pbs_compression(self):
        parser = Qfs3Compression()
        file_name = 'test/samples/LDIABL.PBS.BIN'
        with open(file_name, 'rb') as file:
            compressed = parser.compress(file, os.path.getsize(file_name))
            comp_stream = io.BytesIO(compressed)
            uncompressed = parser.uncompress(comp_stream, len(compressed))
            print('CHECKING OUTPUT....')
            with open('test/samples/LDIABL.PBS.BIN', 'rb') as fsh_file:
                fsh = fsh_file.read()
                self.assertEqual(len(fsh), len(uncompressed))
                for i in range(len(fsh)):
                    if i % 10000 == 0:
                        print(f"{i}/{len(fsh)}")
                    self.assertEqual(fsh[i], uncompressed[i])