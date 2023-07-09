import unittest

from library import require_file


class TestSerializeDeserialize(unittest.TestCase):

    def test_tri_should_remain_the_same(self):
        tri_map = require_file('test/samples/AL1.TRI')
        output = tri_map.to_bytes()
        with open('test/samples/AL1.TRI', 'rb') as bdata:
            original = bdata.read()
            self.assertEqual(len(original), len(output))
            for i, x in enumerate(original):
                self.assertEqual(x, output[i], f"Wrong value at index {i}")

    def test_fsh_should_remain_the_same(self):
        fsh = require_file('test/samples/VERTBST.FSH')
        output = fsh.to_bytes()
        with open('test/samples/VERTBST.FSH', 'rb') as bdata:
            original = bdata.read()
            # self.assertEqual(len(original), len(output))
            for i, x in enumerate(original):
                self.assertEqual(x, output[i], f"Wrong value at index {i}")

    def test_cfm_should_remain_the_same(self):
        car_fam = require_file('test/samples/LDIABL.CFM')
        output = car_fam.to_bytes()
        with open('test/samples/LDIABL.CFM', 'rb') as bdata:
            original = bdata.read()
            # self.assertEqual(len(original), len(output))
            for i, x in enumerate(original):
                self.assertEqual(x, output[i], f"Wrong value at index {i}")

    def test_ffn_should_remain_the_same(self):
        font_res = require_file('test/samples/MAIN24.FFN')
        output = font_res.to_bytes()
        with open('test/samples/MAIN24.FFN', 'rb') as bdata:
            original = bdata.read()
            self.assertEqual(len(original), len(output))
            for i, x in enumerate(original):
                self.assertEqual(x, output[i], f"Wrong value at index {i}")

    def test_ffn_can_be_reconstructed_from_files(self):
        font_res = require_file('test/samples/MAIN24.FFN')
        import tempfile
        from serializers import get_serializer
        serializer = get_serializer(font_res.block)
        self.assertTrue(serializer.setup_for_reversible_serialization())
        with tempfile.TemporaryDirectory() as tmp:
            serializer.serialize(font_res, tmp)
            serializer.deserialize(tmp, font_res)
        output = font_res.to_bytes()
        with open('test/samples/MAIN24.FFN', 'rb') as bdata:
            original = bdata.read()
            self.assertEqual(len(original), len(output))
            for i, x in enumerate(original):
                self.assertEqual(x, output[i], f"Wrong value at index {i}")
