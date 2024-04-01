from io import BufferedReader

import struct
import heapq

from library.utils import read_short, read_int
from library.utils.asm_runner import AsmRunner
from resources.eac.compressions.base import BaseCompressionAlgorithm

# Credits:
# Original research / file format by AndyGura (github.com/AndyGura)
# Uncompression code rewritten by Geoshock5 based on code from AndyGura
# Compression method by Geoshock5 (github.com/Geoshock5)
#
# Huffman node/table code from Aashish Barnwal / GeeksForGeeks (https://www.geeksforgeeks.org/huffman-coding-greedy-algo-3/)
# License: CCBY-SA

# Node class. Used by Huffman tree code.
class node:  
    def __init__(self, freq, symbol, left=None, right=None): 
        # frequency of symbol 
        self.freq = freq 
  
        # symbol name (character) 
        self.symbol = symbol 
  
        # node left of current node 
        self.left = left 
  
        # node right of current node 
        self.right = right 
  
        # tree direction (0/1) 
        self.huff = '' 
  
    def __lt__(self, nxt): 
        return self.freq < nxt.freq 

# utility function to print huffman 
# codes for all symbols in the newly 
# created Huffman tree 
def printNodes(node, val=''):

    # huffman code for current node 
    newVal = val + str(node.huff)

    # if node is not an edge node 
    # then traverse inside it 
    if(node.left): 
        printNodes(node.left, newVal) 
    if(node.right): 
        printNodes(node.right, newVal) 

        # if node is edge node then 
        # display its huffman code 
    if(not node.left and not node.right): 
        print(f"{node.symbol} -> {newVal}") 

# Derived function from printNodes. Needed to populate the Huffman dictionary.
def setNodes(ndict, node, val=''):

    # huffman code for current node 
    newVal = val + str(node.huff)

    # if node is not an edge node 
    # then traverse inside it 
    if(node.left): 
        setNodes(ndict, node.left, newVal) 
    if(node.right): 
        setNodes(ndict, node.right, newVal)

        # if node is edge node then 
        # record its huffman code 
    if(not node.left and not node.right): 
        ndict[node.symbol] = newVal

# Bitpack algorithm used in RefPack. Short vals <=3 are 3-bits long with a leading bit set to 1 (i.e. 0 = 100, 3 = 111).
# Values >3 use leading zeroes to flag the length in bits of what follows (i.e. len = 3 + num. of 0s)
def pack_bits(val: int):
    if val<4:
        len = 3
        return len, (val+4) # return with flag bit 1xx set
    else:
        len = (val+4).bit_length()
        return len, (val+4) # return with flag bit 0...01xx set

class Qfs3Compression(BaseCompressionAlgorithm, AsmRunner):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, asm_virtual_memory_size=2 * 1024, **kwargs)
        self.output_length = 0
        self.index_table_0 = []
        self.unk_table = [0] * 256    # values to output to uncompressed file
        self.unk_table_2 = [0] * 256  # length of huffman tree code
        self.unk_table_3 = [0] * 256  # comparators / codes for each level.  Used to check the length in bits of the current value
        self.accumulator = 0
        self.available_acc_bits = 0
        self.huff_dict = [None] * 257

    def append_to_output(self, buffer, value):
        buffer.extend(value.to_bytes(1, 'little'))

    def reuse_output_byte(self, buffer, length):
        value = buffer[-1]
        for i in range(length):
            self.append_to_output(buffer, value)

    def read_next(self, buffer):
        self.accumulator = (read_short(buffer, 'big') | (self.accumulator << 16))

    def accumulate_if_needed(self, buffer):
        if self.available_acc_bits < 0:
            self.read_next(buffer)
            self.esi = self.accumulator << -self.available_acc_bits
            self.available_acc_bits += 16

    def uncompress(self, buffer: BufferedReader, input_length: int) -> bytes:
        uncompressed: bytearray = bytearray()

        table_110 = [0] * 16

        self.define_variable('var_14C', -0x14C, 4)
        self.define_variable('var_154', -0x154, 4)

        index_table_0_capacity = 0
        unk_counter = 0

        self.available_acc_bits = 0
        file_header = self.accumulator = read_short(buffer, 'big')
        self.read_next(buffer)
        self.esi = self.accumulator << 16
        # if compressed size presented
        if file_header & 0x100:
            self.edx = read_int(buffer, 'big')
            self.available_acc_bits = 8
            self.esi = self.edx << 8
            self.accumulator = self.edx
            # reset size presented flag bit
            file_header = file_header & 0xFEFF
        self.ebx = self.esi >> 0x18
        self.available_acc_bits -= 8
        self.esi = self.esi << 8
        self.accumulate_if_needed(buffer)
        self.eax = self.esi >> 16
        self.available_acc_bits -= 16
        self.esi = self.esi << 16
        self.output_length = self.eax
        self.accumulate_if_needed(buffer)
        self.edx = self.output_length = self.output_length | (self.ebx << 16)
        self.ebx = 1
        self.eax = self.esi >> 0x18
        self.available_acc_bits -= 8
        self.esi = self.esi << 8
        unk_1 = self.al
        self.accumulate_if_needed(buffer)
        unk_90 = 1
        unk_3 = 15
        unk_2 = 4
        self.ebp = 0
        while True:
            self.ebp = self.ebp << 1
            self.ecx = index_table_0_capacity
            self.eax = self.ebp - self.ecx
            self.edx = unk_2
            assert self.edx % 4 == 0  # logic below will work in different way if not
            self.unk_table_3[int(self.edx / 4)] = self.eax
            if self.get_register_signed_value('esi') >= 0:
                self.eax = 2
                if (self.esi >> 16) == 0:
                    self.edx = 0
                    while self.edx == 0:
                        self.edx = self.esi >> 0x1F
                        self.eax += 1
                        self.available_acc_bits -= 1
                        self.esi = self.esi << 1
                        self.accumulate_if_needed(buffer)
                else:
                    while self.get_register_signed_value('esi') >= 0:
                        self.esi = self.esi << 1
                        self.eax += 1
                    self.edx = self.eax - 1
                    self.available_acc_bits -= self.edx
                    self.esi = self.esi << 1
                    self.accumulate_if_needed(buffer)
                if self.get_register_signed_value('eax') <= 16:
                    self.edx = self.esi >> (0x20 - self.eax)
                    self.ecx = self.al
                    self.available_acc_bits -= self.eax
                    self.esi = self.esi << self.cl
                    self.accumulate_if_needed(buffer)
                else:
                    self.ebx = self.eax - 16
                    self.edx = self.esi >> (0x20 - self.ebx)
                    self.ecx = self.bl
                    self.available_acc_bits -= self.ebx
                    self.esi = self.esi << self.cl
                    self.accumulate_if_needed(buffer)
                    self.ebx = self.esi >> 16
                    self.available_acc_bits -= 16
                    self.esi = self.esi << 16
                    self.accumulate_if_needed(buffer)
                    self.edx = (self.edx << 16) | self.ebx
                self.edx += (1 << self.al)
            else:
                self.edx = self.esi >> 0x1D
                self.available_acc_bits -= 3
                self.esi = self.esi << 3
                self.accumulate_if_needed(buffer)
            self.edx -= 4
            self.eax = unk_2
            assert self.eax % 4 == 0
            table_110[int(self.eax / 4)] = self.edx
            self.eax = index_table_0_capacity
            self.ebp += self.edx
            self.eax += self.edx
            self.ecx = 0
            index_table_0_capacity = self.eax
            if self.edx != 0:
                self.cl = unk_3 & 0xFF
                self.eax = self.ebp << self.cl
                self.ecx = self.eax & 0xFFFF
            self.ebx = unk_3 - 1
            self.eax = unk_2 = unk_2 + 4
            unk_3 = self.ebx
            self.ebx = unk_90 + 1
            self.set_value('[esp+eax+550h+var_154]', self.ecx)
            unk_90 = self.ebx
            if self.edx == 0:
                continue
            if self.ecx == 0:
                break
        self.ecx = 0xFFFFFFFF
        unk_8 = self.ebx - 1
        self.set_value('[esp+ebx*4+550h+var_154]', self.ecx)
        self.ebx = 0
        huff_table_0_iter_ptr = 0
        self.eax = 0xFF
        unk_4 = 0
        if index_table_0_capacity > 0:
            self.index_table_0 = [None] * index_table_0_capacity
            while True:
                if self.get_register_signed_value('esi') >= 0:
                    self.edx = 2
                    if (self.esi >> 16) == 0:
                        while True:
                            self.ebp = self.esi
                            self.edx += 1
                            self.available_acc_bits -= 1
                            self.ebp = self.ebp >> 0x1F
                            self.esi = self.esi << 1
                            self.accumulate_if_needed(buffer)
                            if self.ebp != 0:
                                break
                    else:
                        while self.get_register_signed_value('esi') >= 0:
                            self.esi = self.esi << 1
                            self.edx += 1
                        self.ebx = self.edx - 1
                        self.available_acc_bits -= self.ebx
                        self.esi = self.esi << 1
                        self.accumulate_if_needed(buffer)
                    self.ebx = self.accumulator << 8
                    if self.get_register_signed_value('edx') <= 16:
                        self.ebx = self.esi >> (0x20 - self.edx)
                        self.cl = self.dl
                        self.available_acc_bits -= self.edx
                        self.esi = self.esi << self.cl
                        self.accumulate_if_needed(buffer)
                    else:
                        unk_4 = self.edx - 16
                        self.ecx = 0x20 - unk_4
                        self.ebx = self.esi >> self.cl
                        self.cl = unk_4
                        self.esi = self.esi << self.cl
                        self.available_acc_bits -= unk_4
                        self.accumulate_if_needed(buffer)
                        self.ecx = self.esi >> 16
                        self.available_acc_bits -= 16
                        self.esi = self.esi << 16
                        unk_7 = self.ecx
                        self.accumulate_if_needed(buffer)
                        self.ecx = unk_7
                        self.ebx = (self.ebx << 16) | self.ecx
                    self.cl = self.dl
                    self.edx = 1 << self.cl
                    self.ebx += self.edx
                else:
                    self.ebx = self.esi >> 0x1D
                    self.available_acc_bits -= 3
                    self.esi = self.esi << 3
                    self.accumulate_if_needed(buffer)
                self.ebx -= 3
                while self.ebx != 0:
                    self.al += 1
                    self.edx = self.al
                    if self.edx not in self.index_table_0:
                        self.ebx -= 1
                self.edx = self.al
                self.edx = huff_table_0_iter_ptr
                self.ecx = index_table_0_capacity
                self.ebx = self.edx + 1
                self.index_table_0[self.edx] = self.al
                huff_table_0_iter_ptr = self.ebx
                if self.get_register_signed_value('ebx') >= self.get_register_signed_value('ecx'):
                    break
        self.unk_table_2 = [0x40] * 256
        self.edx = 0
        self.ebx = 0
        self.ecx = unk_8
        unk_5 = 1
        unk_9 = 0
        if self.ecx >= 1:
            self.ecx = 7
            self.eax = 4
            unk_2 = self.ecx
            unk_3 = self.eax
            while True:
                self.eax = unk_3
                assert self.eax % 4 == 0
                self.eax = table_110[int(self.eax / 4)]
                self.ebp = unk_5
                unk_6 = self.eax
                if self.get_register_signed_value('ebp') >= 9:
                    break
                self.ecx = unk_2
                self.ebp = 1 << self.cl
                break_outer = False
                continue_outer = False
                while True:
                    self.eax = unk_6 - 1
                    unk_6 = self.eax
                    if self.eax == 0xFFFFFFFF:
                        self.ebp = unk_3 = unk_3 + 4
                        self.eax = unk_2 = unk_2 - 1
                        self.ecx = unk_5 + 1
                        unk_5 = self.ecx
                        self.ebp = unk_8
                        if self.get_register_signed_value('ecx') <= self.get_register_signed_value('ebp'):
                            continue_outer = True
                            break
                        else:
                            break_outer = True
                            break
                    else:
                        index_table_0_value = self.index_table_0[unk_counter]
                        unk_counter += 1
                        unk_0 = unk_5
                        self.ecx = index_table_0_value
                        self.eax = unk_1
                        if self.eax == self.ecx:
                            self.eax = unk_5
                            unk_9 = self.eax
                            unk_0 = 0x60
                        self.eax = 0
                        if self.get_register_signed_value('ebp') <= 0:
                            continue
                        while self.get_register_signed_value('eax') < self.get_register_signed_value('ebp'):
                            self.edx += 1
                            self.cl = index_table_0_value
                            self.ebx += 1
                            self.unk_table[self.edx - 1] = self.cl
                            self.cl = unk_0
                            self.eax += 1
                            self.unk_table_2[self.ebx - 1] = self.cl
                if break_outer:
                    break
                if continue_outer:
                    continue
                break
        while True:
            if len(uncompressed) > self.output_length:
                raise Exception('Uncompress algorythm writes more that file length')
            self.eax = self.unk_table_2[self.esi >> 0x18]
            self.available_acc_bits -= self.eax
            while self.available_acc_bits >= 0:
                self.edx = self.esi >> 0x18
                for _ in range(4):
                    self.append_to_output(uncompressed, self.unk_table[self.edx])
                    self.esi = self.esi << self.al
                    self.edx = self.esi >> 0x18
                    self.eax = self.unk_table_2[self.edx]
                    self.available_acc_bits -= self.eax
                    if self.available_acc_bits < 0:
                        break
            self.available_acc_bits += 0x10
            if self.available_acc_bits >= 0:
                self.append_to_output(uncompressed, self.unk_table[(self.esi >> 0x18)])
                self.read_next(buffer)
                self.esi = self.accumulator << (0x10 - self.available_acc_bits)
                continue
            self.available_acc_bits += self.eax - 0x10
            if self.eax == 0x60:
                self.eax = unk_9
            else:
                self.eax = 8
                self.edx = self.esi >> 16
                self.ecx = 0x20
                while True:
                    self.eax += 1
                    self.ebp = self.get_value('[esp+ecx+550h+var_14C]')[0]
                    self.ecx += 4
                    if self.edx < self.ebp:
                        break
            self.ecx = 0x20 - self.eax
            self.edx = self.esi >> self.cl
            self.cl = self.al
            self.available_acc_bits -= self.eax
            self.esi = self.esi << self.cl
            self.ecx = self.unk_table_3[self.eax]
            self.eax = self.edx - self.ecx
            self.al = self.index_table_0[self.eax]
            if self.al != unk_1:
                if self.available_acc_bits >= 0:
                    self.append_to_output(uncompressed, self.al)
                    continue
            self.accumulate_if_needed(buffer)
            if self.al != unk_1:
                self.append_to_output(uncompressed, self.al)
                continue
            if self.get_register_signed_value('esi') >= 0:
                self.eax = 2
                if (self.esi >> 16) == 0:
                    self.ebp = 0
                    while unk0 == 0:
                        unk0 = self.esi >> 0x1F
                        self.eax += 1
                        self.available_acc_bits -= 1
                        self.esi = self.esi << 1
                        self.accumulate_if_needed(buffer)
                else:
                    while True:
                        self.esi = self.esi << 1
                        self.eax += 1
                        if self.get_register_signed_value('esi') < 0:
                            break
                    self.ecx = self.eax - 1
                    self.available_acc_bits -= self.ecx
                    self.esi = self.esi << 1
                    self.accumulate_if_needed(buffer)
                if self.get_register_signed_value('eax') <= 16:
                    self.ecx = 0x20 - self.eax
                    self.ebp = self.esi >> self.cl
                    self.available_acc_bits -= self.eax
                    self.cl = self.al
                    fill_bytes_length = self.ebp
                    self.esi = self.esi << self.cl
                    self.accumulate_if_needed(buffer)
                    self.cl = self.al
                    self.eax = (1 << self.cl) + fill_bytes_length
                else:
                    self.ecx = self.eax - 16
                    self.ebp = self.esi >> (0x20 - self.ecx)
                    self.esi = self.esi << self.cl
                    self.available_acc_bits -= self.ecx
                    unk_4 = self.ecx
                    fill_bytes_length = self.ebp
                    self.accumulate_if_needed(buffer)
                    self.ecx = self.esi >> 16
                    self.available_acc_bits -= 16
                    self.esi = self.esi << 16
                    unk3 = self.ecx
                    self.accumulate_if_needed(buffer)
                    self.ecx = fill_bytes_length << 16
                    self.eax = (1 << self.al) + (self.ecx | unk3)
                self.eax = self.eax - 4
                fill_bytes_length = self.eax
            else:
                self.eax = self.esi >> 0x1D
                self.available_acc_bits -= 3
                self.esi = self.esi << 3
                fill_bytes_length = self.eax
                self.accumulate_if_needed(buffer)
                fill_bytes_length -= 4
            self.ebp = fill_bytes_length
            if self.ebp == 0:
                self.ebp = self.esi >> 0x1F
                self.available_acc_bits -= 1
                self.esi = self.esi << 1
                self.accumulate_if_needed(buffer)
                if self.ebp != 0:
                    if file_header == 0x34FB:
                        value_b = value_a = 0
                        for i in range(self.output_length):
                            value_a += uncompressed[i]
                            value_b += value_a
                            uncompressed[i] = value_b & 0xFF
                    elif file_header == 0x32FB:
                        value = 0
                        for i in range(self.output_length):
                            value += uncompressed[i]
                            uncompressed[i] = value & 0xFF
                    break
                else:
                    self.eax = self.esi >> 0x18
                    self.available_acc_bits -= 8
                    self.esi = self.esi << 8
                    self.accumulate_if_needed(buffer)
                    self.append_to_output(uncompressed, self.eax & 0xFF)
            else:
                self.reuse_output_byte(uncompressed, self.ebp)
        return uncompressed
    
    def uncompress_new(self, buffer: BufferedReader, input_length: int) -> bytes:
        uncompressed: bytearray = bytearray()

        # A simplified version of the decompression algorithm by Geoshock5 for TNFS PBS files.  Used to help reverse-engineer the compression algorithm.
        #
        # QFS3 file structure overview for TNFS PBS files:
        # Header: similar to standard RefPack.
        # Some QFS3 versions potentially have compression bitflag, which affects the header structure, but it's not seen in TNFS so not documented here.
        # 1. 0x30FB magic (2b)
        # 2. File size of the uncompressed output in bytes, big-endian (3b). In some files, compressed size follows in same format.
        # 3. Number of characters in the uncompressed file (1b)
        #
        # Chunk 1: Bit-packed data defining the number of leaves on each level of the Huffman tree (i.e. val = n values of length = 1 bit, 2 bit, 3 bit ... 16-bit)
        #
        # Chunk 2: Character table. Chunk 3 is a canonical Huffman tree, using consecutive values for each node (0~255).  This second chunk translates that into a lookup table,
        # presumably sorted in something approximating probability order.  The last value (matching header item 3) is an escape character used to flag special functions,
        # e.g. repeat the last byte x times, End Of File, or do some other special things based on the following bit(s).
        # This chunk is encoded in bitpacked form as an offset from the previous value (i.e. [val+1] = [val] + offset.  But also, if we've previously seen a value between val and val+1 in this list, don't decrement the offset counter.
        #
        # Chunk 3: These are the bitpacked values of the compressed file.  Read them in one-by-one and then use the dictionary from chunk 2 to write them to the output, processing any special bytes as required.
        #
        # Jon T / Geoshock5 / March 2024
        
        # 1: read the file header
        file_header = struct.unpack(">h", buffer.read(2))[0]
        out_size = int.from_bytes(struct.unpack("<BBB", buffer.read(3)),"big")
        char_count = struct.unpack("B", buffer.read(1))[0]

        # 2: read in chunk 1
        length = 1          # current length of Huff code in bits
        huff_key = 0
        value_count = 0     # current sum of values in tree
        val = 0
        val_check = 0

        buf = int.from_bytes(struct.unpack("<BBBB",buffer.read(4)),"big")
        
        while (length<16):
            huff_key = huff_key << 1
            self.huff_table_codes[length] = huff_key - value_count
            val, buf, self.sub_ptr = self.count_bits(buf, buffer, self.sub_ptr)
            # subtract 4 (for some reason - it happens a lot...)
            val -= 4
            self.huff_chars_per_level.append(val)
            huff_key += val
            value_count += val
            val_check = 0
            if(val != 0):
                # do the sum check to see if we're complete
                val_check = huff_key << (16 - length)
                val_check = val_check & 0xFFFF
            
            length += 1
            self.depth_table.append(val_check)

            if(val==0):
                continue
            if(val_check==0):
                break

        self.depth_table[-1] = 0xFFFFFFFF
        # print("Tree complete.  Total values: " + str(value_count) + ", max. length: " + str(length)) # DEBUG
        length -= 1

        # 3: read in the character list
        char_table = [None] * value_count
        current_char=0
        char_value = 0xFF

        while True:
            val, buf, self.sub_ptr = self.count_bits(buf, buffer, self.sub_ptr)
            val -=3 # TODO: why is this only 3 and not 4 per regular RefPack?
            while (val!=0):
                char_value = (char_value + 1) & 0xFF
                if char_value not in char_table:
                    val -=1
            char_table[current_char] = char_value
            current_char += 1
            if (current_char >= len(char_table)):
                break

        # 4: generate the lookup tables for output as needed
        # This was part of the original algorithm but is not used here as it makes the decompression longer and much more complex to understand.
        # Computationally this is probably more expensive but clock cycles are cheap these days vs the 90s ;-)

        # 5: read the input and write the uncompressed output
        huff_key = 0
        read_val_len = 1
        break_outer = 0     # Used for debug only

        while True:
            if len(uncompressed) > out_size:
                print("Error decompressing - length exceeds header!")
                break
            if len(uncompressed) == out_size:
                print("Reached compression length. Finished...")
                break
            if(break_outer==1):
                break
            while True:
                huff_key = buf >> 16
                read_val_len = 1
                while ((huff_key) > self.depth_table[read_val_len]):
                    read_val_len += 1
                val = buf >> 32-read_val_len
                val = (val - self.huff_table_codes[read_val_len]) & 0xFF
                if(val>len(char_table)): # DEBUG: check for errors - should not happen
                    print("Error reading file! Aborting...")
                    break_outer = 1
                    break
                if(char_table[val]!=char_count):
                    uncompressed.append(char_table[val] & 0xFF)
                    buf = (buf << read_val_len) & 0xFFFFFFFF
                    self.sub_ptr += read_val_len
                    self.sub_ptr, buf = self.refill_buf(buf, buffer, self.sub_ptr)
                    continue
                else:
                    fill_length = 0
                    buf = (buf << read_val_len) & 0xFFFFFFFF
                    self.sub_ptr += read_val_len
                    self.sub_ptr, buf = self.refill_buf(buf, buffer, self.sub_ptr)
                    fill_length, buf, self.sub_ptr = self.count_bits(buf, buffer, self.sub_ptr)
                    fill_length -= 4
                    if (fill_length == 0):
                        # read next bit. If 0, do stuff, else output next byte
                        fill_length = buf>>31
                        buf = buf << 1
                        self.sub_ptr += 1
                        if(fill_length!=0):
                            # TODO: Add for completeness.  In the original algorithm but doesn't actually do anything for car specs.
                            print("End of file flag reached. Breaking...")
                            break
                        else:
                            uncompressed.append((buf >> 24) & 0xFF)
                            buf = (buf << 8) & 0xFFFFFFFF
                            self.sub_ptr += 8
                            self.sub_ptr, buf = self.refill_buf(buf, buffer, self.sub_ptr)
                    else:
                        self.reuse_output_byte(uncompressed, fill_length)
        return(uncompressed)
    
    def compress(self, buffer: BufferedReader, input_length: int) -> bytes:
        compressed: bytearray = bytearray()
        
        # 0 set up variables
        byte_count = [0] * 257

        buffer.seek(0,2)
        in_len = buffer.tell()
        buffer.seek(0,0)

        # 1 count the bits in the source file, splitting out dupes
        last_val = None
        for _ in range(int(in_len)):
            val = buffer.read(1)[0]
            if(val!=last_val):
                byte_count[val] += 1
            else:
                byte_count[256] += 1
            last_val = val

        # 2 sort by count and then by number using a priority queue
        nodes = [] # initialize priority queue
        
        for i in range(257):
            if(byte_count[i]>0):
                heapq.heappush(nodes,node(byte_count[i],i))

        count_val = len(nodes)

        while len(nodes)>1:
            # sort all the nodes in ascending order
            # based on their frequency
            left = heapq.heappop(nodes)
            right = heapq.heappop(nodes)
        
            # assign directional value to these nodes 
            left.huff = 0
            right.huff = 1
        
            # combine the 2 smallest nodes to create 
            # new node as their parent 
            newNode = node(left.freq+right.freq, left.symbol+right.symbol, left, right) 
        
            heapq.heappush(nodes, newNode) 
        
        # Huffman Tree is ready! Create a dictionary
        setNodes(self.huff_dict, nodes[0])

        # 3 Rebuild the codes canonically
        # 4 generate dictionary of Huff values
        canon = []
        for n in range(len(self.huff_dict)):
            if(self.huff_dict[n]!=None):
                heapq.heappush(canon,(len(self.huff_dict[n]),n))

        code = 0
        last_length = 1
        length = 1
        self.huff_chars_per_level = [0] * 16

        while len(canon)>0:
            key = heapq.heappop(canon)
            length = key[0]
            if(length!=last_length):
                code += 1
                code = code << length-last_length
            elif(length>1):
                code+=1
            self.huff_dict[key[1]] = code
            self.huff_chars_per_level[length] += 1
            last_length = length

        # need escape val to be a character not used in the file
        escape_char=0
        while(self.huff_dict[escape_char]!=None):
            escape_char += 1

        # 5 write the output file
        # 5.1 header
        compressed.extend([0x30,0xFB]) # standard RefPack magic bytes
        compressed.extend(in_len.to_bytes(3, byteorder='big')) # uncompressed length - write as big-endian
        compressed.append(escape_char) # special escape char for repeated byte functions, EOF etc.
        # 5.2 tree depth
        outstr = ''
        for d in range(1,length+1):
            outlen, outval = pack_bits(self.huff_chars_per_level[d])
            if outlen>3:
                for _ in range(int(outlen-3)):
                    outstr += '0'
            outstr += str(bin(outval))[2:]
        while(len(outstr)>8):
            outval = int(outstr[:8],2)
            outstr=outstr[8:]
            compressed.extend(outval.to_bytes(1,byteorder='little'))

        # 5.3 alphabet
        last_val = 0xFF
        chars_added = []

        for n in range(len(self.huff_dict)):
            if(self.huff_dict[n]!=None):
                heapq.heappush(canon,(self.huff_dict[n],n))

        while(len(canon)>0):
            offset=0
            seek = last_val
            val = heapq.heappop(canon)[1]

            if(val==256):
                val=escape_char

            while seek!=val:
                seek = (seek+1) & 0xFF
                if (seek not in chars_added):
                    offset += 1
            chars_added.append(val)                        
            if(val==count_val):
                special_char_done = True
            last_val=val
            offset-=1 # HACK: Not clear why.  Normal RefPack is offset 4.
            outlen, outval = pack_bits(offset)
            if(outlen>3):
                while(outlen>3):
                    outstr+='0'
                    outlen-=1
            outstr+=str(bin(outval))[2:]

        while(len(outstr)>8):
            outval = int(outstr[:8],2)
            outstr=outstr[8:]
            compressed.extend(outval.to_bytes(1,byteorder='little'))

        # 5.4 encoded output
        buffer.seek(0,0)
        last_val = None
        repeat_len=0
        while(buffer.tell() < in_len):
            val = buffer.read(1)[0]
            
            if(val==last_val):
                repeat_len+=1
            else:
                if(repeat_len>0):
                    outstr+=str(bin(self.huff_dict[256]))[2:]
                    # handle repeat bits
                    outlen, outval = pack_bits(repeat_len)
                    repeat_len=0
                    if(outlen>3):
                        while(outlen>3):
                            outstr+='0'
                            outlen-=1
                    outstr+=str(bin(outval))[2:]
                outval = self.huff_dict[val]
                outstr+=str(bin(outval))[2:]
            last_val = val

        # 6: append EOF flag
        outstr+=str(bin(self.huff_dict[256]))[2:]
        outstr += str(bin(pack_bits(0)[1]))[2:]
        outstr += '1'

        # 7: Write out bit buffer to file
        while(len(outstr)>8):
            outval = int(outstr[:8],2)
            outstr=outstr[8:]
            compressed.extend(outval.to_bytes(1,byteorder='little'))
        
        if(len(outstr)>0):
            while len(outstr)<8:
                outstr+='0'
            outval = int(outstr[:8],2)
            compressed.extend(outval.to_bytes(1,byteorder='little'))
            outstr=outstr[8:]

        return(compressed)