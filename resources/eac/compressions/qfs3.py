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
# Huffman node/table code from Aashish Barnwal / GeeksForGeeks (https://www.geeksforgeeks.org/huffman-coding-greedy-algo-3/)
#
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
        self.unk_table_3 = [0] * 256  # maybe has different size. TODO: understand how this table is built.  Created during loop 1. More important for binary pattern than decimal number.
        self.accumulator = 0
        self.available_acc_bits = 0

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

        table_110 = [0] * 16                                        # this is the value of the edx in each part of chunk 1

        self.define_variable('var_14C', -0x14C, 4)
        self.define_variable('var_154', -0x154, 4)

        index_table_0_capacity = 0                                  # this is the sum of all the values in table_110
        unk_counter = 0

        # set up counters for each operation
        counter_1 = 0
        counter_2 = 0
        counter_3 = 0
        counter_4 = 0

        self.available_acc_bits = 0
        file_header = self.accumulator = read_short(buffer, 'big') # 2b read - read the header
        self.read_next(buffer)                                     # 2b read - read the output size (little-endian)
        self.esi = self.accumulator << 16                          # insert stream into esi, minus header
        # if compressed size presented                             -- stock TNFS PBS are not compressed (0x30FB) so this will not happen
        if file_header & 0x100:                                    
            self.edx = read_int(buffer, 'big')
            self.available_acc_bits = 8
            self.esi = self.edx << 8
            self.accumulator = self.edx
            # reset size presented flag bit
            file_header = file_header & 0xFEFF
        self.ebx = self.esi >> 0x18                                # ebx will be little-endian value of size = 000D
        self.available_acc_bits -= 8                               # cancel compressed value if compressed
        self.esi = self.esi << 8                                   # step through 1 byte
        self.accumulate_if_needed(buffer)                          # this will read unless compressed - so for TNFS, yes
        self.eax = self.esi >> 16                                  # esi should be blank so = 0000
        self.available_acc_bits -= 16                              # should be 16 by now, so = 0
        self.esi = self.esi << 16                                  # esi should be blank so = 0000
        self.output_length = self.eax                              # Set output length.  eax = first 3 bytes = length header in file
        self.accumulate_if_needed(buffer)
        self.edx = self.output_length = self.output_length | (self.ebx << 16) # size is a 3-byte value. Use bitmask to adjust.
        self.ebx = 1
        self.eax = self.esi >> 0x18
        self.available_acc_bits -= 8
        self.esi = self.esi << 8
        unk_1 = self.al                                     # store for later.  Not used in loop 1 but this equals the sum of values in table_110
        self.accumulate_if_needed(buffer)
        count_sing = 1
        val_shift = 15                                          # gives a max of 16 values. Decrements with every iteration of loop 1. Used to bitshift ebp when edx > 0
        count_quad = 4                                           # 0b100
        self.ebp = 0
        while True:                                         # Loop 1: this loop appears to be reading the first chunk of data 3 bytes at a time and generating tables from it.  table_110 and unk_table_3 are built.  Some values are set with set_value but not yet clear where/why
                                                            # eax = length of tree branch? ebx holds counters for inc/dec. edx holds command bits for mem storage. esi = input buffer
            counter_1 += 1
            self.ebp = self.ebp << 1                        # TODO: evaluate evolution of ebp across loops
            self.ecx = index_table_0_capacity               # set to the sum of all flag bits so far
            self.eax = self.ebp - self.ecx                  # TODO: work out what this value represents. Could be some kind of mask?
            self.edx = count_quad                                # unk_2 is always a multiple of 4 as it increments by 4 each run
            assert self.edx % 4 == 0                        # logic below will work in different way if not
            self.unk_table_3[int(self.edx / 4)] = self.eax  # store the value above (bp - ecx) in a lookup table. How many we've had vs could have had?
            if self.get_register_signed_value('esi') >= 0:  # only switch if the leading bit is 0
                self.eax = 2                                # start eax at 2
                if (self.esi >> 16) == 0:                   # check we aren't low on data
                    self.edx = 0
                    while self.edx == 0:                    # keep reading until we hit a 1 in our input buffer
                        self.edx = self.esi >> 0x1F
                        self.eax += 1                       # count how many zeroes, minimum 2 (eax = 2 initially)
                        self.available_acc_bits -= 1        # move sub-pointer up 1 bit
                        self.esi = self.esi << 1
                        self.accumulate_if_needed(buffer)
                else:
                    while self.get_register_signed_value('esi') >= 0:   # move 1 bit at a time and count the zeroes until we hit a 1
                        self.esi = self.esi << 1            # move through inbuf 1 bit at a time until we hit 1
                        self.eax += 1                       # add 1 to the count for every zero
                    self.edx = self.eax - 1                 # total zeroes read, minus 1(?)
                    self.available_acc_bits -= self.edx      # update how much we read
                    self.esi = self.esi << 1                # move up 1 bit
                    self.accumulate_if_needed(buffer)
                if self.get_register_signed_value('eax') <= 16: # checking against 2 + number of zeroes in a row.  Why does this need to be signed?
                    self.edx = self.esi >> (0x20 - self.eax)    # esi shifted by 32 minus the value.
                    self.ecx = self.al                          # 
                    self.available_acc_bits -= self.eax         # subtract what we've just processed
                    self.esi = self.esi << self.cl              # clear the data we just processed from the input stream
                    self.accumulate_if_needed(buffer)           # read more if we're low
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
                self.edx += (1 << self.al)          # table value = 2^eax. Append it to edx.
            else:                                   # for values where the if is <=0
                self.edx = self.esi >> 0x1D         # read the top 3 bits of the buffer (shift 0x1D = 29 bits)
                self.available_acc_bits -= 3
                self.esi = self.esi << 3            # advance the buffer 3 bits
                self.accumulate_if_needed(buffer)
            self.edx -= 4                           # Arriving from the normal loops, this will subtract. Arriving from the else, this would cancel bit 0 of the 3-bit buffer.
            self.eax = count_quad
            assert self.eax % 4 == 0                # should always be true as unk_2 increments by 4 per run
            table_110[int(self.eax / 4)] = self.edx # holding the value of edx from this run
            self.eax = index_table_0_capacity       # retrieve total values in the tree so far
            self.ebp += self.edx                    # add the value of the 2 bits to the base pointer
            self.eax += self.edx                    # add the 2 bits to eax
            self.ecx = 0                            # zero ecx ready for next operation
            index_table_0_capacity = self.eax       # add the total values counted from this run to the total
            if self.edx != 0:                       # were there any values recorded this run?
                self.cl = val_shift & 0xFF              # Bitmask the counter
                self.eax = self.ebp << self.cl      # shift the bp by ecx bits
                self.ecx = self.eax & 0xFFFF        # drop the clamped result back into ecx.  Eventually this will overflow.
            self.ebx = val_shift - 1                # decrement unk_3 but drop it in ebx
            self.eax = count_quad = count_quad + 4        # increment unk_2 by 4 and drop it in eax
            val_shift = self.ebx                    # formally decrement unk_3.  Is this a counter of e.g. number of chunks / instructions?
            self.ebx = count_sing + 1               # increment unk_90 but drop it in ebx
            self.set_value('[esp+eax+550h+var_154]', self.ecx)  # what is this value? varies run-to-run. ecx is the output - where is it going?
            count_sing = self.ebx                   # formally increment unk_90
            if self.edx == 0:                   
                continue                        # not finished yet - jump to start of loop
            if self.ecx == 0:                   # check if we've accumulated all values (255) yet
                break                           # end of reading
        # Set things up for second loop
        self.ecx = 0xFFFFFFFF
        unk_8 = self.ebx - 1                    # record the final count of steps for later.  Needed for loop 3.
        self.set_value('[esp+ebx*4+550h+var_154]', self.ecx)
        self.ebx = 0
        huff_table_0_iter_ptr = 0
        self.eax = 0xFF
        unk_4 = 0
        if index_table_0_capacity > 0:                          # Loop 2
            self.index_table_0 = [None] * index_table_0_capacity
            while True:                                         # iterate through, populating index_table_0.  Count runs in ebx.  Hold total run count in ecx and break at num limit.
                counter_2 += 1
                if self.get_register_signed_value('esi') >= 0:  # if leading bit = 0
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
                        while self.get_register_signed_value('esi') >= 0: # keep reading esi bit-by-bit until we hit a one
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
        self.unk_table_2 = [0x40] * 256 # create a 256-length table, value 64
        self.edx = 0
        self.ebx = 0
        self.ecx = unk_8
        unk_5 = 1
        unk_9 = 0                       # used in loop 3 - and loop 4 in case of break flag (0x60h / 96d)
        if self.ecx >= 1:               # should be > 1 as this value came from loop 1
            self.ecx = 7
            self.eax = 4
            count_quad = self.ecx
            val_shift = self.eax
            while True:                                         # populate unk_table and unk_table_2
                counter_3 += 1
                self.eax = val_shift
                assert self.eax % 4 == 0
                self.eax = table_110[int(self.eax / 4)]         # read back in from table_110
                self.ebp = unk_5
                unk_6 = self.eax
                if self.get_register_signed_value('ebp') >= 9:
                    break
                self.ecx = count_quad
                self.ebp = 1 << self.cl
                break_outer = False
                continue_outer = False
                while True:
                    self.eax = unk_6 - 1
                    unk_6 = self.eax
                    if self.eax == 0xFFFFFFFF:
                        self.ebp = val_shift = val_shift + 4
                        self.eax = count_quad = count_quad - 1
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
        while True:                                             # write the file out using the lookup tables
            counter_4 += 1
            if len(uncompressed) > self.output_length:
                raise Exception('Uncompress algorythm writes more that file length')
            self.eax = self.unk_table_2[self.esi >> 0x18]       # should be little-endian - clamp index to 256. Len to read from inbuf / esi
            self.available_acc_bits -= self.eax
            while self.available_acc_bits >= 0:
                self.edx = self.esi >> 0x18                     # get 1 byte from front of stream
                for _ in range(4):
                    self.append_to_output(uncompressed, self.unk_table[self.edx])   # read value of unk_table from instream lookup and output
                    self.esi = self.esi << self.al
                    self.edx = self.esi >> 0x18
                    self.eax = self.unk_table_2[self.edx]   # length of value that was read from input
                    self.available_acc_bits -= self.eax
                    if self.available_acc_bits < 0:
                        break                                   # stop and get more input
            self.available_acc_bits += 0x10                     # add 2 bytes to available acc. bits tally
            if self.available_acc_bits >= 0:
                self.append_to_output(uncompressed, self.unk_table[(self.esi >> 0x18)]) # read value of unk_table from instream lookup
                self.read_next(buffer)                          # read 2 bytes
                self.esi = self.accumulator << (0x10 - self.available_acc_bits) # advance by 16 minus available bits
                continue                                        # skip the rest and restart the loop
            self.available_acc_bits += self.eax - 0x10
            if self.eax == 0x60:                                # 96 = specific break code?
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
        # This was part of the original algorithm but is not used here as it makes the decompression much more complex to understand.
        # Computationally probably more expensive but clock cycles are cheap these days vs the 90s ;-)

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
                            print("Function not implemented yet or end of file flag reached!")
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

        checksum = 0    # Refer to car physics documentation for info.  Not used in compression algorithm but the game will break the physics if this doesn't match.

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
            
            if(buffer.tell()<=1880):
                checksum += val # calculate for editing physics

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

        # print ("Checksum value: " + str(hex(checksum))) # DEBUG: Since we're reading it anyway, print sum of bytes 0-1880 for modding purposes.

        return(compressed)