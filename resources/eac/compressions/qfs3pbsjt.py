from io import BufferedReader
import struct

from resources.eac.compressions.base import BaseCompressionAlgorithm

class Qfs3Compression(BaseCompressionAlgorithm):

    def __init__(self, *args, **kwargs):
        output_length=0
        # 0: set up necessary variables
        self.sub_ptr = 0                 # current offset in bits

        self.huff_table_codes = [0] * 16       # 0.1 - header code for each level
        self.huff_chars_per_level = [0]        # 0.2 table of code count on each level
        self.depth_table = [0]               # 0.3 table with max. val. at each level
        
    def append_to_output(self, buffer, value):
        buffer.extend(value.to_bytes(1, 'little'))

    def reuse_output_byte(self, buffer, length):
        value = buffer[-1]
        for i in range(length):
            self.append_to_output(buffer, value)

    def refill_buf(self, buf, buffer, sub_ptr):
        while (sub_ptr >= 8):
            sub_ptr -= 8
            buf = buf | int.from_bytes(buffer.read(1), "big") << sub_ptr
            buf = buf & 0xFFFFFFFF
        return sub_ptr, buf

    def count_bits(self, buf, buffer, sub_ptr):
        if(buf >> 31 == 0):
            len_val = 2
            while True:
                # count the zeroes
                len_val += 1
                buf = (buf << 1) & 0xFFFFFFFF
                sub_ptr += 1
                sub_ptr, buf = self.refill_buf(buf, buffer, sub_ptr) #TODO: collapse these into a single function
                if(buf>>31 == 1):
                    #sub_ptr -= 1
                    break
            buf = (buf << 1) & 0xFFFFFFFF
            sub_ptr += 1 #len_val
            val = buf >> (32-len_val)
            val += 1 << len_val
            buf = (buf << len_val) & 0xFFFFFFFF
            sub_ptr += len_val
            sub_ptr, buf = self.refill_buf(buf, buffer, sub_ptr)
        else:
            val = buf >> 29
            sub_ptr += 3
            buf = (buf << 3) & 0xFFFFFFFF
            sub_ptr, buf = self.refill_buf(buf, buffer, sub_ptr)
        return val, buf, sub_ptr

    def uncompress(self, buffer: BufferedReader, input_length: int) -> bytes:
        uncompressed: bytearray = bytearray()

        # 1: read the file header
        file_header = struct.unpack(">h", buffer.read(2))[0]
        out_size = int.from_bytes(struct.unpack("<BBB", buffer.read(3)),"big")
        char_count = struct.unpack("B", buffer.read(1))[0]

        # 2: read in the Huff tree from the header
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
        print("Tree complete.  Total values: " + str(value_count) + ", max. length: " + str(length))
        length -= 1

        # 3: read in the character list
        char_table = [None] * value_count
        current_char=0
        char_value = 0xFF

        while True:
            val, buf, self.sub_ptr = self.count_bits(buf, buffer, self.sub_ptr)
            val -=3
            while (val!=0):
                char_value = (char_value + 1) & 0xFF
                if char_value not in char_table:
                    val -=1
            char_table[current_char] = char_value
            current_char += 1
            if (current_char >= len(char_table)):
                break

        # 4: generate the lookup tables for output as needed
        # This was part of the original algorithm but is not used here.
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
                print("Reached compression length")
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
                            print("Error: Not implemented yet or end of file!")
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
        byte_count = [0] * 256
        sorted_bytes = [0] * 256

        buffer.seek(0,2)
        in_len = buffer.tell()
        buffer.seek(0,0)

        # 1 count the bits in the source file
        for _ in range(int(in_len)):
            val = buffer.read(1)[0]
            byte_count[val] += 1

        # 2 sort by count and then by number
        sorted_bytes = byte_count
        sorted_bytes.sort(reverse=1)

        # 3 generate the tree
        current_val = 0
        length = 1
        val_length = 1
        huff_key = 0

        while True:
            vals_rem_on_level = (2**length) - 1
            val_limit = (2**(val_length-1))
            while((current_val<val_limit) & (vals_rem_on_level>0)):
                huff_key += 1
                vals_rem_on_level -= 1
                # move to next value
                current_val+=1
            if(vals_rem_on_level==0):
                length +=1
            if(current_val>=val_limit):
                length +=1
                val_length += 1
            huff_key << 1
            continue

        # 4 generate dictionary of Huff values

        # 5 write the output file
        # 5.1 header
        # 5.2 tree depth
        # 5.3 alphabet
        # 5.4 encoded output

        return(compressed)