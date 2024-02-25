from io import BufferedReader

from library.utils import read_short, read_int
from library.utils.asm_runner import AsmRunner
from resources.eac.compressions.base import BaseCompressionAlgorithm


class Qfs3Compression(BaseCompressionAlgorithm, AsmRunner):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, asm_virtual_memory_size=2 * 1024, **kwargs)
        self.output_length = 0
        self.index_table_0 = []
        self.out_value_table = [0] * 256    # values to output to uncompressed file
        self.out_len_table = [0] * 256  # length of huffman tree code. Used to know how far to step sub-pointer while reading.
        self.unk_table_3 = [0] * 16  # maybe has different size. TODO: understand how this table is built.  Created during loop 1. More important for binary pattern than decimal number.
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

        # set up counters for each operation
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
        char_count = self.al                                     # store for later. Count of characters in uncompressed file.  Not used in loop 1 but this equals the sum of values in table_110
        self.accumulate_if_needed(buffer)
        count_sing = 1
        val_shift = 15                                          # gives a max of 16 values. Decrements with every iteration of loop 1. Used to bitshift ebp when edx > 0
        count_quad = 4                                           # 0b100
        self.ebp = 0
        while True:                                         # Loop 1: this loop appears to be reading the first chunk of data 3 bytes at a time and generating tables from it.  table_110 and unk_table_3 are built.  Some values are set with set_value but not yet clear where/why
                                                            # eax = length of tree branch? ebx holds counters for inc/dec. edx holds command bits for mem storage. esi = input buffer
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
            self.eax = index_table_0_capacity = index_table_0_capacity + self.edx       # retrieve total values in the tree so far
            self.ebp += self.edx                    # add the value of the 2 bits to the base pointer
            self.ecx = 0                            # zero ecx ready for next operation
            if self.edx != 0:                       # were there any values recorded this run?
                self.cl = val_shift = val_shift & 0xFF              # Bitmask the counter
                self.eax = self.ebp << val_shift      # shift the bp by ecx bits
                self.ecx = self.eax & 0xFFFF        # drop the clamped result back into ecx.  Eventually this will overflow.
            
            val_shift = val_shift - 1                # decrement ecx shift
            self.eax = count_quad = count_quad + 4        # increment count*4 by 4 and drop it in eax
            self.ebx = count_sing = count_sing + 1               # increment count
            self.set_value('[esp+eax+550h+var_154]', self.ecx)  # what is this value? varies run-to-run. ecx is the output - where is it going?

            if self.edx == 0:                   
                continue                        # not finished yet - jump to start of loop
            if self.ecx == 0:                   # check if we've accumulated all values (255) yet
                break                           # end of reading
        # Set things up for second loop
        self.ecx = 0xFFFFFFFF
        count_loop1 = count_sing - 1                    # record the final count of steps for later.  Needed for loop 3.
        self.set_value('[esp+ebx*4+550h+var_154]', self.ecx)
        self.ebx = 0
        huff_table_0_iter_ptr = 0               # keep a count of iterations in loop 2
        self.eax = 0xFF
        unk_4 = 0
        if index_table_0_capacity > 0:                          # Loop 2
            self.index_table_0 = [None] * index_table_0_capacity
            while True:                                         # iterate through, populating index_table_0.  Count runs in ebx.  Hold total run count in ecx and break at num limit.
                counter_2 += 1
                if self.get_register_signed_value('esi') >= 0:  # if leading bit = 0
                    self.edx = 2
                    if (self.esi >> 16) == 0:
                        while True:                             # read through until we hit a 1
                            self.ebp = self.esi                 # open inbuf
                            self.edx += 1                       # count bits read (zeroes) in edx
                            self.available_acc_bits -= 1        # reduce accum counter
                            self.ebp = self.ebp >> 0x1F         # check last bit of ebp
                            self.esi = self.esi << 1            # move sub-pointer
                            self.accumulate_if_needed(buffer)
                            if self.ebp != 0:                   # stop if we hit 1
                                break
                    else:
                        while self.get_register_signed_value('esi') >= 0: # keep reading esi bit-by-bit until we hit a one
                            self.esi = self.esi << 1
                            self.edx += 1
                        self.ebx = self.edx - 1
                        self.available_acc_bits -= self.ebx
                        self.esi = self.esi << 1
                        self.accumulate_if_needed(buffer)
                    # self.ebx = self.accumulator << 8
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
                        # unk_7 = self.ecx
                        self.accumulate_if_needed(buffer)
                        # self.ecx = unk_7
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
                self.index_table_0[huff_table_0_iter_ptr] = self.al
                huff_table_0_iter_ptr += 1
                if huff_table_0_iter_ptr >= index_table_0_capacity:
                    break
        
        # Loop 3 - create output value and length table.  Write out the main characters first, then go back for the special characters from Loop 1.
        in_table_idx = 0                # The index to read from index_table_0. Increments after every read.
        self.out_len_table = [0x40] * 256 # create a 256-length table, value 64. This is the length of each output value in the input stream.
        out_table_idx = 0               # The index of out table to write to.  Increments after every value.
        huff_table_ptr_len = 1          # Unknown counter. Used in loop 3 only.  Start 1, max 8.  Length in bits of the Huffman table pointer being read from source in loop 4.
        unk_9 = 0                       # used in loop 3 - and loop 4 in case of break flag (0x60h / 96d)
        
        # Run outer and inner loops to populate output tables.  Expected output:
        #   - Write the increasing value of index_table_0 to out_value_table, table_110[idx] value times.
        #   - Write the increasing value of huff_table_ptr_len to out_len_table, table_110[idx] value times.
        if count_loop1 >= 1:               # should be > 1 as this value came from loop 1
            val_shift = 4
            count_loop3 = 7
            while True:                                         # populate out_value_table and out_len_table
                counter_3 += 1                                  # TODO: remove, unnecessary.  Used only for debug at this time.
                assert val_shift % 4 == 0
                count_loop3_inner = table_110[int(val_shift / 4)]         # read back in from table_110
                if huff_table_ptr_len >= 9:                                  # run max 8 times
                    break
                break_outer = False
                continue_outer = False
                while True:
                    count_loop3_inner -= 1
                    if count_loop3_inner < 0:           # when we hit zero
                        val_shift = val_shift + 4       # move to next value in table_110
                        count_loop3 -= 1                # decrement outer loop counter
                        huff_table_ptr_len += 1
                        if huff_table_ptr_len <= val_shift:
                            continue_outer = True
                            break
                        else:
                            break_outer = True
                            break
                    else:
                        index_table_0_value = self.index_table_0[in_table_idx]
                        in_table_idx += 1
                        out_len_val = huff_table_ptr_len
                        if char_count == index_table_0_value:
                            unk_9 = huff_table_ptr_len
                            out_len_val = 0x60
                        if ((1 << count_loop3) <=0):
                            continue
                        for x in range(0,(1 << count_loop3)):
                            self.out_value_table[out_table_idx] = index_table_0_value
                            self.out_len_table[out_table_idx] = out_len_val

                            out_table_idx += 1
                if break_outer:
                    break
                if continue_outer:
                    continue
                break

        # Loop 4: write the file out using the lookup tables
        while True:                                             
            counter_4 += 1
            if len(uncompressed) > self.output_length:
                raise Exception('Uncompress algorythm writes more that file length')
            self.eax = self.out_len_table[self.esi >> 0x18]       # should be little-endian - clamp index to 256. Len to read from inbuf / esi
            self.available_acc_bits -= self.eax
            while self.available_acc_bits >= 0:
                self.edx = self.esi >> 0x18                     # get 1 byte from front of stream
                for _ in range(4):
                    self.append_to_output(uncompressed, self.out_value_table[self.edx])   # read value of out_value_table from instream lookup and output
                    self.esi = self.esi << self.al
                    self.edx = self.esi >> 0x18
                    self.eax = self.out_len_table[self.edx]   # length of value that was read from input
                    self.available_acc_bits -= self.eax
                    if self.available_acc_bits < 0:
                        break                                   # stop and get more input
            self.available_acc_bits += 0x10                     # add 2 bytes to available acc. bits tally
            if self.available_acc_bits >= 0:
                self.append_to_output(uncompressed, self.out_value_table[(self.esi >> 0x18)]) # read value of out_value_table from instream lookup
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
            if self.al != char_count:
                if self.available_acc_bits >= 0:
                    self.append_to_output(uncompressed, self.al)
                    continue
            self.accumulate_if_needed(buffer)
            if self.al != char_count:
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
        print("C1=" + str(count_sing - 1))
        print("C2=" + str(counter_2))
        print("C3=" + str(counter_3))
        print("C4=" + str(counter_4))
        return uncompressed
