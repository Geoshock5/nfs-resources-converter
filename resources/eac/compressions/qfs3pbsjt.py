from io import BufferedReader
import struct
import heapq

from resources.eac.compressions.base import BaseCompressionAlgorithm

# Credits:
# Original research / file format by AndyGura (github.com/AndyGura)
# Uncompression code rewritten by Geoshock5 based on code from AndyGura
# Compression method by Geoshock5 (github.com/Geoshock5)
# Huffman node/table code from Aashish Barnwal / GeeksForGeeks (https://www.geeksforgeeks.org/huffman-coding-greedy-algo-3/)
#
# License: CCBY-SA

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

def pack_bits(val: int):
    if val<4:
        len = 3
        return len, (val+4)
    else:
        len = (val+4).bit_length()
        return len, val+4

class Qfs3Compression(BaseCompressionAlgorithm):

    def __init__(self, *args, **kwargs):
        output_length=0
        # 0: set up necessary variables
        self.sub_ptr = 0                 # current offset in bits

        self.huff_table_codes = [0] * 16       # 0.1 - header code for each level
        self.huff_chars_per_level = [0]        # 0.2 table of code count on each level
        self.depth_table = [0]               # 0.3 table with max. val. at each level

        self.huff_dict = [None] * 257
        
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

        # File structure overview:
        # Header: similar to standard RefPack
        # 1. 0x30FB magic (2b)
        # 2. File size of the uncompressed output in bytes, big-endian (3b)
        # 3. Number of characters in the uncompressed file (1b)
        #
        # Chunk 1: Bit-packed data defining the number of leaves on each level (1 bit, 2 bit, 3 bit...)
        #
        # Chunk 2: Character table. The table in chunk 1 uses alphabetical / consecutive values for each node (0~255).  This second chunk translates that into a lookup table,
        # presumably sorted in something approximating probability order.  The last value (matching header item 3) is an escape character used to flag special functions,
        # e.g. repeat the last byte x times or do some other special things based on the following bit(s).
        # This is encoded in bitpacked form as an offset from the previous value (i.e. [val+1] = [val] + offset.  But also, if we've previously seen a value between val and val+1 in this list, don't decrement the offset.
        #
        # Chunk 3: These are the bitpacked values of the compressed file.  Read them in one-by-one and then use the dictionary from chunk 2 to write them to the output, processing any special bytes as required.
        #
        # The code below seems to work OK with TNFS files.  It's not optimal but clean enough to make the algorithm more transparent.
        # The next step will be to learn how to compress the original file back into this format, hopefully.
        # Jon T / Geoshock5 / March 2024
        
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
            val -=3 # why is this only 3 and not 4 per regular RefPack?
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

        print("256 -> " + str(self.huff_dict[256]) + ", len: " + str(len(self.huff_dict[256])))

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
        compressed.extend(in_len.to_bytes(3, byteorder='big')) # uncompressed length - TODO: write as big-endian
        compressed.append(escape_char) # special escape char for repeated bytes etc
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

        #while(len(outstr)<8):   # pack the last line to check compression works
        #    outstr += '0'

        #outval = int(outstr[:8],2)
        #compressed.extend(outval.to_bytes(1,byteorder='little'))

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
            offset-=1 # HACK: Not clear why.  Normal RefPack is offset 4.  Correcting the output to fix.
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

        # TODO: fix this
        outstr+=str(bin(self.huff_dict[256]))[2:]
        # outstr+='1' # this is incorrect.  Should be bitpacked
        outstr += str(bin(pack_bits(0)[1]))[2:]
        outstr += '1'

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
