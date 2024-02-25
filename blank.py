from io import BufferedReader, IOBase
from resources.eac.compressions.qfs3pbs import Qfs3Compression as Qfs

inbuf = open("TSUPRA.PBS","rb")
outbuf = open("decoded.pbs","wb")

uncom = Qfs()

outbuf.write(uncom.uncompress(inbuf, inbuf.__sizeof__))
outbuf.close()