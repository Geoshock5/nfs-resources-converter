from io import BufferedReader


def read_utf_bytes(buffer: BufferedReader, len: int) -> str:
    return buffer.read(len).decode('utf-8')


def read_int(buffer: BufferedReader, byteorder='little') -> int:
    return int.from_bytes(buffer.read(4), byteorder=byteorder)


def read_3int(buffer: BufferedReader, byteorder='little') -> int:
    return int.from_bytes(buffer.read(3), byteorder=byteorder)


def read_signed_int(buffer: BufferedReader) -> int:
    return int.from_bytes(buffer.read(4), byteorder='little', signed=True)


def read_short(buffer: BufferedReader, byteorder='little') -> int:
    value = buffer.read(2)
    value = value.ljust(2, b'\0')
    return int.from_bytes(value, byteorder=byteorder)

def read_signed_short(buffer: BufferedReader) -> int:
    return int.from_bytes(buffer.read(2), byteorder='little', signed=True)


def read_byte(buffer: BufferedReader) -> int:
    return int.from_bytes(buffer.read(1), byteorder='little')


def read_signed_byte(buffer: BufferedReader) -> int:
    return int.from_bytes(buffer.read(1), byteorder='little', signed=True)


def read_vector3_as_list(buffer, size_koeff = 1):
    coords = [read_signed_short(buffer) * size_koeff for _ in range(3)]
    return [coords[0], coords[2], coords[1]]
