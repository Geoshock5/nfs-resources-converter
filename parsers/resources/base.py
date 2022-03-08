from abc import ABC, abstractmethod
from io import BufferedReader
from typing import List


class BaseResource(ABC):

    resources: List['BaseResource'] = []

    name = None
    parent = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.resources = []

    @abstractmethod
    def read(self, buffer: BufferedReader, length: int, path: str = None) -> int:
        # returns how many bytes were used
        pass

    @abstractmethod
    def save_converted(self, path: str):
        pass
