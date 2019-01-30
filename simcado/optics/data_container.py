from astropy.table import Table
from astropy.io import ascii as ioascii
from astropy.io import fits

from .. import utils


class DataContainer:
    def __init__(self, filename=None, **kwargs):

        filename = utils.find_file(filename)
        self.meta = {"filename" : filename}
        self.meta.update(kwargs)

        self.headers = []
        self.table = None
        self._file = None

        if filename is not None:
            if self.is_fits:
                self._load_fits()
            else:
                self._load_ascii()

    def _load_ascii(self):
        self.table = ioascii.read(self.meta["filename"])
        hdr_dict = utils.convert_table_comments_to_dict(self.table)
        if isinstance(hdr_dict, dict):
            self.headers += [hdr_dict]
        else:
            self.headers += [None]

        self.meta.update(hdr_dict)

    def _load_fits(self):
        self._file = fits.open(self.meta["filename"])
        for ext in self._file:
            self.headers += [ext.header]

        self.meta.update(dict(self._file[0].header))

    def get_data(self, ext=0, layer=None):
        data_set = None
        if self.is_fits:
            if isinstance(self._file[ext], fits.BinTableHDU):
                data_set = Table.read(self._file[ext], format="fits")
            else:
                if self._file[ext].data is not None:
                    data_dims = len(self._file[ext].data.shape)
                    if data_dims == 3 and layer is not None:
                        data_set = self._file[ext].data[layer]
                    else:
                        data_set = self._file[ext].data
        else:
            data_set = self.table

        return data_set

    @property
    def is_fits(self):
        flag = False
        if self.meta["filename"] is not None:
            if self.meta["filename"].split(".")[-1].lower() in "fits":
                flag = True

        return flag
