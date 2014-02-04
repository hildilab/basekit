from __future__ import with_statement
from __future__ import division

import struct
import collections


MrcHeader = collections.namedtuple( "MrcHeader", [
    "nx", "ny", "nz", "mode", 
    "nxstart", "nystart", "nzstart",
    "mx", "my", "mz",
    "xlen", "ylen", "zlen",
    "alpha", "beta", "gamma",
    "mapc", "mapr", "maps",
    "amin", "amax", "amean",
    "ispg", "next", "createid", 
    "extra_data1", "nint", "nreal", 
    "extra_data2", "imod_stamp", "imod_flags",
    "idtype", "lens", "nd1", "nd2", "vd1", "vd2",
    "current_tilt_x", "current_tilt_y", "current_tilt_z",
    "original_tilt_x", "original_tilt_y", "original_tilt_z",
    "xorg", "yorg", "zorg",
    "cmap", "stamp", "rms", "nlabl",
    "label1", "label2", "label3", "label4", "label5", 
    "label6", "label7", "label8", "label9", "label10"
])


def get_mrc_header( mrc_file ):
    """ http://bio3d.colorado.edu/imod/doc/mrc_format.txt
    """
    byteorder = {
        0x1111: '>',    # big,
        0x4144: '<'     # little
    }
    with open( mrc_file, "rb" ) as fp:
        header = fp.read( 1024 )
    stamp = struct.unpack( "i", header[212:216] )[0]
    if header[208:212]!="MAP ":
        raise Exception("can only read new style mrc files")
    if stamp not in byteorder:
        raise Exception("could not deduce byteorder")
    h = struct.unpack(
        (
            byteorder.get( stamp )+
            "3i"    # number of cols, rows, sections
            "i"     # mode
            "3i"    # xyz start
            "3i"    # grid size
            "3f"    # cell size
            "3f"    # cell angles
            "3i"    # map col row section
            "3f"    # min max mean
            "2i"    # space group, no bytes in ext header
            "h"     # create id
            "30s"   # extra data
            "2h"    # nint, nreal
            "20s"   # extra data
            "2i"    # imodStamp, imodFlag
            "6h"    # idtype, lens, nd1, nd2, vd1, vd2
            "6f"    # tiltangles
            "3f"    # origin of image
            "4s"    # "MAP "
            "i"     # First two bytes have 
                    # 17 and 17 for big-endian or 
                    # 68 and 65 for little-endian
            "f"     # RMS deviation of densities from mean density
            "i"     # Number of labels with useful data
            +("80s"*10)  # 10 labels of 80 charactors
        ),
        header
    )
    return MrcHeader._make( h )


def getMrc( mrc_file, param ):
    header = get_mrc_header( mrc_file )
    for name, value in zip(header._fields, header):
        if name==param:
            return value



