# automatically generated by the FlatBuffers compiler, do not modify

# namespace: NetEncoding

import flatbuffers


class JSONEncodedArray(object):
    __slots__ = ["_tab"]

    @classmethod
    def GetRootAsJSONEncodedArray(cls, buf, offset):
        n = flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, offset)
        x = JSONEncodedArray()
        x.Init(buf, n + offset)
        return x

    # JSONEncodedArray
    def Init(self, buf, pos):
        self._tab = flatbuffers.table.Table(buf, pos)

    # JSONEncodedArray
    def Data(self, j):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            a = self._tab.Vector(o)
            return self._tab.Get(
                flatbuffers.number_types.Uint8Flags, a + flatbuffers.number_types.UOffsetTFlags.py_type(j * 1)
            )
        return 0

    # JSONEncodedArray
    def DataAsNumpy(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            return self._tab.GetVectorAsNumpy(flatbuffers.number_types.Uint8Flags, o)
        return 0

    # JSONEncodedArray
    def DataLength(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            return self._tab.VectorLen(o)
        return 0


def JSONEncodedArrayStart(builder):
    builder.StartObject(1)


def JSONEncodedArrayAddData(builder, data):
    builder.PrependUOffsetTRelativeSlot(0, flatbuffers.number_types.UOffsetTFlags.py_type(data), 0)


def JSONEncodedArrayStartDataVector(builder, numElems):
    return builder.StartVector(1, numElems, 1)


def JSONEncodedArrayEnd(builder):
    return builder.EndObject()
