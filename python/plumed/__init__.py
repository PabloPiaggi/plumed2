import ctypes as ct
from ctypes.util import find_library
import sys

import numpy

## Load the plumed libary. This is machine independent, so
## it will use the right extension depending on the architecture
## (.dylib on mac, .so on linux, ...)
## Remember to source the 'sourceme.sh' file first.
## The previous code was not working on linux. Until it is solved
## in a cleaner fashion, I will leave linux as a special case.
if sys.platform=="linux" or sys.platform=="linux2":
	_libplumed = ct.CDLL("libplumed.so")
else:
	libraryname = find_library('libplumed')
	if libraryname is None:
    		raise ImportError("Unable to find the PLUMED library.")
	_libplumed = ct.CDLL(libraryname)

class PlumedError(Exception):
    """
    Raised when a Plumed-related error happens.
    """
    pass

class Plumed(object):
    def __init__(self,start=False):
        """
        Init an empty Plumed class. For efficiency, don't start the
        plumed environment, by default.
        :param start: if True, also start the plumed environment
        """
        self._p = None
        # The type of integer in the compiled version of plumed
        # will be a numpy dtype, either int32 or int64
        # It is stored the first time start_plumed is called, and
        # it is NOT cleared at stop_plumed (for efficiency reasons,
        # assuming the library is not recompiled with different int
        # size while the python script is running...)
        # self._int_type is the numpy integer type
        self._int_type = None
        # self._c_int_type is the python ctypes integer type
        self._c_int_type = None
        if start:
            self.start_plumed()
        # Will cache here values to avoid garbage collection
        # At the moment used only by self._cmd_old
        self._values_cache = []

    def __del__(self):
        """
        When the Python object is deleted or goes out of scope, stop the
        plumed environment if the user did not do it already.
        """
        if self._p is not None:
            self.stop_plumed()

    def start_plumed(self):
        """
        Start the plumed environment. 
        :raise PlumedError: if the Plumed environment has already been started.
        """
        if self._p is None:
            # Call the plumed_create() C function, but 
            # casting the result as a void pointer (important, on some
            # platforms it would cut the 64-bit address to 32 bits only
            f = _libplumed.plumed_create
            f.restype = ct.c_void_p
            self._p = f()

            # Here I detect the integer size (# of bytes)
            if self._int_type is None:
                # I know that Plumed returns this as a unsigned 8bit integer
                # (to be machine-independent: at this point, we still do not
                # know the integer size!!)
                int_prec_c = ct.c_uint8()
                # Manual call to plumed_cmd to avoid the logic inside self.cmd()
                _libplumed.plumed_cmd(ct.c_void_p(self._p), ct.c_char_p(bytes(
                    "getIntegerPrecision".encode(encoding='UTF-8',errors='strict'))), 
                    ct.byref(int_prec_c))
                if int_prec_c.value == 4:
                    self._int_type = numpy.int32
                    self._c_int_type = ct.c_int32
                elif int_prec_c.value == 8:
                    self._int_type = numpy.int64
                    self._c_int_type = ct.c_int64
                else:
                    raise PlumedError("Internal errror: unrecognized value of the integer "
                                      "size ({})".format(int_prec_c.value))
        else:
            raise PlumedError("The Plumed environment has already been started!")

    def stop_plumed(self):
        """
        Stop the plumed environment. 
        :raise PlumedError: if the Plumed environment has been already stopped,
           or was never started.
        """
        if self._p is not None:
            # Must cast self._p to a void pointer
            _libplumed.plumed_finalize(ct.c_void_p(self._p))
            self._p = None
            # Remove cache to empty memory (garbage-collector will take care of this)
            self._values_cache = []
        else:
            raise PlumedError("The Plumed environment has not been started!")

    def cmd(self, key, val=None):
        """
        Call the Plumed 'cmd' command with a key and a value.
        :param key: a string with the command name
        :param value: a value with the correct type expected by Plumed. Type is checked
          for conversion to the correct C type, but no check is done on the type expected
          by Plumed.
          Only int, float, strings, and numpy arrays of floats or ints are currently supported.
        :raise ValueError: if the value is not among the expected ones
        :return: None for the time being.
        .. note:: Some technical notes:
           
          - it is not needed to flatten the array. But if you do, you need to cache
            the flattened array in a class method (see the implementation of _cmd_old)
            to avoid that the garbage collector of python thinks the variable is out of
            scope and reuses the memory before Plumed can read the value
          - it is important to call, instead, ``numpy.ascontiguousarray`` to have the
            array contiguous in memory and in C order (important e.g. if the user
            passes ``arr[::2]`` or ``arr[:,1]``, for instance).
          - It is not correct to pass arr.ctypes.data, see
            http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.ndarray.ctypes.html
            Use instead ``arr.ctypes.data_as(ct.POINTER(ct.c_double))``
            Actually here I directly use ``arr.ctypes.data_as(ct.c_void_p)``
            since Plumed accepts void pointers and I can reduce the number of type checks in
            python, in this way.
        .. todo:: Things to check:
           - tuples and lists are not accepted for the time being. If needed to simplify
             the user interface, convert them transparently to arrays with something like::
             
               if (isinstance(val, (list,tuple)):
                   internal_val = numpy.array(val)

           - Should probably also check the size of floats/doubles (assumed to be 64 bits here)
           - I'm not sure that the way we are treating now arrays with ascontiguousarray,
             and strings by creating bytes(val.encode(...)), are safe from garbage collection.
             If they are not, we need to cache them.
        """
        if self._p is None:
            raise PlumedError("The plumed environment has not been started yet. "
                             "Call the start_plumed() method first.")

        if self._int_type is None:
            raise PlumedError("Error! The size of integer is not set! This is an internal error.")
        
        # Handle different data types
        # I try to minimize the number of checks as they are quite expensive in python
        if val is None:
            value = None
        elif isinstance(val, (int, long)):
            # Cast to the correct int type
            value=ct.byref(self._c_int_type(self._int_type(val)))
        elif isinstance(val, float):
            value=ct.byref(ct.c_double(val))
        elif isinstance(val, basestring):
            value=ct.c_char_p(bytes(val.encode(encoding='UTF-8',errors='strict')))
        elif isinstance(val, numpy.ndarray):
            if numpy.issubdtype(val.dtype, float):
                value = numpy.ascontiguousarray(val).ctypes.data_as(ct.c_void_p)
            elif numpy.issubdtype(val.dtype, int):
                # If it's integer, cast it first to the correct type
                value = numpy.ascontiguousarray(val.astype(self._int_type)).ctypes.data_as(
                    ct.c_void_p)
            else:
                raise ValueError("Unknown array type ({})".format(val.dtype.name))
        else:
            raise ValueError("Unknown value type ({})".format(str(type(val))))

        _libplumed.plumed_cmd(ct.c_void_p(self._p), ct.c_char_p(bytes(
                    key.encode(encoding='UTF-8',errors='strict'))), value)

    def grab(self, key):
        # API assumption: data is always float
        
        ndims_c = self._c_int_type()
        shape = numpy.zeros((10,), dtype=self._int_type) # maximum # of dimensions is 10
        shape_c = shape.ctypes.data_as(ct.POINTER(self._c_int_type))
        
        _libplumed.plumed_grabdimension(ct.c_void_p(self._p), ct.c_char_p(bytes(
                    key.encode(encoding='UTF-8',errors='strict'))),
                    ct.byref(ndims_c), shape_c)

        ndims = ndims_c.value
        # print "NDIMS:", ndims
        # print "SHAPE_FULL:", shape
        # print "SHAPE:", shape[:ndims]

        if ndims == 0:
            # Single variable
            value_c = ct.c_float()
            
            _libplumed.plumed_grabdata(ct.c_void_p(self._p), ct.c_char_p(bytes(
                key.encode(encoding='UTF-8',errors='strict'))),
                ct.byref(value_c))

            return value_c.value
        else:
            # Array
            value = numpy.zeros(shape[:ndims], dtype=float)
            value_c = value.ctypes.data_as(ct.c_void_p)
            _libplumed.plumed_grabdata(ct.c_void_p(self._p), ct.c_char_p(bytes(
                key.encode(encoding='UTF-8',errors='strict'))),
                value_c)
            return value
            
        
        
    def _cmd_old(self, key, val=None):       
        """
        Some old implementation, do not use. To be removed eventually.
        Left here to give a look to some issues about caching and avoiding the
        garbage-collector to delete the arrays you are passing.
        """

        ### FLAGS TO CHANGE THE BEHAVIOR

        # If True, flattens arrays, and stores them inside self._values_cache to avoid
        # garbage-collection. Actually, this is better if you don't want
        # the array to be changed by Plumed.
        # Both behaviors should work properly and avoid garbage collection
        # with_caching makes copies so maybe can be avoided
        with_caching = False
        # Simpler casting of arrays
        simpler_casting = True

        if self._p is None:
            raise ValueError("Plumed was not initialized")
	# Handle different data types
        if val is None:
            value = None
        elif isinstance(val, (int, long)):
            value=ct.byref(self._c_int_type(val))
        elif isinstance(val, float):
            value=ct.byref(ct.c_double(val))
        elif isinstance(val, basestring):
            value=ct.c_char_p(val)
        elif isinstance(val, numpy.ndarray):
            if simpler_casting:
                if val.dtype not in [float, int]:
                    raise ValueError("Unknown array type ({})".format(str(flval.dtype)))
                # ascontiguous array is important! e.g. if you pass arr[::2]
                # not sure if this is enough to avoid garbage-collection though...
                value = numpy.ascontiguousarray(val).ctypes.data_as(ct.c_void_p)
            else:
                if with_caching:
                    #flval = numpy.ascontiguousarray(val.flatten().copy())
                    flval = numpy.ascontiguousarray(val).flatten()
                    # Important! Otherwise this will be garbage-collected and 
                    # weird numbers will be received by the C code!!!
                    self._values_cache.append(flval)
                    if flval.dtype == float:
                        # NOT THE RIGHT WAY, SEE http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.ndarray.ctypes.html
                        #value=flval.ctypes.data
                        value = flval.ctypes.data_as(ct.POINTER(ct.c_double))
                    elif val.dtype == int:
                        ## WARNING! PROBABLY, IF IN C++ INT ARE INT32, WE SHOULD CONVERT THIS TO INT32!!
                        ## WITH .astype(np.int32)
                        value = flval.ctypes.data_as(ct.POINTER(self._c_int_type))
                    else:
                        raise ValueError("Unknown array type ({})".format(str(flval.dtype)))
                    if key == 'setPositions':
                        debug_print = True
                else:
                    if val.dtype.name == "float64":
                        # NOT THE RIGHT WAY, SEE http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.ndarray.ctypes.html
                        #value=flval.ctypes.data
                        value = numpy.ascontiguousarray(val).ctypes.data_as(ct.POINTER(ct.c_double))
                    elif val.dtype == int:
                        value = numpy.ascontiguousarray(val).ctypes.data_as(ct.POINTER(self._c_int_type))
                    else:
                        raise ValueError("Unknown array type ({})".format(val.dtype.name))
        else:
            raise ValueError("Unknown value type ({})".format(str(type(val))))
        # Add the possibility to input tuples or lists.

        _libplumed.plumed_cmd(ct.c_void_p(self._p), ct.c_char_p(key), value)
        #if debug_print:
        #    print 'A', key, numpy.ctypeslib.as_array(value,shape=(flval.size,))     
