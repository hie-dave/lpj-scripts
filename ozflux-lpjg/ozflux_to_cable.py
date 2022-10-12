#!/usr/bin/env python
# coding: utf-8

""" 
The aim of this notebook is to convert a standard OzFlux L6 netcdf file to 
a CABLE input met netcdf file """


from netCDF4 import Dataset
import netCDF4
import os
import copy
import pandas as pd    # dataframes
import numpy as numpy  # multi-dimensional arrays
import numpy.ma as ma
from shutil import copyfile
import numpy as np
import urllib
import itertools
import datetime

from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import calendar

#rootgrp.close()



""" helpful functions """

def ncdump(nc_fid, verb=True): # from: http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
            - note unicode is just str() in py3
        """
        try:
            print ("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print ('\t\t%s:' % ncattr,repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print ("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print ("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print ('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print ("NetCDF dimension information:")
        for dim in nc_dims:
            print ("\tName:", dim)
            print ("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print ("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print ('\tName:', var)
                print ("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print ("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars



""" way forward - use Peter's code

see: 'convert_v27tov28' function for reading in a nc series and exporting it
-use 'nc_read_series' to open as a dataseries
-see 'nc_open_write' and 'nc_write_series' for writing new netcdf

made this py3 compatible 

"""


class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}
        self.globalattributes["Functions"] = ""
        self.mergeserieslist = []
        self.averageserieslist = []
        self.returncodes = {"value":0,"message":"OK"}

def get_datetimefromnctime(ds,time,time_units):
    """
    Purpose:
     Create a series of datetime objects from the time read from a netCDF file.
    Usage:
     pfp_utils.get_datetimefromnctime(ds,time,time_units)
    Side effects:
     Creates a Python datetime series in the data structure
    Author: PRI
    Date: September 2014
    """
    ts = int(ds.globalattributes["time_step"])
    nRecs = int(ds.globalattributes["nc_nrecs"])
    dt = netCDF4.num2date(time,time_units)
    ds.series[str("DateTime")] = {}
    ds.series["DateTime"]["Data"] = dt
    ds.series["DateTime"]["Flag"] = numpy.zeros(nRecs)
    ds.series["DateTime"]["Attr"] = {}
    ds.series["DateTime"]["Attr"]["long_name"] = "Datetime in local timezone"
    ds.series["DateTime"]["Attr"]["units"] = "None"


def GetSeries(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    """ Returns the data, QC flag and attributes of a series from the data structure."""
    # number of records
    if "nc_nrecs" in ds.globalattributes:
        nRecs = int(ds.globalattributes["nc_nrecs"])
    else:
        nRecs = len(ds.series[ThisOne]["Data"])
    # check the series requested is in the data structure
    if ThisOne in ds.series.keys():
        # series is in the data structure
        if isinstance(ds.series[ThisOne]['Data'],list):
            # return a list if the series is a list
            Series = list(ds.series[ThisOne]['Data'])
        elif isinstance(ds.series[ThisOne]['Data'],numpy.ndarray):
            # return a numpy array if series is an array
            Series = ds.series[ThisOne]['Data'].copy()
        # now get the QC flag
        if 'Flag' in ds.series[ThisOne].keys():
            # return the QC flag if it exists
            Flag = ds.series[ThisOne]['Flag'].copy()
        else:
            # create a QC flag if one does not exist
            Flag = numpy.zeros(nRecs,dtype=numpy.int32)
        # now get the attribute dictionary
        if "Attr" in ds.series[ThisOne].keys():
            Attr = GetAttributeDictionary(ds,ThisOne)
        else:
            Attr = MakeAttributeDictionary()
    else:
        # make an empty series if the requested series does not exist in the data structure
        #logger.warning("GetSeries: "+ThisOne+" not found, making empty series ...")
        Series,Flag,Attr = MakeEmptySeries(ds,ThisOne)
    # tidy up
    if ei==-1: ei = nRecs - 1
    if mode=="truncate":
        # truncate to the requested start and end indices
        si = max(0,si)                  # clip start index at 0
        ei = min(nRecs,ei)              # clip end index to nRecs
        Series = Series[si:ei+1]        # truncate the data
        Flag = Flag[si:ei+1]            # truncate the QC flag
    elif mode=="pad":
        # pad with missing data at the start and/or the end of the series
        if si<0 and ei>nRecs-1:
            # pad at the start
            Series = numpy.append(float(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series)
            Flag = numpy.append(numpy.ones(abs(si),dtype=numpy.int32),Flag)
            # pad at the end
            Series = numpy.append(Series,float(c.missing_value)*numpy.ones((ei-(nRecs-1)),dtype=numpy.float64))
            Flag = numpy.append(Flag,numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si<0 and ei<=nRecs-1:
            # pad at the start, truncate the end
            Series = numpy.append(float(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series[:ei+1])
            Flag = numpy.append(numpy.ones(abs(si),dtype=numpy.int32),Flag[:ei+1])
        elif si>=0 and ei>nRecs-1:
            # truncate at the start, pad at the end
            Series = numpy.append(Series[si:],float(c.missing_value)*numpy.ones((ei-(nRecs-1)),numpy.float64))
            Flag = numpy.append(Flag[si:],numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si>=0 and ei<=nRecs-1:
            # truncate at the start and end
            Series = Series[si:ei+1]
            Flag = Flag[si:ei+1]
        else:
            msg = 'GetSeries: unrecognised combination of si ('+str(si)+') and ei ('+str(ei)+')'
            raise ValueError(msg)
    elif mode=="mirror":
        # reflect data about end boundaries if si or ei are out of bounds
        if si<0 and ei>nRecs-1:
            # mirror at the start
            Series = numpy.append(numpy.fliplr([Series[1:abs(si)+1]])[0],Series)
            Flag = numpy.append(numpy.fliplr([Flag[1:abs(si)+1]])[0],Flag)
            # mirror at the end
            sim = 2*nRecs-1-ei
            eim = nRecs-1
            Series = numpy.append(Series,numpy.fliplr([Series[sim:eim]])[0])
            Flag = numpy.append(Flag,numpy.fliplr([Flag[sim:eim]])[0])
        elif si<0 and ei<=nRecs-1:
            # mirror at start, truncate at end
            Series = numpy.append(numpy.fliplr([Series[1:abs(si)+1]])[0],Series[:ei+1])
            Flag = numpy.append(numpy.fliplr([Flag[1:abs(si)+1]])[0],Flag[:ei+1])
        elif si>=0 and ei>nRecs-1:
            # truncate at start, mirror at end
            sim = 2*nRecs-1-ei
            eim = nRecs
            Series = numpy.append(Series[si:],numpy.fliplr([Series[sim:eim]])[0])
            Flag = numpy.append(Flag[si:],numpy.fliplr([Flag[sim:eim]])[0])
        elif si>=0 and ei<=nRecs-1:
            # truncate at the start and end
            Series = Series[si:ei+1]
            Flag = Flag[si:ei+1]
        else:
            msg = 'GetSeries: unrecognised combination of si ('+str(si)+') and ei ('+str(ei)+')'
            raise ValueError(msg)
    else:
        raise ValueError("GetSeries: unrecognised mode option "+str(mode))
    return Series,Flag,Attr

    
def nc_read_series(ncFullName,checktimestep=True,fixtimestepmethod="round"):
    """
    Purpose:
     Reads a netCDF file and returns the meta-data and data in a DataStructure.
     The returned data structure is an instance of pfp_io.DataStructure().
     The data structure consists of:
      1) ds.globalattributes
         A dictionary containing the global attributes of the netCDF file.
      2) ds.series
         A dictionary containing the variable data, meta-data and QC flag
         Each variable dictionary in ds.series contains;
         a) ds.series[variable]["Data"]
            A 1D numpy float64 array containing the variable data, missing
            data value is -9999.
         b) ds.series[variable]["Flag"]
            A 1D numpy int32 array containing the QC flag data for this variable.
         c) ds.series[variable]["Attr"]
            A dictionary containing the variable attributes.
    Usage:
     nc_name = pfp_io.get_filename_dialog(path="../Sites/Whroo/Data/Processed/")
     ds = pfp_io.nc_read_series(nc_name)
     where nc_name is the full name of the netCDF file to be read
           ds is the returned data structure
    Side effects:
     This routine checks the time step of the data read from the netCDF file
     against the value of the global attribute "time_step", see pfp_utils.CheckTimeStep.
     If a problem is found with the time step (duplicate records, non-integral
     time steps or gaps) then pfp_utils.FixTimeStep is called to repair the time step.
     Fixing non-integral timne steps requires some user input.  The options are to
     quit ([Q]), interpolate ([I], not implemented yet) or round ([R]).  Quitting
     causes the script to exit and return to the command prompt.  Interpolation
     is not implemented yet but will interpolate the data from the original time
     step to a regular time step.  Rounding will round any non-itegral time steps
     to the nearest time step.
    Author: PRI
    Date: Back in the day
    """
    #logger.info(" Reading netCDF file "+ntpath.split(ncFullName)[1])
    netCDF4.default_encoding = 'latin-1'
    ds = DataStructure()
    # check to see if the requested file exists, return empty ds if it doesn't
#     if ncFullName[0:4]!="http":
#         if not pfp_utils.file_exists(ncFullName,mode="quiet"):
#             logger.error(' netCDF file '+ncFullName+' not found')
#             raise Exception("nc_read_series: file not found")
    # file probably exists, so let's read it
    ncFile = netCDF4.Dataset(ncFullName,'r')
    # disable automatic masking of data when valid_range specified
    ncFile.set_auto_mask(False)
    # now deal with the global attributes
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
    # get a list of the variables in the netCDF file (not their QC flags)
    varlist = [x for x in ncFile.variables.keys() if "_QCFlag" not in x]
    for ThisOne in varlist:
        # skip variables that do not have time as a dimension
        dimlist = [x.lower() for x in ncFile.variables[ThisOne].dimensions]
        if "time" not in dimlist: continue
        # create the series in the data structure
        ds.series[str(ThisOne)] = {}
        # get the data and the QC flag
        data,flag,attr = nc_read_var(ncFile,ThisOne)
        ds.series[ThisOne]["Data"] = data
        ds.series[ThisOne]["Flag"] = flag
        ds.series[ThisOne]["Attr"] = attr
    ncFile.close()
    # make sure all values of -9999 have non-zero QC flag
    # NOTE: the following was a quick and dirty fix for something a long time ago
    #       and needs to be retired
    #pfp_utils.CheckQCFlags(ds)
    # get a series of Python datetime objects
    if "time" in ds.series.keys():
        time,f,a = GetSeries(ds,"time")
        get_datetimefromnctime(ds,time,a["units"])
    else:
        get_datetimefromymdhms(ds)
    # round the Python datetime to the nearest second
#     pfp_utils.round_datetime(ds,mode="nearest_second")
    # check the time step and fix it required
#     if checktimestep:
#         if pfp_utils.CheckTimeStep(ds):
#             pfp_utils.FixTimeStep(ds,fixtimestepmethod=fixtimestepmethod)
#             # update the Excel datetime from the Python datetime
#             pfp_utils.get_xldatefromdatetime(ds)
#             # update the Year, Month, Day etc from the Python datetime
#             pfp_utils.get_ymdhmsfromdatetime(ds)
    # tell the user when the data starts and ends
    ldt = ds.series["DateTime"]["Data"]
    msg = " Got data from "+ldt[0].strftime("%Y-%m-%d %H:%M:%S")+" to "+ldt[-1].strftime("%Y-%m-%d %H:%M:%S")
#    logger.info(msg)
    return ds

def nc_read_var(ncFile,ThisOne):
    """ Reads a variable from a netCDF file and returns the data, the QC flag and the variable
        attribute dictionary.
    """
    # check the number of dimensions
    nDims = len(ncFile.variables[ThisOne].shape)
    if nDims not in [1,3]:
        msg = "nc_read_var: unrecognised number of dimensions ("+str(nDims)
        msg = msg+") for netCDF variable "+ ThisOne
        raise Exception(msg)
    if nDims==1:
        # single dimension
        data = ncFile.variables[ThisOne][:]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data,dummy = MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    elif nDims==3:
        # 3 dimensions
        data = ncFile.variables[ThisOne][:,0,0]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        # may not be needed after adding ncFile.set_auto_mask(False) in nc_read_series().
        if numpy.ma.isMA(data): data,dummy = MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:,0,0]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    # force float32 to float64
    if data.dtype=="float32": data = data.astype(numpy.float64)
    # check for Year, Month etc as int64, force to int32 if required
    if ThisOne in ["Year","Month","Day","Hour","Minute","Second"]:
        if data.dtype=="int64": data = data.astype(numpy.int32)
    # get the variable attributes
    vattrlist = ncFile.variables[ThisOne].ncattrs()
    attr = {}
    if len(vattrlist)!=0:
        for vattr in vattrlist:
            attr[vattr] = getattr(ncFile.variables[ThisOne],vattr)
    return data,flag,attr

def GetSeriesasMA(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    """
    Purpose:
     Returns a data series and the QC flag series from the data structure.
    Usage:
     data,flag,attr = pfp_utils.GetSeriesasMA(ds,label,si=0,ei=-1)
    where the arguments are;
      ds    - the data structure (dict)
      label - label of the data series in ds (string)
      si    - start index (integer), default 0
      ei    - end index (integer), default -1
    and the returned values are;
      data - values for the requested series in ds
             (numpy masked array, float64)
      flag - QC flag for the requested series in ds
             (numpy masked array, int32)
      attr - attribute dictionary for series
    Example:
     The code snippet below will return the incoming shortwave data values
     (Fsd) and the associated QC flag (f) as numpy masked arrays;
      ds = pfp_io.nc_read_series("HowardSprings_2011_L3.nc")
      Fsd,f,a = pfp_utils.GetSeriesasMA(ds,"Fsd")
    Author: PRI
    """
    Series,Flag,Attr = GetSeries(ds,ThisOne,si=si,ei=ei,mode=mode)
    Series,WasND = SeriestoMA(Series)
    return Series,Flag,Attr

def get_datetimefromymdhms(ds):
    ''' Creates a series of Python datetime objects from the year, month,
    day, hour, minute and second series stored in the netCDF file.'''
    SeriesList = ds.series.keys()
    if ('Year' not in SeriesList or 'Month' not in SeriesList or 'Day' not in SeriesList or
        'Hour' not in SeriesList or 'Minute' not in SeriesList or 'Second' not in SeriesList):
#        logger.info(' get_datetimefromymdhms: unable to find all datetime fields required')
        return
#    logger.info(' Getting the date and time series')
    year = ds.series["Year"]["Data"]
    month = ds.series["Month"]["Data"]
    day = ds.series["Day"]["Data"]
    hour = ds.series["Hour"]["Data"]
    minute = ds.series["Minute"]["Data"]
    second = ds.series["Second"]["Data"]
    dt = [datetime.datetime(yr,mn,dy,hr,mi,se) for yr,mn,dy,hr,mi,se in zip(year,month,day,hour,minute,second)]
    ds.series["DateTime"] = {}
    ds.series["DateTime"]["Data"] = numpy.array(dt)
    ds.series["DateTime"]["Flag"] = numpy.zeros(len(dt))
    ds.series["DateTime"]["Attr"] = {}
    ds.series["DateTime"]["Attr"]["long_name"] = "Datetime in local timezone"
    ds.series["DateTime"]["Attr"]["units"] = "None"
    return

def MAtoSeries(Series):
    """
    Convert a masked array to a numpy ndarray with masked elements set to c.missing_value.
    Useage:
     Series, WasMA = MAtoSeries(Series)
     where:
      Series (input)    is the data series to be converted.
      WasMA  (returned) is a logical, True if the input series was a masked array.
      Series (output)   is the input series convered to an ndarray with c.missing_value values
                        for missing data.
    """
    WasMA = False
    if numpy.ma.isMA(Series):
        WasMA = True
        Series = numpy.ma.filled(Series,float(c.missing_value))
    return Series, WasMA

def GetAttributeDictionary(ds,ThisOne):
    attr = {}
    # if series ThisOne is in the data structure
    if ThisOne in ds.series.keys():
        attr = ds.series[ThisOne]['Attr']
    else:
        attr = MakeAttributeDictionary()
    return copy.deepcopy(attr)

def SeriestoMA(Series):
    """
    Convert a numpy ndarray to a masked array.
    Useage:
     Series, WasND = SeriestoMA(Series)
     where:
      Series (input)    is the data series to be converted.
      WasND  (returned) is a logical, True if the input series was an ndarray
      Series (output)   is the input series convered to a masked array.
    """
    WasND = False
    if Series.dtype == "float64":
        if not numpy.ma.isMA(Series):
            WasND = True
            Series = numpy.ma.masked_where(abs(Series-numpy.float64(c.missing_value)) < c.eps, Series)
    return Series, WasND



def create_file_from_source(src_file, trg_file):
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        #trg.createVariable(name, var.dtype, var.dimensions)
        
        # Copy the variable attributes
        #trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        #trg.variables[name][:] = src.variables[name][:]
        if var in in_var_list:
            print (var)
        
    # Save the file
    trg.close()

#create_file_from_source(os.path.join(i_path, i_file), os.path.join(o_path, o_file))
#
#
#
#
#with netCDF4.Dataset("in.nc") as src, netCDF4.Dataset("out.nc", "w") as dst:
#    # copy global attributes all at once via dictionary
#    dst.setncatts(src.__dict__)
#    # copy dimensions
#    for name, dimension in src.dimensions.items():
#        dst.createDimension(
#            name, (len(dimension) if not dimension.isunlimited() else None))
#    # copy all file data except for the excluded
#    for name, variable in src.variables.items():
#        if name in in_var_list:
#            x = dst.createVariable(name, variable.datatype, variable.dimensions)
#            dst[name][:] = src[name][:]
#            # copy variable attributes all at once via dictionary
#            dst[name].setncatts(src[name].__dict__)
#
#
#
#import netCDF4 as nc
#
#src = Dataset(os.path.join(i_path, i_file), "r")
#trg = Dataset(os.path.join(o_path, o_file), mode='w')
#
#
## copy global attributes all at once via dictionary
#trg.setncatts(src.__dict__)
## copy dimensions
#for name, dimension in src.dimensions.items():
#    trg.createDimension(
#        name, (len(dimension) if not dimension.isunlimited() else None))
## copy all file data except for the excluded
#for name, variable in src.variables.items():
#    if name in in_var_list:
#        x = trg.createVariable(name, variable.datatype, variable.dimensions)
#        trg[name][:] = src[name][:]
#        # copy variable attributes all at once via dictionary
#        trg[name].setncatts(src[name].__dict__)






###########################
# # Working version #######
###########################
## TODO:
    
generate_site_list      = False    # print site list after the loop?
full_years_only         = True     # take full years only?
add_end_start_timesteps = True     # replicate first/last timesteps if they are missing?
                                   # timeseries is expected to start at 00:00 and end at 23:00
multiple_tiles          = True     # 1 (False) or more tiles (True)?
OzFlux_default          = True     # use default OzFlux processing (True) or site-specific processing 
                                   # as provided by the site PIs (False)? ATTENTION: only set up for CumberlandPlain at the moment!!!! check LAI calculation if applied to other sites
Elise_2019              = False    # special case: EC data from Elise (Dec. 2019) with specific requirements
Will_2019               = False    # process Tumbarumba 2019 and add to existing file

## constants
day2sec     = 86400.
Kelvin      = 273.15
k           = 0.5     # extinction coefficient in the LAI calculations (Trudinger et al. 2016, Biogeosciences)
fcorr_patch = 0.65    # correction factor related to fraction of woody and grass (nadir vs. projected)


## paths
basepath='C:/Users/kna016/Documents/CABLE/'

## stuff related to cover fraction file
#fcover   = pd.read_csv(basepath + 'cover_fract/ozflux28.gimms3g_fcovfperfrec.csv',index_col='PointName')  # old one with 28 sites
fcover   = pd.read_csv(basepath + 'cover_fract/ozflux30.gimms3g_fcovfperfrec.csv',index_col='PointName') # new one with 30 sites
fcov_index = [k for k, i in enumerate(fcover.index) if '9999' in i and 'fcov' in i]
fper_index = [k for k, i in enumerate(fcover.index) if '9999' in i and 'fper' in i]
frec_index = [k for k, i in enumerate(fcover.index) if '9999' in i and 'frec' in i]
fcover_sites = list(fcover.columns)


## initialize lists
canopy_height = []
latitude      = []
longitude     = []
tower_height  = []
startyear     = []
endyear       = []
vegetation    = [] 


## OzFlux sites
if Elise_2019:
  all_sites = ('CumberlandPlain',)
elif Will_2019:   
  all_sites = ('Tumbarumba',)
else:
  all_sites = 'AdelaideRiver', 'AliceSpringsMulga', 'Boyagin', 'Calperum', 'CapeTribulation', 'CowBay', 'CumberlandPlain', 'DalyPasture', 'DalyUncleared', 'Dargo', \
  'DryRiver', 'Emerald', 'FoggDam', 'Gingin', 'GreatWesternWoodlands', 'HowardSprings', 'Litchfield', 'Loxton', 'Nimmo', 'Otway', 'RedDirtMelonFarm', \
  'Ridgefield', 'RiggsCreek', 'RobsonCreek', 'Samford', 'SturtPlains', 'TiTreeEast', 'Tumbarumba', 'WallabyCreek', 'Warra', 'Whroo', 'WombatStateForest', \
  'Yanco'
#all_sites = ('Boyagin','WallabyCreek')
#all_sites = ('Boyagin',) 

elevation = list(itertools.repeat(-9999,len(all_sites)))
#elevation[-6] = 1200


## only used if LAI is not provided
LAI = list(itertools.repeat(-9999,len(all_sites)))
#LAI[-6] = 1.8


### 0) load site list and check if all sites are included
sitelist = pd.read_csv('C:/Users/kna016/Documents/CABLE/OzFLUX_sitelist_v6.txt',sep="\t",index_col=0)

## short excercise: extract all start and end dates from the sites
#start_date = []
#end_date = []
#data_sites = []
#
#for s,site in enumerate(all_sites):
#  print('starting site nr. ' + str(s+1) + ' (' + site + ')')
#  
#  outfile='E:/CSIRO/CABLE/OzFluxL6/' + site + '.nc'
#  url='http://dap.ozflux.org.au/thredds/fileServer/ozflux/sites/' + site + '/L6/default/' + site + '_L6.nc'  # note: 'fileServer' instead of 'catalog' (was not working with 'dodsC')      
#  try:
#    urllib.request.urlretrieve(url,outfile)   
#  
#    ds = nc_read_series(outfile,checktimestep=False)
#    start_date.append(ds.globalattributes['start_date'])
#    end_date.append(ds.globalattributes['end_date'])
#    data_sites.append(site)
#  
#  except urllib.error.HTTPError:
#    print("no data available online for site '" + site + "'! Continuing with next site.")
#
#
#start_end_date_df = pd.DataFrame(index=data_sites)
#
#start_end_date_df['start_date'] = start_date
#start_end_date_df['end_date'] = end_date
#
#### export dataframe to file
#start_end_date_df.to_csv('E:/CSIRO/CABLE/OzFLUX_startend_dates_June_2020.csv')
#start_end_date_df.to_csv('E:/CSIRO/CABLE/OzFLUX_startend_dates_June_2020.txt')




data_sites = []

for s,site in enumerate(all_sites):
   
  print('starting site nr. ' + str(s+1) + ' (' + site + ')')
  
  
  ### 1) download nc files from the Ozflux webpage
  now = datetime.datetime.now()
  outdir = 'C:/Users/kna016/Documents/CABLE/OzFluxL6/' + str(now.year) + '_' + str(now.month) + '_' + str(now.day) + '/'
  #outfile='D:/CSIRO/CABLE/OzFluxL6/' + site + '.nc'
  #outfile='C:/Users/kna016/Documents/CABLE/forcing_files/Cumberland_Elise/L6_merged_filled.nc'
  #outfile='C:/Users/kna016/Documents/CABLE/OzFluxL6/2019_6_17/CumberlandPlain_L6.nc'
  try:
    os.mkdir(outdir)
  except FileExistsError:
    pass
  
  
  metfile_present = True
  iveg
  if not Elise_2019:
    if Will_2019:
      outfile = 'C:/Users/kna016/Documents/CABLE/OzFluxL6/Tumbarumba_2019_L6_all_Gill_EP_wFilter_GUI_v2_noEvG.nc'
    else:
      l = 6
      while l > 3:
        ls = str(l)
        try:
          outfile=outdir + site + '_L' + ls + '.nc'
          if OzFlux_default:
            url='http://dap.ozflux.org.au/thredds/fileServer/ozflux/sites/' + site + '/L' + ls + '/default/' + site + '_L' + ls + '.nc'  # note: 'fileServer' instead of 'catalog' (was not working with 'dodsC') 
          else:
            url='http://dap.ozflux.org.au/thredds/fileServer/ozflux/sites/' + site + '/L' + ls + '/site_pi/' + site + '_L' + ls + '.nc'
          urllib.request.urlretrieve(url,outfile)  # this is python 3 code. Use urllib.urlretrieve for python 2
          break
      
        except urllib.error.HTTPError:
          if l > 4:
            print("no L" + ls + " data available online for site '" + site + "'! Trying L" + str(l-1) + " next.")
          else:
            print("no data available online for site '" + site + "'! Continue with next site.")
            metfile_present=False
          l=l-1
     
     
  if metfile_present:
    ### 2) load nc file and write into a structured 'OzFlux' object and retrieve meteorological variables 
    ds = nc_read_series(outfile,checktimestep=False)
    
    # delete nc_nrecs global attribute if not consistent with the actual length
    if ds.globalattributes["nc_nrecs"] != len(ds.series['time']['Data']):
      del ds.globalattributes["nc_nrecs"]
    
    
    l6_Ws,l6_Ws_flag,Ws_attr = GetSeries(ds,"Ws")
    l6_time,l6_time_flag,time_attr = GetSeries(ds,"time")
    l6_Fsd,l6_Fsd_flag,Fsd_attr = GetSeries(ds,"Fsd")
    l6_Ta,l6_Ta_flag,Ta_attr = GetSeries(ds,"Ta")e
    l6_Precip,l6_Precip_flag,Precip_attr = GetSeries(ds,"Precip")
    if OzFlux_default and not Will_2019:
      l6_q,l6_q_flag,q_attr = GetSeries(ds,"SH")       # 'q' in previous version
    else:
      l6_q,l6_q_flag,q_attr = GetSeries(ds,"q")
    l6_ps,l6_ps_flag,ps_attr = GetSeries(ds,"ps")
    l6_Fld,l6_Fld_flag,Fld_attr = GetSeries(ds,"Fld")
    
    
    ## throw out redundant attributes (these may cause problems on raijin)
#    del_entries = ('instrument','ancillary_variables','serial_number','height','coverage_L2',
#                   'standard_name','coverage_L6','coverage_L5','coverage_L4','coverage_L3','valid_range',
#                   'rangecheck_upper','rangecheck_lower')
#    for key in del_entries:
#      del Ws_attr[key]
#      del Fsd_attr[key]
#      del Ta_attr[key]
#      del Precip_attr[key]
#      if key=='coverage_L2' or key=='rangecheck_upper' or key=='rangecheck_lower':  
#        if OzFlux_default:
#          del q_attr[key]
#      del ps_attr[key]
#      del Fld_attr[key]
#    
#    del time_attr['calendar'], time_attr['long_name'], time_attr['standard_name']
    
    if Elise_2019:
      ## set missing value and delete attributes
      Ws_attr['missing_value']=-9999.0
      Fsd_attr['missing_value']=-9999.0
      Ta_attr['missing_value']=-9999.0
      Precip_attr['missing_value']=-9999.0
      q_attr['missing_value']=-9999.0
      ps_attr['missing_value']=-9999.0
      Fld_attr['missing_value']=-9999.0
      
      del_entries = ('_FillValue', 'nrecs')
      for key in del_entries:
          del Ws_attr[key]
          del Fsd_attr[key]
          del Ta_attr[key]
          del Precip_attr[key]
          del ps_attr[key]
          del Fld_attr[key]
          
      del q_attr['_FillValue']    
      del time_attr['_FillValue']    
      
      ## set missing value to mean of the other years
      # Note: the following is poor coding, ideally one would integrate
      # over time index using pandas...
      nr_missing = sum(l6_Fsd < -9998)
      
      vars_correct = ('Ws','Fsd','Ta','Precip','ps','Fld')
      for varname in vars_correct:
        var = globals()['l6_' + varname]
      
        year1 = var[0:17520]     # 2014
        year2 = var[17520:35040]
        year3 = var[35040:52608]
        year4 = var[52608:70128]
        year5 = var[70128:87648]
      
        vr =  range((len(year1) - nr_missing),len(year1))
        vr1 = range((len(year3) - nr_missing),len(year3))
        meanmiss = (year1[vr] + year2[vr] + year3[vr1] + year4[vr] + year5[vr])/5
      
        var[len(var)-nr_missing:len(var)] = meanmiss
      
        globals()['l6_' + varname] = var
        
      l6_Ws_flag[numpy.where(l6_Ws_flag < -9998)] = -1.0
      l6_Fsd_flag[numpy.where(l6_Fsd_flag < -9998)] = -1.0
      l6_Ta_flag[numpy.where(l6_Ta_flag < -9998)] = -1.0
      l6_Precip_flag[numpy.where(l6_Precip_flag < -9998)] = -1.0
      l6_q_flag[numpy.where(l6_q_flag < -9998)] = -1.0
      l6_ps_flag[numpy.where(l6_ps_flag < -9998)] = -1.0
      l6_Fld_flag[numpy.where(l6_Fld_flag < -9998)] = -1.0
      
    else:
      Ws_attr['missing_value']=float(Ws_attr['missing_value'])
      Fsd_attr['missing_value']=float(Fsd_attr['missing_value'])
      Ta_attr['missing_value']=float(Ta_attr['missing_value'])
      Precip_attr['missing_value']=float(Precip_attr['missing_value'])
      q_attr['missing_value']=float(q_attr['missing_value'])
      ps_attr['missing_value']=float(ps_attr['missing_value'])
      Fld_attr['missing_value']=float(Fld_attr['missing_value'])
    
    ## TODO: change long_name attribute to a more reasonable value!  
      
    
    
    
    ## checking for NANs
#    if OzFlux_default:
#      new_list = ['Ws','Fsd','Ta','Precip','SH','ps','Fld']
#      for var in new_list:
#        data_vals, flags, time_att = GetSeries(ds,var)
#        if len(data_vals[np.isnan(data_vals)]) > 0:
#          print (var, 'found nans!!')
#        print (var, max(data_vals),min(data_vals))
      
    
    ### 3) Extract start- and end year (at least in Tumbarumba, first time step is 1 am, not midnight)    
    start_year = ds.globalattributes['start_date'][0:4]
    end_year = ds.globalattributes['end_date'][0:4]
    start_day =  ds.globalattributes['start_date'][5:10]
    end_day = ds.globalattributes['end_date'][5:10]
    start_hour = ds.globalattributes['start_date'][-8:-3]
    end_hour = ds.globalattributes['end_date'][-8:-3]
  
    ## extract date from original nc (might be a better way)
    ncd = netCDF4.Dataset(outfile)
    time_orig = ncd.variables['time']
    date = netCDF4.num2date(time_orig[:],time_orig.units,time_orig.calendar)
    timestep = (date[1] - date[0]).seconds  # halfhourly or hourly?
  
  
    # Convert time from days to seconds and get length of ts
    sec_l6_time = l6_time*day2sec-l6_time[0]*day2sec
    sec_l6_time = np.array(sec_l6_time,dtype=int)  # starts at 0
      
    startdate_orig = date[0]
    startdate = startdate_orig  # needed for time axis in nc file (see end of loop)
    enddate_orig = date[len(date) -1]
    
  
  
    ### 4) if necessary add timestep to work with full years
    # what does CABLE expect here?
    if add_end_start_timesteps:
      start_added_one = False
      start_deleted_one = False
      end_added_one = False
      end_added_two = False
      valid_days_start = list(range(0,len(date)))
      valid_days_end = list(range(0,len(date)))
      ctr=0
      if (start_hour != '00:30' and timestep == 1800) or (start_hour != '01:00' and timestep == 3600.):  # normal case: beginning of the day
        #if (start_day == 1):  # otherwise don't bother for now
        if start_hour == '00:00':
          sec_l6_time = np.append(sec_l6_time[0]-timestep,sec_l6_time)
          #date = np.delete(date,0)
          del valid_days_start[0]
          del valid_days_end[0]
          ctr=1
          print('one timestep deleted at the beginning of the timeseries')
          start_deleted_one = True
          startdate = startdate_orig + datetime.timedelta(seconds=timestep)
          # elif ### case beginning added missing, but doesn't exist right now!
        else:
          print("unusual start hour! First day of time series is omitted.")
          ## does not work now for full days, but ok for full years:
          valid_days_start = list(np.where(date > datetime.datetime(startdate_orig.year,startdate_orig.month,startdate_orig.day) + datetime.timedelta(days=1,seconds=1))[0])    ## check here: next day starts at 00:00, but should be 00:30
        
      if end_hour != '00:00':
        if end_hour == '23:30' or end_hour == '23:00' or end_hour == '22:00':  # 22:30 could be added, but doesn't exist right now
          if (full_years_only==False and end_day != '12-31') or end_day == '12-31':
            if (end_hour == '22:00' and timestep == 3600.) or (end_hour == '23:00' and timestep == 1800.):
              sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
              sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
              date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep))
              date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep*2))
              valid_days_start.append(len(valid_days_end)+ctr)
              valid_days_end.append(len(valid_days_end)+ctr)
              valid_days_start.append(len(valid_days_end)+ctr)
              valid_days_end.append(len(valid_days_end)+ctr)
              print('two timesteps added at the end of the timeseries')
              end_added_one = True
              end_added_two = True
            elif (end_hour == '23:00' and timestep == 3600.) or (end_hour == '23:30' and timestep == 1800.):
              sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
              date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep))
              valid_days_start.append(len(valid_days_end)+ctr)
              valid_days_end.append(len(valid_days_end)+ctr)
              print('one timestep added at the end of the timeseries')
              end_added_one = True
        else:       # delete the whole day
          print('unusual end hour! Last day of the time series is omitted.')
          #valid_days_end = list(np.where(date <= datetime.datetime(enddate_orig.year,enddate_orig.month,enddate_orig.day))[0])
          valid_days_end = list(np.where(date < datetime.datetime(enddate_orig.year,enddate_orig.month,enddate_orig.day))[0])
          
## old version, assuming that the time stamp marks the beginning of the time step:        
#      if start_hour != '00:00':
#        if (full_years_only==False and start_day != '01-01') or start_day == '01-01':   # don't add the half hour because year will be deleted later
#          if start_hour == '01:00' or start_hour == '00:30':
#            sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
#            date = numpy.append(startdate_orig - datetime.timedelta(seconds=timestep),date)
#            valid_days_start.append(len(valid_days_start))
#            valid_days_end.append(len(valid_days_end))
#            print('one timestep added at the beginning of the timeseries')
#            start_added_one = True
#            startdate = startdate_orig - datetime.timedelta(seconds=timestep)
#          else:
#            print("unusual start hour! First day of time series is omitted.")
#            valid_days_start = list(np.where(date >= startdate_orig + datetime.timedelta(days=1))[0])
#      if end_hour != '23:30':
#        if (full_years_only==False and end_day != '12-31') or end_day == '12-31':
#          if end_hour == '22:00' and timestep == 3600.:
#            sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
#            date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep))
#            valid_days_start.append(len(valid_days_start))
#            valid_days_end.append(len(valid_days_end))
#            print('one timestep added at the end of the timeseries')
#            end_added_one = True
#          elif end_hour == '23:00' and timestep == 1800.:
#            sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
#            date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep))
#            valid_days_start.append(len(valid_days_start))
#            valid_days_end.append(len(valid_days_end))
#            print('one timestep added at the end of the timeseries')
#            end_added_one = True
#          elif end_hour == '22:30' and timestep == 1800.:
#            sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
#            sec_l6_time = numpy.append(sec_l6_time,sec_l6_time[-1] + timestep)
#            date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep))
#            date = numpy.append(date,enddate_orig + datetime.timedelta(seconds=timestep*2))
#            valid_days_start.append(len(valid_days_start))
#            valid_days_start.append(len(valid_days_start))
#            valid_days_end.append(len(valid_days_end))
#            valid_days_end.append(len(valid_days_end))
#            print('two timesteps added at the end of the timeseries')
#            end_added_one = True
#            end_added_two = True
#          else:       # delete the whole day
#            print('unusual end hour! Last day of the time series is omitted.')
#            valid_days_end = list(np.where(date < datetime.datetime(enddate_orig.year,enddate_orig.month,enddate_orig.day))[0])
          
    
   
    valid_days = list(set(valid_days_start) & set(valid_days_end))
     
    # update start day and month       
    start_day   = date[valid_days[0]].day
    start_month = date[valid_days[0]].month
    end_day     = date[valid_days[-1]].day
    end_month   = date[valid_days[-1]].month
    
    
    ## checks:
    # valid_days_start[0]   
    # valid_days_start[-1] 
    # valid_days_end[0]   
    # valid_days_end[-1] 
    # valid_days[0]   
    # valid_days[-1] 
    
    
      
    ### 5) take full years only
    valid_years = list(range(0,len(date)))
    if full_years_only:
      valid_years_start = list(range(0,len(date)))
      valid_years_end = list(range(0,len(date)))
      if (str(start_day) + '-' + str(start_month) != '1-1'):  # delete first year
      #if start_day != 1 and start_month != 1:  # delete first year
        print('first year is incomplete and is omitted.')
        valid_years_start = list(np.where(date >= datetime.datetime(int(start_year) + 1,1,1,0,1))[0])
        #start_added_one = False
      if (str(end_day) + '-' + str(end_month) != '31-12' and str(end_day) + '-' + str(end_month) != '1-1'): # delete last year
      #if (end_day != 31 and end_month != 12) and (end_day != 1 and end_month != 1):    # delete last year
        print('last year is incomplete and is omitted.')
        #valid_years_end = list(np.where(date < datetime.datetime(int(end_year) - 1,12,31,23,59))[0])
        valid_years_end = list(np.where(date < datetime.datetime(int(end_year),1,1,0,0,1))[0])
        #end_added_one = False
        #end_added_two = False
     
      valid_years = list(set(valid_years_start) & set(valid_years_end))  # intersection of ind1 and ind2
      valid_years.sort()  # no idea why it wasn't sorted for one site
  
    valid = list(set(valid_days) & set(valid_years))
    valid.sort()
        
    ## some checks:
    # valid_years_start[0]   
    # valid_years_start[-1] 
    # valid_years_end[0]   
    # valid_years_end[-1] 
    # valid_years[0]   
    # valid_years[-1] 

    # date[valid[0]]
    # date[valid[-1]]
    # len(valid)
    
    if len(valid) < 1:
      print('no full years available!')
    else:  
      data_sites.append(site)
      
      # note: sec_l6_time has to start at 0!
      sec_l6_time = sec_l6_time[valid]
      ts_len = len(sec_l6_time)
  
      
      
      
  
      ### 6) extract other relevant attributes
      if generate_site_list:
        try:
          canopy_height_string = ds.globalattributes['canopy_height']
            
          if canopy_height_string[-2:] == 'cm':
            canopy_height.append(float(canopy_height_string[0:2])/100)
          else:
            try:
              canopy_height.append(float(canopy_height_string[0:2]))
            except ValueError:
              if len(canopy_height_string) > 0:
                  try:
                    canopy_height.append(float(canopy_height_string[0:1]))
                  except ValueError:
                    print("canopy height is not available")
                    canopy_height.append(None)
              else:
                print("canopy height is not available")
                canopy_height.append(None)
          
        except KeyError:
          print("canopy height is not available")
          canopy_height.append(None)
            
          
                
        try:
          tower_height.append(float(ds.globalattributes['tower_height'][0:2]))
        except (ValueError, KeyError):
          print("tower height is not available")
          tower_height.append(None)
        
        if site == 'Gingin':
          latitude.append(-31.375)
          longitude.append(115.65)
        else:    
          latitude.append(float(ds.globalattributes['latitude']))
          longitude.append(float(ds.globalattributes['longitude']))
      
      
        try:
          vegetation.append(ds.globalattributes['vegetation'])
        except KeyError:
          try:
            vegetation.append(ds.globalattributes['Vegetation'])
          except KeyError:
            print('Vegetation type is not available')
            vegetation.append(None)
      
      
      
      ### 7) create empty nc file from scratch
      if full_years_only == True and add_end_start_timesteps == True:
        path_met = 'C:/Users/kna016/Documents/CABLE/forcing_files/fullyears_fulldays/'
      elif full_years_only == True and add_end_start_timesteps == False:
        path_met = 'C:/Users/kna016/Documents/CABLE/forcing_files/fullyears'
      elif full_years_only == False and add_end_start_timesteps == True:
        path_met = 'C:/Users/kna016/Documents/CABLE/forcing_files/fulldays'
      elif full_years_only == False and add_end_start_timesteps == False:
        path_met = 'C:/Users/kna016/Documents/CABLE/forcing_files/no_change'
        
    
        
      # outfile_met = site + '_' + start_year + '_' + end_year + '.nc'
      str_startyear = str(date[valid][0].year)
      if full_years_only:
        str_endyear = str(date[valid][-1].year - 1)
      else:
        str_endyear = str(date[valid][-1].year)
      
      if multiple_tiles:
        if OzFlux_default:
          outfile_met = site + '_' + str_startyear + '_' + str_endyear + '_2tiles.nc' 
        else:
          outfile_met = site + '_' + 'site_pi' + '_' + str_startyear + '_' + str_endyear + '_2tiles.nc' 
      else:
        if OzFlux_default:
          outfile_met = site + '_' + str_startyear + '_' + str_endyear + '.nc' 
        else:
          outfile_met = site + '_' + 'site_pi' + '_' + str_startyear + '_' + str_endyear + '.nc' 
          
      
      startyear.append(int(str_startyear))
      endyear.append(int(str_endyear))

      
      #nc_file = Dataset(os.path.join(path_met, outfile_met), "w", format="NETCDF3_CLASSIC")
      nc_file = Dataset(os.path.join(path_met, outfile_met), "w", format="NETCDF4")
        
      # 7.1) create dimensions
      nc_file.createDimension("x", 1)
      nc_file.createDimension("y", 1)
      nc_file.createDimension("time", ts_len) # this will need to match the length of the L6 file
      #nc_file.createDimension("time", None)  # this should create an 'UNLIMITED' time axis
      nc_file.createDimension("z", 1)
      if multiple_tiles:
        nc_file.createDimension("patch",2)
      
      
      # 7.2) create variables
      x_val = nc_file.createVariable("x","float64",("x",))
      y_val = nc_file.createVariable("y","float64",("y",))
      lat   = nc_file.createVariable("latitude","float32",("y","x",))
      lon   = nc_file.createVariable("longitude","float32",("y","x",))
      times = nc_file.createVariable("time","int32",("time",))
      fsd   = nc_file.createVariable("SWdown","float32",("time","y","x",))
      z_val = nc_file.createVariable("z","float64",("z",))
      tair  = nc_file.createVariable("Tair","float32",("time","z","y","x",))
      precip= nc_file.createVariable("Rainf","float32",("time","y","x",))
      qair  = nc_file.createVariable("Qair","float32",("time","z","y","x",))
      ws    = nc_file.createVariable("Wind","float32",("time","z","y","x",))
      ps    = nc_file.createVariable("PSurf","float32",("time","z","y","x",))
      fld    = nc_file.createVariable("LWdown","float32",("time","y","x",))
      fsd_qc = nc_file.createVariable("SWdown_qc","int32",("time","y","x",))
      tair_qc= nc_file.createVariable("Tair_qc","int32",("time","y","x",))
      precip_qc = nc_file.createVariable("Rainf_qc","int32",("time","y","x",))
      qair_qc= nc_file.createVariable("Qair_qc","int32",("time","y","x",))
      ws_qc  = nc_file.createVariable("Wind_qc","int32",("time","y","x",))
      ps_qc  = nc_file.createVariable("PSurf_qc","int32",("time","y","x",))
      fld_qc = nc_file.createVariable("LWdown_qc","int32",("time","y","x",))
      elev = nc_file.createVariable("elevation","float32",("y","x",))
      ref_h    = nc_file.createVariable("reference_height","float32",("y","x",))
      if multiple_tiles:
        iveg_val   = nc_file.createVariable("iveg","int32",("patch","y","x",))
        patch_frac = nc_file.createVariable("patchfrac","float32",("patch","y","x",))
        LAI_val    = nc_file.createVariable("LAI","float32",("time","patch","y","x",))  # important: time must be first dimension!!
      else:
        iveg_val   = nc_file.createVariable("iveg","int32",("y","x",))
        patch_frac = nc_file.createVariable("patchfrac","float32",("y","x",))
        LAI_val    = nc_file.createVariable("LAI","float32",("time","y","x",))
      
      
      
      # 7.3) write data as a variable instance after the variable has been created
      ## requires L6 data to be loaded (can add as a df using Peter's tools)
        
      # append missing start and end values (if applicable), for now, just copy the first and/or last one
      if end_added_one:
        l6_Fsd         = numpy.append(l6_Fsd,l6_Fsd[-1])
        l6_Fsd_flag    = numpy.append(l6_Fsd_flag,l6_Fsd_flag[-1])
        l6_Ta          = numpy.append(l6_Ta,l6_Ta[-1])
        l6_Ta_flag     = numpy.append(l6_Ta_flag,l6_Ta_flag[-1])
        l6_Precip      = numpy.append(l6_Precip,l6_Precip[-1])
        l6_Precip_flag = numpy.append(l6_Precip_flag,l6_Precip_flag[-1])
        l6_q           = numpy.append(l6_q,l6_q[-1])
        l6_q_flag      = numpy.append(l6_q_flag,l6_q_flag[-1])
        l6_Ws          = numpy.append(l6_Ws,l6_Ws[-1])
        l6_Ws_flag     = numpy.append(l6_Ws_flag,l6_Ws_flag[-1])
        l6_ps          = numpy.append(l6_ps,l6_ps[-1])
        l6_ps_flag     = numpy.append(l6_ps_flag,l6_ps_flag[-1])
        l6_Fld         = numpy.append(l6_Fld,l6_Fld[-1])
        l6_Fld_flag    = numpy.append(l6_Fld_flag,l6_Fld_flag[-1])
      
      if end_added_two:
        l6_Fsd         = numpy.append(l6_Fsd,l6_Fsd[-1])
        l6_Fsd_flag    = numpy.append(l6_Fsd_flag,l6_Fsd_flag[-1])
        l6_Ta          = numpy.append(l6_Ta,l6_Ta[-1])
        l6_Ta_flag     = numpy.append(l6_Ta_flag,l6_Ta_flag[-1])
        l6_Precip      = numpy.append(l6_Precip,l6_Precip[-1])
        l6_Precip_flag = numpy.append(l6_Precip_flag,l6_Precip_flag[-1])
        l6_q           = numpy.append(l6_q,l6_q[-1])
        l6_q_flag      = numpy.append(l6_q_flag,l6_q_flag[-1])
        l6_Ws          = numpy.append(l6_Ws,l6_Ws[-1])
        l6_Ws_flag     = numpy.append(l6_Ws_flag,l6_Ws_flag[-1])
        l6_ps          = numpy.append(l6_ps,l6_ps[-1])
        l6_ps_flag     = numpy.append(l6_ps_flag,l6_ps_flag[-1])
        l6_Fld         = numpy.append(l6_Fld,l6_Fld[-1])
        l6_Fld_flag    = numpy.append(l6_Fld_flag,l6_Fld_flag[-1])
      
      
      if start_added_one:
        l6_Fsd         = numpy.append(l6_Fsd[0],l6_Fsd)
        l6_Fsd_flag    = numpy.append(l6_Fsd_flag[0],l6_Fsd_flag)
        l6_Ta          = numpy.append(l6_Ta[0],l6_Ta)
        l6_Ta_flag     = numpy.append(l6_Ta_flag[0],l6_Ta_flag)
        l6_Precip      = numpy.append(l6_Precip[0],l6_Precip)
        l6_Precip_flag = numpy.append(l6_Precip_flag[0],l6_Precip_flag)
        l6_q           = numpy.append(l6_q[0],l6_q)
        l6_q_flag      = numpy.append(l6_q_flag[0],l6_q_flag)
        l6_Ws          = numpy.append(l6_Ws[0],l6_Ws)
        l6_Ws_flag     = numpy.append(l6_Ws_flag[0],l6_Ws_flag)
        l6_ps          = numpy.append(l6_ps[0],l6_ps)
        l6_ps_flag     = numpy.append(l6_ps_flag[0],l6_ps_flag)
        l6_Fld         = numpy.append(l6_Fld[0],l6_Fld)
        l6_Fld_flag    = numpy.append(l6_Fld_flag[0],l6_Fld_flag)

    
  
      # 7.4) fill nc with values
      times[:] = sec_l6_time   # note the change to sec from days (note that 'valid' statement was made before already)
      fsd[:]   = l6_Fsd[valid]
      fsd_qc[:]= l6_Fsd_flag[valid]
      tair[:]  = l6_Ta[valid] + Kelvin   # degC to K
      tair_qc  = l6_Ta_flag[valid]
      precip[:]= l6_Precip[valid]/timestep   # mm to mm/s
      precip_qc[:]=l6_Precip_flag[valid]
      qair[:]  = l6_q[valid]
      qair_qc  = l6_q_flag[valid]
      ws[:]    = l6_Ws[valid]
      ws_qc[:] = l6_Ws_flag[valid]
      ps[:]    = l6_ps[valid] * 1000.    # kPa to Pa
      ps_qc[:] = l6_ps_flag[valid]
      fld[:]   = l6_Fld[valid]
      fld_qc[:]= l6_Fld_flag[valid]
      elev[:]  = elevation[s]
      x_val[:] = 1.0
      y_val[:] = 1.0
      z_val[:] = 1.0
      if generate_site_list:
        lat[:] = latitude[-1]
        lon[:] = longitude[-1]
      else:
        if site == 'Gingin':
          lat[:] = -31.375
          lon[:] = 115.65
        else:
          lat[:] = float(ds.globalattributes['latitude'])
          lon[:] = float(ds.globalattributes['longitude'])
      try:  
        ref_h[:] = sitelist.loc[site,'reference_height']  
      except KeyError:
        pass
        
      # LAI_val[:]   = np.ones(ts_len)*LAI[s]
      
     
      NVIS5_class = sitelist.loc[site,'NVIS5_Group']
      C4_fract = sitelist.loc[site,'C4_fraction']
          
      if C4_fract >= 0.5:
          grass_patch = 7
      else:
          grass_patch = 6
          
      if NVIS5_class == 7:
          forest_patch = 1
      else:
          forest_patch = 2  

        
        
      ## LAI (Trudinger et al. 2016, Biogeosciences)
      # 1) load fPAR data (values interpreted to represent the middle of the month)
      site_nr = fcover_sites.index(site)
        
      fw = np.array(fcover.iloc[fper_index,site_nr])
      fg = np.array(fcover.iloc[frec_index,site_nr])
          
        
      # 2) calculate monthly LAI (Trudinger et al. 2016, Biogeosciences)
      LAIw_monthly = -1/k * np.log(1-fw)
      LAIg_monthly = -1/k * np.log(numpy.maximum(1-(fg/(1 - fg)), 0.01))   # boundaries need to be better defined here!
          
          
      # 3) interpolate to halfhourly/hourly values
      ## for cubic interpolation, add three outer months to cut them later
      LAIw_monthly2 = np.append(LAIw_monthly[-3:],LAIw_monthly)
      LAIw_monthly2 = np.append(LAIw_monthly2,LAIw_monthly[0:3])
      LAIg_monthly2 = np.append(LAIg_monthly[-3:],LAIg_monthly)
      LAIg_monthly2 = np.append(LAIg_monthly2,LAIg_monthly[0:3])
        
      LAIw = []
      LAIg = []
          
      years = range(int(str_startyear),int(str_endyear)+1)   # why +1 here??
      for yr in years:
      # print(yr)
        if calendar.isleap(yr):
            dpm = np.array([31,30,31,31,29,31,30,31,30,31,31,30,31,30,31,31,29,31])   # Oct - March over three years 
        else:
            dpm = np.array([31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,31,28,31])   # Oct - March over three years
              
        x = (numpy.cumsum(dpm)-15.5) * 24 * 3600 
        x2= np.arange(x[0],x[-1],timestep) 
        
        fun_interp_w = interp1d(x,LAIw_monthly2,kind='cubic')     # gives the function to interpolate
        fun_interp_g = interp1d(x,LAIg_monthly2,kind='cubic')
            
        LAIw_timestep = fun_interp_w(x2)
        LAIg_timestep = fun_interp_g(x2)
        
        # avoid negative values (potentially caused by interpolation)
        LAIw_timestep = np.clip(LAIw_timestep,0.01,12)
        LAIg_timestep = np.clip(LAIg_timestep,0.01,12)
            
        # now cut the additional months again and concatenate vectors
        first = (numpy.cumsum(dpm)[3]-15.5)*24*3600
        last = (numpy.cumsum(dpm)[-3]-15.5)*24*3600
            
        LAIw_timestep = LAIw_timestep[(x2 >= first) & (x2 < last)]
        LAIg_timestep = LAIg_timestep[(x2 >= first) & (x2 < last)]
            
        LAIw = numpy.append(LAIw,LAIw_timestep)
        LAIg = numpy.append(LAIg,LAIg_timestep)
          
        #LAI_total = LAIw + LAIg  # incorrect
        frac_forest0 = sitelist.loc[site,'forest_fraction']
        frac_grass0 = sitelist.loc[site,'grass_fraction']
          
        frac_forest = min(frac_forest0/fcorr_patch,1.0)
        frac_grass  = 1.0 - frac_forest 
          
        LAI_total = LAIw * frac_forest + LAIg * frac_grass
          
        #plt.plot(x,LAIg_monthly2,'o',x2,fun_interp_g(x2),'-')
          
      if multiple_tiles:
          LAI_val[:,0,:,:] = LAIw
          LAI_val[:,1,:,:] = LAIg
            
          ## iveg
          iveg_val[:,:] = [forest_patch,grass_patch]
          
          ## patchfrac
          patch_frac[0] = frac_forest
          patch_frac[1] = frac_grass
          
          
      else:    # add LAI only (one patch only)
          LAI_val[:,:,:] = LAI_total
          iveg_val[:] = [forest_patch]
          patch_frac[:] = 1.0
      
      
        # retreieve non ts values
        #src = Dataset(os.path.join(i_path, i_file), "r")
        
      ## 7.5) copy and modify variable attributes
      nc_file.variables['time'].setncatts(time_attr)
      nc_file.variables['SWdown'].setncatts(Fsd_attr)
      nc_file.variables['SWdown_qc'].setncatts(Fsd_attr)
      nc_file.variables['Tair'].setncatts(Ta_attr)
      nc_file.variables['Tair_qc'].setncatts(Ta_attr)
      nc_file.variables['Rainf'].setncatts(Precip_attr)
      nc_file.variables['Rainf_qc'].setncatts(Precip_attr)
      nc_file.variables['Qair'].setncatts(q_attr)
      nc_file.variables['Qair_qc'].setncatts(q_attr)
      nc_file.variables['Wind'].setncatts(Ws_attr)
      nc_file.variables['Wind_qc'].setncatts(Ws_attr)
      nc_file.variables['PSurf'].setncatts(ps_attr)
      nc_file.variables['PSurf_qc'].setncatts(ps_attr)
      nc_file.variables['LWdown'].setncatts(Fld_attr)
      nc_file.variables['LWdown_qc'].setncatts(Fld_attr)
      
      # modify time attribute
      #nc_file.variables['time'].units='seconds since 2017-01-01 00:00:00.0'
      
      
      ### shift time step one forward because JSBACH marks the beginning of the time step, not the end
      startdate = startdate - datetime.timedelta(seconds=timestep)
      
      
      date_names = ('month','day','hour','minute')
      date_elements =  {}
      for i in date_names:
        name  = 's' + i
        value = str(getattr(startdate,i))
        if len(value) < 2:
          value = '0' + value
        
        date_elements[name]=value
        
#        globals()[name] = value
      

      nc_file.variables['time'].units='seconds since ' + str(startdate.year) + '-' + date_elements['smonth'] + '-' + date_elements['sday'] + \
                                      ' ' + date_elements['shour']  + ':' + date_elements['sminute'] + ':00'
      del nc_file.variables['time'].calendar
        
      # modify other variables' units
      nc_file.variables['elevation'].units='m'
      nc_file.variables['reference_height'].units='m'
      nc_file.variables['latitude'].units='degrees_north'
      nc_file.variables['longitude'].units='degrees_east'
      nc_file.variables['Rainf'].units='mm/s'
      nc_file.variables['PSurf'].units='Pa'
      nc_file.variables['Tair'].units='K'
      nc_file.variables['LWdown_qc'].units='-'  # not sure if necessary, just in case
      nc_file.variables['PSurf_qc'].units='-'
      nc_file.variables['Wind_qc'].units='-'
      nc_file.variables['Qair_qc'].units='-'
      nc_file.variables['Rainf_qc'].units='-'
      nc_file.variables['Tair_qc'].units='-'
      nc_file.variables['SWdown_qc'].units='-'
      nc_file.variables['x'].units=''
      nc_file.variables['y'].units=''
      nc_file.variables['z'].units=''
  
  
      # batch
      # src = Dataset(os.path.join(i_path, i_file), "r")
        
      # for name, variable in src.variables.items():
      #     if name in rename_zip_list:
        
      #         # copy variable attributes all at once via dictionary
      #         ww_nc[name].setncatts(src[name].__dict__)
        
      ## 8.6) close nc file
      nc_file.close()

## end site loop
  
  
  
  
  
### create dataframe with basic site properties
if generate_site_list:
  site_attr_df = pd.DataFrame(index=data_sites)
  #site_attr_df['site'] = data_sites
  site_attr_df['startyear'] = startyear 
  site_attr_df['endyear'] = endyear
  site_attr_df['latitude'] = latitude
  site_attr_df['longitude'] = longitude
  site_attr_df['canopy_height'] = canopy_height
  site_attr_df['tower_height'] = tower_height
  site_attr_df['vegetation'] = vegetation


  ### export dataframe to file
  site_attr_df.to_csv('C:/Users/kna016/Documents/CABLE/OzFLUX_sitelist_v1.txt',sep="\t")





  ### combine site list with other information
  sitelist = pd.read_csv('C:/Users/kna016/Documents/CABLE/OzFLUX_sitelist_v1.txt',sep="\t",index_col=0) # as created above

  ## NVIS5 vegetation classification
  NVIS5 = pd.read_csv('C:/Users/kna016/Documents/CABLE/cover_fract/ozflux30.nvis5_group_type.csv',sep=",")
  NVIS5_Group = NVIS5.iloc[4,1:]


  sitelist = sitelist.join(NVIS5_Group,how='left')
  sitelist.rename(columns={4:'NVIS5_Group'},inplace=True)
  sitelist['NVIS5_Group'] = sitelist['NVIS5_Group'].astype('int')


  sitelist.to_csv('C:/Users/kna016/Documents/CABLE/OzFLUX_sitelist_v2.txt',sep="\t")



#NVIS5.info  # pandas.DataFrame.info
#NVIS5.index


## In[199]:
#
#
## map output var names to L6 var names (doesn't include x,y,z,elevation,reference_height,LAI)
#in_var_list = ['latitude','longitude','time','Fsd','Ta',
#              'Precip','q','Ws','ps','Fld','Fsd_QCFlag','Ta_QCFlag','Precip_QCFlag',
#              'q_QCFlag','Ws_QCFlag','ps_QCFlag','Fld_QCFlag']
#
#in_var_list_ts = ['time','Fsd','Ta',
#              'Precip','q','Ws','ps','Fld','Fsd_QCFlag','Ta_QCFlag','Precip_QCFlag',
#              'q_QCFlag','Ws_QCFlag','ps_QCFlag','Fld_QCFlag']
#
## no qc
#in_var_list_ts_no_qc = ['time','Fsd','Ta','Precip','q','Ws','ps','Fld']
#
#
## In[182]:
#
#
#in_var_list
#
#
## In[189]:
#
#
#in_var_list_not_ts = ['latitude','longitude']
#
#
## In[200]:
#
#
#rename_zip_list = ['LWdown','LWdown_qc','SWdown','SWdown_qc','Rainf','Rainf_qc','Tair','Tair_qc',
#                  'Wind','Wind_qc','latitude','longitude','PSurf','PSurf_qc','Qair','Qair_qc','time']
#
#rename_zip_list_ts = ['LWdown','LWdown_qc','SWdown','SWdown_qc','Rainf','Rainf_qc','Tair','Tair_qc',
#                  'Wind','Wind_qc','PSurf','PSurf_qc','Qair','Qair_qc','time']
#
## does not include qc vars as these are read in with the data
#rename_zip_var_list_ts = ['LWdown','SWdown','Rainf','Tair','Wind','PSurf','Qair','time']
#
#zip_list = list(zip(sorted(in_var_list_ts_no_qc),rename_zip_var_list_ts))
#zip_list






