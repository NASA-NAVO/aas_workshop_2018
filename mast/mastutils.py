#
# Imports
#

import requests
import io
import numpy as np
from astropy.table import Table, unique

#
# Support for VOTABLEs as astropy tables
#

def astropy_table_from_votable_response(response):
    """
    Takes a VOTABLE response from a web service and returns an astropy table.
    
    Parameters
    ----------
    response : requests.Response
        Response whose contents are assumed to be a VOTABLE.
        
    Returns
    -------
    astropy.table.Table
        Astropy Table containing the data from the first TABLE in the VOTABLE.
    """
    
    # The astropy table reader would like a file-like object, so convert
    # the response content a byte stream.  This assumes Python 3.x.
    # 
    # (The reader also accepts just a string, but that seems to have two 
    # problems:  It looks for newlines to see if the string is itself a table,
    # and we need to support unicode content.)
    file_like_content = io.BytesIO(response.content)
    
    # The astropy table reader will auto-detect that the content is a VOTABLE
    # and parse it appropriately.
    aptable = Table.read(file_like_content)
    
    # String values in the VOTABLE are stored in the astropy Table as bytes instead 
    # of strings.  To makes accessing them more convenient, we will convert all those
    # bytes values to strings.
    stringify_table(aptable)
    
    return aptable

def find_column_by_ucd(table, ucd):
    """
    Given an astropy table derived from a VOTABLE, this function returns
    the first Column object that has the given Universal Content Descriptor (UCD).  
    The name field of the returned value contains the column name and can be used 
    for accessing the values in the column.
    
    Parameters
    ----------
    table : astropy.table.Table
        Astropy Table which was created from a VOTABLE (as if by astropy_table_from_votable_response).
    ucd : str
        The UCD identifying the column to be found.
    
    Returns
    -------
    astropy.table.Column
        The first column found which had the given ucd.  None is no such column found.
        
    Example
    -------
    col = find_column_by_ucd(my_table, 'VOX:Image_Title')
    print ('1st row title value is:', my_table[col.name][0])
    """
    
    # Loop through all the columns looking for the UCD
    for key in table.columns:
        col = table.columns[key]
        ucdval = col.meta.get('ucd')
        if (ucdval is not None):
            if (ucd == ucdval):
                return col
    
    return None

#
# Wrappers for Virtual Observatory queries
#

def tap_query(query, tap_service):
    """
    Executes the specified TAP query at the specified TAP service and returns
    the result as an astropy Table.  This function assumes that tap_service is
    a synchronous TAP service URL.
    
    Parameters
    ----------
    query : str
        ADQL (SQL-like) to execute on the given TAP service.
    tap_service : str
        URL endpoint for the desired synchronous TAP service.  This URL will be
        the endpoint for the TAP service identified in the NAVO registry, followed by '/sync'.
    
    Returns
    -------
    astropy.table.Table
        An astropy table created from the contents of the VOTABLE response from the TAP query.
           
    """
    tap_params = {
        "request": "doQuery",
        "lang": "ADQL",
        "query": query
    }
    results = requests.post(tap_service, data=tap_params)
    aptable = astropy_table_from_votable_response(results)
    return aptable

def sia_query(access_url, ra, dec, size, image_formats='ALL', print_url=False):
    """
        
    Parameters
    ----------
    access_url : str
        The URL endpoint for the SIA service.
    ra : str or number
        The Right Ascension (ICRS degrees) of the central search point
    dec: str or number
        The Declination (ICRS degrees) of the central search point
    size: str or number
        The diameter of the search circle (cone) in degrees
    image_formats: str
        A comma-delimited list of acceptable image MIME types 
        (choices are usually 'image/fits', 'image/jpeg', 'image/png').
        May also choose:
            'ALL' - All formats acceptable (the default)
            'GRAPHIC' - Allows jpeg png and gif
            'METADATA' - Denotes a Metadata Query: no images are requested; only metadata 
            (i.e., SIA result FIELD definitions) should be returned.
    print_url : bool
        Whether or not to print the query URL.  
    
    Returns
    -------
    astropy.table.Table
        An astropy table created from the contents of the VOTABLE response from the SIA query.
    """
    params = {
        'POS': sval(ra) + ',' + sval(dec),
        'SIZE': sval(size),
        'FORMAT': image_formats
    }
    results = requests.get(access_url, params)
    if (print_url):
        print('Queried: ' + results.url)
    aptable = astropy_table_from_votable_response(results)
    return aptable

#
# Functions to help replace bytes with strings in astropy tables that came from VOTABLEs
#

def sval(val):
    """
    Returns a string value for the given object.  When the object is an instanceof bytes,
    utf-8 decoding is used.
    
    Parameters
    ----------
    val : object
        The object to convert
    
    Returns
    -------
    string
        The input value converted (if needed) to a string
    """
    if (isinstance(val, bytes)):
        return str(val, 'utf-8')
    else:
        return str(val)

# Create a version of sval() that operates on a whole column.
svalv = np.vectorize(sval)

def sval_whole_column(single_column):
    """
    Returns a new column whose values are the string versions of the values
    in the input column.  The new column also keeps the metadata from the input column.
    
    Parameters
    ----------
    single_column : astropy.table.Column
        The input column to stringify
        
    Returns
    -------
    astropy.table.Column
        Stringified version of input column
    """
    new_col = svalv(single_column)
    new_col.meta = single_column.meta
    return new_col

def stringify_table(t):
    """
    Substitutes strings for bytes values in the given table.
    
    Parameters
    ----------
    t : astropy.table.Table
        An astropy table assumed to have been created from a VOTABLE.
    
    Returns
    -------
    astropy.table.Table
        The same table as input, but with bytes-valued cells replaced by strings.
    """
    # This mess will look for columns that should be strings and convert them.
    if (len(t) is 0):
        return   # Nothing to convert
    
    scols = []
    for col in t.columns:
        colobj = t.columns[col]
        if (colobj.dtype == 'object' and isinstance(t[colobj.name][0], bytes)):
            scols.append(colobj.name)

    for colname in scols:
        t[colname] = sval_whole_column(t[colname])
        

