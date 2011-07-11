This is a Python add-on module to read in RIF data.
It contains two functions written in ANSI C and wrapped with f2py.
To install you need to have Python, Numpy and MinGW (C/C++ and FORTRAN compilers) installed.
You also need to have the MinGW bin directory (usually C:\MinGW\bin) added to your PATH variable.

To compile the module (against Python 2.6) run:

C:\RIF_DAQ\Python\rif>C:\Python26\Scripts\f2py.py -c --fcompiler=gnu95 --compiler=mingw32 -lmsvcr90 -IC:\RIF_DAQ\Python\rif -LC:\RIF_DAQ\Python\rif -lfftw3-3 _rif.pyf _rif.c

and copy rif.pyd to a location where Python will find it (possibly the same directory containing
your data but preferably your local site-packages directory.

Alternatively you can set the PYTHONPATH variable to include the directory used for building the module.
On windows this is done by:

 To augment PYTHONPATH, run regedit and navigate to  HKEY_LOCAL_MACHINE\SOFTWARE\Python\PythonCore  and
then select the folder for the python version you wish to use. Inside this is a folder labelled PythonPath,
with one entry that specifies the paths where the default install stores modules. Right-click on PythonPath
and choose to create a new key. You may want to name the key after the project whose module locations it will
specify; this way, you can easily compartmentalize and track your path modifications.
Your new key will have one string value entry, named (Default). Right-click on it and modify its value data;
this should be text in the same format as the Path environment variable discussed above--absolute directory paths,
separated by semicolons. If one project will use modules from several directories, add them all to this list.
(Don't bother attempting to add more string value entries to your new key, or to the original PythonPath key, since they will be ignored.)
Once these new registry entries are in place, your scripts' import statements should work fine. 
