# projections
conversion tools for Lambert conformal projections used with reanalysis datasets
Paul Leahy UCC 2017,2019

To use:
From the Matlab command prompt type:

p=projection(‘mera’);

Edit the script convert_lcc_1sp_latlon_to_EN.m and replace the test location’s latitude and longitude co-ordinates  (“Ballyshannon” currently) to the co-ordinates you want, in degrees.

Run convert_lcc_1sp_latlon_to_EN.m

The x,y values in the MERA grid will be printed.

There is also a similar script for the reverse conversion, convert_lcc_1sp_EN_to_latlon.m

See also:
http://eel.ucc.ie/2017/06/06/handling-data-lambert-conformal-conic-projections/

