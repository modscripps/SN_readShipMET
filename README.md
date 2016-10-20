SN_readShipMET
==============

SN_READSHIPMET - reads SIO ships MET data files<br/>
<br/>
SN_READSHIPMET(FILENAME) returns a MET structure of variables described for MET data files<br/>
<br/>
SN_READSHIPMET(DIRNAME) reads all *.MET files in the directory DIRNAME<br/>
<br/>
SN_READSHIPMET({FILE1, FILE2,...}) reads all files indicated<br/>
<br/>

Created 2012/03/30 San Nguyen<br/>
Updated 2012/06/09 San Nguyen - updated to read MET data for most ships not just the ones on the R/V Revelle<br/>
Updated 2015/09/02 San Nguyen - update to combine MET data without losing any fields added later on<br/>
Updated 2016/10/19 Emily Shroyer- Corrected Reanme Information<br/>
Updated 2016/10/20 San Nguyen - Update to read MET files in the same folder to preserve maximum number of fields that might have been added after the initial recording, not restricting number of fields to defined by the first file.<br/>