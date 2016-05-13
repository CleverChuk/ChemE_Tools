# CHEM_E_TOOLS
Programs for chemical engineering computations.

This is the same PreoProps furnished by Carl Lira for Introductory to Chemical Engineering Thermodynamics;

Recent Fixes By Clever

Added Dependencies:
            You'll need to have the queryFunc for this program to run without errors
            You'll need to run PreosPropsMenu.m to properly initialize the UI

Fixes:
    Fixed the error that occurs when userObject function is used to guess TnP
            
This is the current version of the programs for 
"Introductory Chemical Engineering Thermodynamics"
by J.R. Elliott, Jr., and C.T. Lira
described at http://chethermo.net

ABOUT THE PROGRAM LOCATIONS IN THE DIRECTORIES

The programs are contained in subdirectories by classification. 
The files are explained briefly in appendix A of the textbook, 
or in the file appendixA.pdf which is an online version
of appendix A that can be viewed with the free Adobe Acrobat reader 
availabe from http://get.adobe.com/reader.

	Excel- Spreadsheet files
	Matlab - Matlab files

Note that an online supplement is available for getting started with
MATLAB and Excel Solver, 
http://chethermo.net/sites/default/files/doc/supp/SuppMatlabSolvOptim.pdf
http://chethermo.net/sites/default/files/doc/supp/SuppExcel.pdf

UPDATE INFORMATION

A summary of the most recent updates is available at https://sourceforge.net/p/chethermo/code/commit_browser. Beginning with release 2.25, the digits after the decimal indicate the revision commit.

2.05 12/1/12 - Matlab/Props/props.mat has been updated. When pasting from propx.xlsx, Cp values were truncated to the fixed format of the display. This created slight errors and slight differences between the spreadsheets and Matlab. This has been fixed by changing to general format before copy/paste.

2.04 6/5/12 - Matlab changes. Using struct variables to load Matlab databases, also some minor improvements in various Matlab codes. Added Chap15/PRfug.m. Updated Matlab QuickRef doc.

2.03 5/9/12 - Matlab/Chap17/Kcalc.m was improved for handling cases with missing data. Better printing of property values.

2.02 4/13/12 - Excel/Props.xls and Matlab/Props/props.mat were updated. Matlab/Chap17 files were added. A bug in Matlab/Chap17/Kcalc.m was fixed. The Matlab QuickRef was updated. Other minor changes.

2.01 3/15/12 - The Matlab/Chap08-09 folder has been revised to provide for matching volume.

2.00 1/15/2012 PRMIX.exe and ESD.exe were discontinued. Additional files were added for 
second edition.

1/22/2011 Files were updated for consistency with revised chapter/example numbers.

1/23/2010 Matlab folders where added. The HP and TI programs were discontinued. The PRPURE content was discontinued. Excel/Rxns.xls was updated to use full Cp polynomials for adiabatic reactors. Preos.xls was updated to use the latest value of R from NIST.

8/11/05 Several updates were included:
        1. Residue.xls reference to Solver.xla was changed to one that is more
           common. This will save some users time in specifying the directory
           for Solver.xla.
        2. Steam.xls was updated to correct errors located by user Edward D. Throm.
        3. HP programs were updated to be compatible with the HP48gII.
        4. The appendix was updated to reflect the changes in the HP programs.

2/28/05 The ftp site was discontinued. The files are now available from the web server.

11/12/04 The readme files for the FORTRAN files were improved to provide better instructions to avoid the 'disappearing window' when the program ends.

2/16/04 An error was fixed in prmix.exe in the calculation of the entropy departure
for mixtures.

1/10/00 Apxcomu.pdf was updated to match the version in the text (only formatting
changes). The following spreadsheets were updated improving formula usage of
named ranges and labels: Actcoeff.xls, Gammafit.xls, Preos.xls, Dewcalc.xls.
In addition, functionality of UNIFAC(vle) was increased on Actcoeff.xls to
permit some flexibility in groups listed in the table.
