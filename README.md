# ROMS_3.7_for_ETRAC
Regional Ocean Modeling System (ROMS), version 3.7, with added capability to produce output for the ETRAC software

This ROMS version is based on ROMS 3.7, revision 854, from the official ROMS SVN repository, accessible for all registered ROMS users: http://myroms.org/

A new type of diagnostic output for the nitrogen cycle was implemented, which allows for "perfect" mass balances of nitrogen state variables (within the order of output data precision).
This output is required as input for the Element TRACing Software (ETRAC) by Fabian Grosse

These CPP flags need to be set to be able to enable the new output:
SOLVE3D
BIO_FENNEL
DENITRIFICATION
BIO_SEDIMENT

The new output is enabled setting the following new CPP flags:
TBNT_OUT       => turn on the new output ("TBNT" refers to the term "Trans-Boundary Nutrient Transports"; see, e.g. Blauw et al (2006))
PSRC_UV_BUGFIX => bug fix for avoiding volume transports onto land at horizontal land-sea interfaces with point sources; see ROMS forum post: https://www.myroms.org/forum/viewtopic.php?f=19&t=4643

In addition, it is recommended to enable the following switches:
TBNT_NFDOUBLE   => double precision output for new output
PERFECT_RESTART => allows for perfect continuation of a simulation after a restart (implies that information for 2 time steps are stored in restart files instead of only 1 (default))

Further additions are required in the following ROMS input files (search for "TBNT" in the listed files):
bio_Fennel.in => see /User/External/bio_Fennel_mch078_TBNT.in
varinfo.dat   => see /ROMS/External/varinfo_Fennel_PO4_RDON_TBNT.dat


Besides this new output, the present ROMS version contains several additions implemented by members of the group of Prof. Katja Fennel at Dalhousie University, e.g. phosporus limitation in the bio_fennel biology etc.


References:
Blauw, A., van deWolfshaar, K., andMeuwese, H. (2006). Transboundary nutrient
transports in the North Sea. WL|Delft Hydraulics Reports, Z4188.
url: http://publicaties.minienm.nl/documenten/transboundary-nutrient-transports-in-the-north-sea-model-study
