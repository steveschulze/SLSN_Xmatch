Tool to crossmatch the output from the scanning page on the GROWTH Marshal with the catalogues in the VIZIER database. The following catalogues are included: GAIA DR2, SDSS DR12 and Milliquas.

-Candidates are identified as stars if the entry in the GAIA DR2 has either a parallax or a proper-motion measurement and the significance >3 sigma. Furthermore, candidates are rejected if they are too close to a star (to ensure good observability and avoid artefacts). This criterion is almost identical to that in Hjorth et al. (2012).

-Candidates are identified as a quasar if the SDSS subclass includes the keywords "AGN" or "BROADLINE".

-No filtering is done based on the Milliquas catalogue because it also includes QSO candidates. Check out the Milliquas documentation https://heasarc.gsfc.nasa.gov/W3Browse/all/milliquas.html.

Requirements: astropy, astroquery, python 3.6 (but it should also work in python 2.7)
