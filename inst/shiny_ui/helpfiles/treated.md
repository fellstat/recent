# Treated HIV Case

A binary variable representing whether each subject is on ARV therapy. Treatment status should be determined by ARV biomarkers. If any of the biomarkers are missing, individuals who self-report as on treatment should be considered 'treated'.

This variable is only required or used if the RITA2 algorithm is desired. By default, RITA3 is performed.

Values can be logical (TRUE or FALSE), string valued, or numeric (1, 0). The largest value is assumed to indicate TRUE if numeric, and the last value in an alpha-numeric sorting indicates TRUE if it is a character string. Missing values should be encoded as NA in the data.
