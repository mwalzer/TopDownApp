## mzTab
The TopDownApp handles input and output data with HUPO-PSI standard formats. Development code, dedicated im-/exporters, and documentation can be found here. 

The export to mzTab is mainly achieved through prior dataframe organisation and pandas tsv export. Note that for this workflow purpose, the mzTab header is kept (almost) static. It is demonstrating how a minimal "Summary Top-Down Proteoform Identification Report" can look like.
The main deviation from the current standard is the use of a PRSM section instead of a PSM section, as PrSM data will have slightly differnet implications and columns that minimally should be reported. 

Main points of consideration for new tools that add to the data mzTab will contain is the availability of appropriate cv terms like `MTD software[2] [MS, 1002714, TOPP FLASHDeconv, ]` or `[MS, MS:1002928, TopPIC:spectral E-value, ]` that need to be added/replaced in the header template.
New columns can be added with adherence to the current mzTab standard by adding "opt_<column_name>" to the column mapping.

