
PGATK File Formats
=====================

The ProteoGenomics Analysis ToolKit is based on standard proteomics formats developed by `HUPO-PSI <https://github.com/HUPO-PSI>`_ and Genomics Standard file formats. This section is to highlight in 10 minutes the most important features of those file formats, How they are used in PGATK and you can contribute to their development.

.. note:: It is important to notice that this Help page only provides the fundamentals of each file format used in PGATK, for more details we provide links to the original documentation of the File format.

.. _bed:

BED
-------------------

BED ***(Browser Extensible Data)** format provides a flexible way to define the data lines that are displayed in an annotation track `UCSC Bed Definition <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_. BED lines have three required fields and nine additional optional fields. The number of fields per line must be consistent throughout any single set of data in an annotation track. The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used.



