Introduction
============

.. sidebar:: Individual Tools Documentation
   :subtitle: **It can make your life easier** if you want to explore individual tools:

   * PepGenome: Mapping Peptidoforms to Genome Coordinates
   * PepBed: R package to analyze peptide features in Bed Files


The PRIDE Archive API provides a Restful entry point to retrieve all data from PRIDE Archive database. The API provides hyperlinks the clients can follow to navigate and search
in the resource. Just like a human user of a regular website, who knows the initial URL of a website and then follows hyperlinks to navigate through the site.

A simple call to PRIDE Archive API docs will provide all the entry points in the web service:

.. http:example:: curl python-requests

   GET api-docs HTTP/1.1
   Host: www.ebi.ac.uk/pride/ws/archive/
   Accept: application/json