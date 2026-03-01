
"""
This module contains Ensembl models used as Entities and DAO/Services
"""

import json


class Species:
    def __init__(self, ensembl_species_entry):
        self.__ensembl_species_entry = ensembl_species_entry

    def get_logger(self):
        return self.__logger

    def get_value_for_key_or_default(self, key, default_value='-not_available-'):
        if key in self.get_ensembl_species_entry():
            return self.get_ensembl_species_entry()[key]
        return default_value

    def get_ensembl_species_entry(self):
        return self.__ensembl_species_entry

    def get_division(self):
        return self._get_value_for_key_or_default('division')

    def get_ncbi_taxonomy_id(self):
        return self._get_value_for_key_or_default('taxon_id')

    def get_name(self):
        return self._get_value_for_key_or_default('name')

    def get_ensembl_release(self):
        return self._get_value_for_key_or_default('release')

    def get_display_name(self):
        return self._get_value_for_key_or_default('display_name')

    def get_assembly_accession(self):
        return self._get_value_for_key_or_default('accession')

    def get_strain_collection(self):
        if self._get_value_for_key_or_default('strain_collection', 'null') == 'null':
            return None
        return self._get_value_for_key_or_default('strain_collection')

    def get_common_name(self):
        return self._get_value_for_key_or_default('common_name')

    def get_strain(self):
        return self._get_value_for_key_or_default('strain')

    def get_aliases(self):
        return self._get_value_for_key_or_default('aliases', [])

    def get_groups(self):
        return self._get_value_for_key_or_default('groups', [])

    def get_assembly(self):
        return self._get_value_for_key_or_default('assembly')

    def __str__(self):
        return json.dumps(self.get_ensembl_species_entry())


class SpeciesService:
    def __init__(self, species_data):
        # I've changed this, we store the original species data, and then we offer two different views
        self.__ensembl_species_data_raw = species_data
        self.__ensembl_species_data_dao = None
        self.__index_by_taxonomy_id = None
        self.__index_by_assembly = None

    def __index_data_for_property(self, data, property_getter):
        """
        Given an iterable data container, and a property getter to run on every object of that container, it returns a
        dictionary where the key is the property value for a particular data object part of the data collection
        :param data: iterable of data objects to index
        :param property_getter: property on which the index should be created, this property is a member function to a
        getter for that property
        :return: a dictionary of the given data objects where the key is the indexed property
        """
        self._get_logger().debug("Creating index on getter '{}' for #{} entries".format(property_getter, len(data)))
        return {property_getter(data_item): data_item for data_item in data}

    def __index_data_by_taxonomy_ensembl_special_case(self, data):
        indexed_data = {}
        for data_item in data:
            if data_item.get_ncbi_taxonomy_id() not in indexed_data:
                indexed_data[data_item.get_ncbi_taxonomy_id()] = data_item
            else:
                # Keep the one with aliases, that seems to be the main one for Ensembl
                if data_item.get_aliases():
                    if indexed_data[data_item.get_ncbi_taxonomy_id()].get_aliases():
                        self._get_logger().error("ENSEMBL SPECIES INDEXING ERROR, already existing species entry '{}' "
                                                 "will not be replaced by non-empty aliases entry '{}'"
                                                 .format(str(indexed_data[data_item.get_ncbi_taxonomy_id()]),
                                                         str(data_item)))
                    else:
                        self._get_logger().warning("ENSEMBL SPECIES INDEX REPLACEMENT, of existing entry '{}' "
                                                   "by new entry '{}'"
                                                   .format(str(indexed_data[data_item.get_ncbi_taxonomy_id()]),
                                                           str(data_item)))
                        indexed_data[data_item.get_ncbi_taxonomy_id()] = data_item
                else:
                    self._get_logger().warning("ALREADY INDEXED ENSEMBL SPECIES ENTRY '{}' "
                                               "WILL NOT BE REPLACED by '{}'"
                                               .format(str(indexed_data[data_item.get_ncbi_taxonomy_id()]),
                                                       str(data_item)))
        return indexed_data

    def _get_logger(self):
        return self.__logger

    def _get_species_data_dao(self):
        """
        This method creates the DAO version of the ensembl species data, it could have been implemented as a generator
        but I want those DAO objects to be reused by the Python reference system, so having them in multiple indexes
        does not require multiple times the same amount of memory
        :return: a list of DAO accessors (Species) to the given RAW ensembl species data
        """
        if self.__ensembl_species_data_dao is None:
            self.__ensembl_species_data_dao = \
                [Species(species_entry) for species_entry in self.get_species_data()['species']]
        return self.__ensembl_species_data_dao

    def _get_index_taxonomy_id(self):
        """
        Build the index for species data by taxonomy ID
        :return: ensembl species data indexed by taxonomy ID
        """
        if self.__index_by_taxonomy_id is None:
            self.__index_by_taxonomy_id = \
                self.__index_data_by_taxonomy_ensembl_special_case(self._get_species_data_dao())
            # Generic and beautiful old way of building the index
            # self.__index_data_for_property(self._get_species_data_dao(), Species.get_ncbi_taxonomy_id)
        return self.__index_by_taxonomy_id

    def _get_index_assembly(self):
        """
        Build the index for species data by Assembly
        :return: ensembl species data indexd by assembly
        """
        if not self.__index_by_assembly:
            self.__index_by_assembly = self.__index_data_for_property(self._get_species_data_dao(),
                                                                      Species.get_assembly)
        return self.__index_by_assembly

    def get_species_data(self):
        """
        Get the Ensembl species data unprocessed, as it came from the Ensembl REST API
        :return: raw ensembl species data
        """
        return self.__ensembl_species_data_raw

    def get_species_entry_for_taxonomy_id(self, taxonomy_id):
        """
        Given a taxonomy ID, get its Ensembl species entry
        :param taxonomy_id: taxonomy ID
        :return: the species entry or None if not found
        """
        self._get_logger().debug("get_species_entry_for_taxonomy_id ---> taxonomy id '{}'".format(taxonomy_id))
        if taxonomy_id in self._get_index_taxonomy_id():
            return self._get_index_taxonomy_id()[taxonomy_id]
        return None

    def get_species_entry_for_assembly(self, assembly):
        """
        Given an Assembly, get its Ensembl species entry
        :param assembly: assembly
        :return: the species entry or None if not found
        """
        self._get_logger().debug("get_species_entry_for_assembly '{}'".format(assembly))
        if assembly in self._get_index_assembly():
            return self._get_index_assembly()[assembly]
        return None

    def count_ensembl_species(self):
        return len(self._get_species_data_dao())


if __name__ == '__main__':
    print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
