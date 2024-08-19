from collections import defaultdict
from tynt.core import Filter, _filter_generator

__all__ = ['filters']


class FilterSet:

    def __init__(self, filter_set, filter_names):
        self._filter_set = filter_set
        self._filter_names = filter_names
        self.filters = set()
        for filter_name in filter_names:
            fname = filter_name.replace('.', '_')
            setattr(
                self, fname,
                Filter(filter_set=filter_set, filter_name=filter_name)
            )
            self.filters.add(f'"{fname}"')

    def __repr__(self):
        return (f"<FilterSet: \"{self._filter_set}\"\n  "
                f"Available filters: [{', '.join(sorted(self.filters))}]>")


def assemble_filter_sets():

    filter_sets = defaultdict(list)

    for name in _filter_generator.available_filters():
        filter_set, filter_name = name.split('/')
        filter_sets[filter_set].append(filter_name)

    sets = dict()
    for filter_set, filter_names in filter_sets.items():
        sets[filter_set] = FilterSet(filter_set, filter_names)

    return sets


class DefaultFacilities:
    TWOMASS = None
    SLOAN = None
    Kepler = None
    TESS = None
    HST = None
    JWST = None
    LSST = None
    Keck = None
    WISE = None
    WFIRST = None
    Roman = None
    Spitzer = None
    GAIA = None
    CHEOPS = None


class FilterLookup(DefaultFacilities):
    _accessors = set()

    def __init__(self):

        filter_sets = assemble_filter_sets()

        for name, filter_set in filter_sets.items():
            if name[0].isnumeric():
                name = name.replace('2', 'TWO')

            setattr(self, name, filter_set)
            self._accessors.add(name)


filters = FilterLookup()
