import os

from tynt import FilterGenerator

output_path = os.path.join(os.path.dirname(__file__), 'docs', 'tynt', 'filters.rst')
d = FilterGenerator().available_filters(as_dataframe=True)
html = d.to_html()

out = """
Available Filters
-----------------

The table below lists the filters available by default by calling :py:meth:`~tynt.FilterGenerator.reconstruct`
using the included, light-weight filter archive in ``fft.fits``. All filters provided by the `SVO Filter Profile Service
<http://svo2.cab.inta-csic.es/theory/fps/>`_ are accessible via ``tynt`` with the
:py:meth:`~tynt.FilterGenerator.download_true_transmittance` method, which will access the internet and download
the transmittance curve at its "native" resolution, without the FFT approximation.

If the filters you need are not listed below but you'd like offline access to the ``tynt`` interface, you can
always create your own archive locally, with a custom set of filters, including all available filters. For a
tutorial on downloading your own custom archive, see :doc:`customarchive`.


To access a filter from the table below, combine the columns like this: ``col1/col2.col3``. For example, if you want
the 2MASS J filter in the first row below, you would call the following in ``tynt``:

.. code-block:: python

    from tynt import FilterGenerator

    f = FilterGenerator()
    filt = f.reconstruct('2MASS/2MASS.J')

.. raw:: html

{html}
"""

spaces = 4
indent_html = '\n'.join([spaces * " " + line for line in html.splitlines()])

with open(output_path, 'w') as w:
    w.write(out.format(html=indent_html))
