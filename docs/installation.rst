.. include:: references.txt

************
Installation
************

Requirements
============

**astroplan** works on Linux, Mac OS X and Windows.
It requires Python 3.5+ as well as the following packages:

* `Numpy`_ (1.10 or later)
* `Astropy`_ (v1.3 or later)

For testing:

* `pytest-astropy`_

Installation
============

You can install the stable version of astroplan from PyPI with::

    pip install tynt

Alternatively, you can install the latest developer version of astroplan by
cloning the git repository::

    git clone https://github.com/bmorris3/tynt

...then installing the package with::

    cd tynt
    python setup.py install

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, start up python, and type::

    import tynt
    tynt.test()

If there are no errors, you are good to go!

.. note::
	If you want to run the tests that access the internet, you'll need to
	replace the last line above with ``tynt.test(remote_data=True)`` and
	have an active connection to the internet.