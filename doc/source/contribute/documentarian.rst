Contribute as a documentarian
#############################

.. grid:: 1 2 3 3
    :padding: 2 2 2 2

    .. grid-item-card:: :fa:`pencil` Write documentation
        :link: write-documentation
        :link-type: ref

        Learn how to get started, use, and contribute to the project.

    .. grid-item-card:: :fa:`laptop-code` Add a new example
        :link: write-examples
        :link-type: ref

        Write a new example to showcase the capabilities of PyAnsys Heart.

    .. grid-item-card:: :fa:`book` Build the documentation
        :link: build-documentation
        :link-type: ref

        Build the documentation to see your changes rendered.

.. _write-documentation:

Write documentation
===================

`Sphinx`_ is the tool used to generate PyAnsys Heart documentation. You write most of the content
in `ReStructuredText`_ files. However, some of the content, like the
`examples <../examples/index>`_, use a mix of `ReStructuredText`_ and Python files, thanks to `Sphinx-Gallery`_.
If you are interested in writing examples, see the :ref:`write-examples`.

The documentation is located in the ``doc/source`` directory. The landing page
is declared in the ``doc/source/index.rst`` file. The subdirectories contain
the pages of different sections of the documentation. Finally, the
``doc/source/_static/`` directory contains various assets like images and CSS
files.

The layout of the ``doc/source`` directory is reflected in the URLs of the
online documentation. For example, the
``doc/source/contribute/documentarian.rst`` file renders as the
``https://heart.health.docs.pyansys.com/version/stable/contribute/documentarian.html`` URL.

Thus, if you create a file, it is important to follow these rules:

- Use lowercase letters for file and directory names.
- Use short and descriptive names.
- Use hyphens to separate words.
- Logically organize the hierarchy of the files and directories

You must include all files in the table of contents. Sphinx does not permit any orphan files.
If you do not include a file in the table of contents, Sphinx raises a warning that causes
the build to fail.

You declare the table of contents using a directive like this:

.. code-block:: rst

    .. toctree::
        :hidden:
        :maxdepth: 3

        path-to-file-A
        path-to-file-B
        path-to-file-C
        ...

The path to the file is relative to the directory where the table of contents
is declared.

.. _write-examples:

Write a new example
===================

The `examples <../examples/index>`_ section of the documentation showcases different
capabilities of PyAnsys Heart. Each example is a standalone Python script. You group
related examples into subdirectories. Despite being PY files, they are written in a mix
of `ReStructuredText`_ and Python. This is possible thanks to the `Sphinx-Gallery`_
extension.

Documentarians writing new examples are encouraged to familiarize themselves with
`Structuring Python scripts for Sphinx-Gallery <https://sphinx-gallery.github.io/stable/syntax.html>`_.
Once the PY file for a new example is properly set up, Sphinx-Gallery automatically
generates `Sphinx`_ `ReStructuredText`_  (RST) files from it. The rendering of the resulting
RST file for each example provides links for downloading a IPYNB (Jupyter notebook) and PY file.

Finally, here are some tips for writing examples:

- Start the example with an explanation of the main topic. Try to use as many relevant
  keywords as possible in this section for search engine optimization.

- Include an explanation with each code cell. The explanations should
  be included before, not after, the corresponding code.

- The examples are built with the documentation. During the build process,
  screenshots are inserted in the rendered document. You do not need
  to include the screenshots yourself.

- When creating a child directory that is to include multiple related examples, ensure that
  you include a ``README.txt`` file  with the ReStructuredText content to
  use for the index page for this subsection's examples in the generated documentation.

.. _build-documentation:

Build the documentation
=======================

`Tox`_ is used for automating the build of the documentation.

To install Tox:

.. code-block:: text

    python -m pip install tox

There are different environments for cleaning the build, building the documentation
in different formats such as HTML and PDF, and running the tests.

The following environments are available:

.. jinja:: toxenvs

    .. dropdown:: Documentation environments
        :animate: fade-in
        :icon: three-bars

        .. list-table::
            :header-rows: 1
            :widths: auto

            * - Environment
              - Description
              - Command
            {% for environment in envs %}
            {% set name, description  = environment.split("->") %}
            {% if name.startswith("doc-")%}
            * - {{ name }}
              - {{ description }}
              - python -m tox -e {{ name }}
            {% endif %}
            {% endfor %}
