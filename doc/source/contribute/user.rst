Contribute as a user
####################

Users can contribute in a variety of ways, such as reporting bugs, requesting
new features, testing in-development features, starting discussions, answering
questions, and sharing their work with the community.

.. warning::

    Do not include any proprietary or sensitive information when reporting bugs
    or showcasing your work.

.. grid:: 1 2 3 3
    :padding: 2 2 2 2

    .. grid-item-card:: :fa:`bug` Report bugs
        :link: report-bugs
        :link-type: ref

        Found a bug? Report it here.

    .. grid-item-card:: :fa:`lightbulb` Request a new feature
        :padding: 2 2 2 2
        :link: request-a-new-feature
        :link-type: ref

        Got an idea for a new feature? Share it!

    .. grid-item-card:: :fa:`vial-circle-check` Test a new feature
        :padding: 2 2 2 2
        :link: test-a-new-feature
        :link-type: ref

        Anxious to try out a new feature? Here's how you can do it.

    .. grid-item-card:: :fa:`comments` Start a discussion
        :padding: 2 2 2 2
        :link: start-a-discussion
        :link-type: ref

        Want to discuss something? Start or contribute to a discussion.

    .. grid-item-card:: :fa:`comment-dots` Answer questions
        :padding: 2 2 2 2
        :link: answer-questions
        :link-type: ref

        Help others by answering their questions.

    .. grid-item-card:: :fa:`bullhorn` Share your work
        :padding: 2 2 2 2
        :link: share-your-work
        :link-type: ref

        Share your work with the community.


.. _report-bugs:

Report bugs
===========

If you encounter a bug or an issue while using the project, report it.
Your feedback helps to identify problems and get them resolved.

- Search the `PyAnsys Heart Issues`_ page to see if the issue has already been reported.

- Create an issue if one doesn't already exist.

  - Include a clear description of the issue.
  - Provide steps to reproduce the issue.
  - Mention the version of the project you're using.
  - Include screenshots or logs if possible.

.. _request-a-new-feature:

Request a new feature
=====================

Do you have an idea for a new feature or an improvement? Your suggestions are
welcome. You can request a new feature by creating an issue on the `PyAnsys Heart Issues`_
page.

.. _test-a-new-feature:

Test a new feature
==================

You can test a new feature before it is officially released. To do
so, you can install PyAnsys Heart from the source code by performing the
steps in the following child topics.

Clone the repository
--------------------

Clone and install the repository:

.. code-block:: bash

    git clone https://github.com/ansys/pyansys-heart

Install for users
-----------------

Install the latest version of PyAnsys Heart to test the latest features as
they are being developed, without having to wait for releases.

Set up a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to the project's root directory:

.. code-block::

    cd pyansys-heart

Create a new virtual environment named ``.venv`` to isolate your system's
Python environment:

.. code-block:: text

    python -m venv .venv

Activate this environment:

.. tab-set::

    .. tab-item:: Windows

        .. tab-set::

            .. tab-item:: CMD

                .. code-block:: text

                    .venv\Scripts\activate.bat

            .. tab-item:: PowerShell

                .. code-block:: text

                    .venv\Scripts\Activate.ps1

    .. tab-item:: macOS/Linux/UNIX

        .. code-block:: text

            source .venv/bin/activate

Install the latest version
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install PyAnsys Heart in editable mode:

.. code-block:: text

    python -m pip install .

Verify the installation by checking the version of the library:

.. code-block:: python

    from ansys.health.heart import __version__

    print(f"PyAnsys Heart version is {__version__}.")

.. jinja::

   .. code-block:: text

      >>> PyAnsys Heart version is {{ PYANSYS_HEART_VERSION }}.

.. _start-a-discussion:

Start a discussion
==================

Complex topics might require a discussion. Whether you want to know how to use
PyAnsys Heart for solving your specific problem or you have a suggestion for a new
feature, a discussion is a good place to start. You can open a new discussion
on the `PyAnsys Heart Discussions`_ page.

.. _answer-questions:

Answer questions
================

Another great way to contribute is to help others by answering their questions.
Maintain a positive and constructive attitude while answering questions. If you
don't know the answer, you can still help by pointing the person in the right
direction.

.. _share-your-work:

Share your work
===============

If you have used PyAnsys Heart to create something interesting, share it with the rest
of the community. You can share your work on the `PyAnsys Heart discussions`_ page. Include
a brief description of your work and any relevant links that others might find
useful.
