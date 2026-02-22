# Contributing to AxiSEM3D

AxiSEM3D is a community project that lives by the participation of its
members. That includes you! Our goal is to build an inclusive
and participatory community. We are happy that you are interested in
participating!

## Asking and answering questions about AxiSEM3D

For questions about AxiSEM3D on all levels, please use the
[AxiSEM3D forum](https://community.geodynamics.org/c/axisem).

Please note that everyone answering forum posts is volunteering
their time to help make our community a friendly and helpful
place. Therefore, please keep your inquiries
polite, take some effort to make them easily readable, and
include all useful information as described below. If there is
something you can test on your own, please test it first,
before asking a question on the forum.

There is no guarantee that we can or will help with your problem, but
you are welcome to contact us. Depending on the level of
effort provided a thank you note or an acknowledgment in
a paper is always appreciated. Keep in mind that for significant
and scientifically creative contributions by a community member
(but usually only then), a co-authorship on a publication
is appropriate.

## Reporting bugs and asking questions

Help the community by reporting any bugs that you
may find. We keep track of all open issues related to AxiSEM3D
[here](https://github.com/AxiSEMunity/AxiSEM3D/issues).

Please follow these instructions before opening a new bug report or ask
a question:


- Search in the
  [list of open and closed issues](https://github.com/AxiSEMunity/AxiSEM3D/issues?q=is%3Aissue)
  for a duplicate of your question.
- Search in the [AxiSEM3D forum](https://community.geodynamics.org/c/axisem) for
  a duplicate of your question.
- If you did not find an answer in the previous searches, open a new question:
  - If you suspect you have found a bug, open a new
    [issue](https://github.com/AxiSEMunity/AxiSEM3D/issues/new) and explain your
    problem as described below.
  - If you are not sure how to set up a model, ask a question in the
    [AxiSEM3D forum](https://community.geodynamics.org/c/axisem) as described below.
  - If you are not sure what to do, you can post a question in the
    [AxiSEM3D forum](https://community.geodynamics.org/c/axisem).
- In either case, attach the following information:
  - your input files with a simplified and small model that reproduces the
    issue,
  - one or several screenshots of the full error message you saw on your
    screen,
  - any information that helps us understand why you think this is a bug
    (screenshots, data series, comparisons to reference results, etc.),
  - instructions for how to reproduce the problem.

Without providing the information above, we will be less likely able to help
you, it may take significantly longer until you receive a reply, and we will
just ask you for this information anyway.

## Making AxiSEM3D better

AxiSEM3D is a community project. We encourage contributions of all kinds. 
Much
appreciated contributions include new examples (cookbooks, tests, or benchmarks),
extended documentation (every paragraph helps), and in particular fixing typos
or updating outdated documentation. We also encourage
contributions to the core functionality. If you consider making a
larger contribution to the core functionality, please open a new
[issue](https://github.com/AxiSEMunity/AxiSEM3D/issues/new) first, to discuss
your idea with one of the Principal Developers. This allows us to give you early
feedback and prevents you from spending too much time on a project that might already be
planned or that conflicts with other plans.

### Getting started with git and GitHub

If you are new to using git or the GitHub platform, you may find it
helpful to review [video lecture
32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html), which
should be enough to help you get started with using GitHub or possibly
contributing to AxiSEM3D itself. Alternatively, GitHub provides a helpful
guide on the process of contributing to an open-source project
[here](https://opensource.guide/how-to-contribute/).

### Opening pull requests

To make a change to AxiSEM3D you should:

- Create a
[fork](https://guides.github.com/activities/forking) (through GitHub) of
the code base.
- Create a separate
[branch](https://guides.github.com/introduction/flow/) (sometimes called a
feature branch) on which you do your modifications.
- You can propose that your branch be merged into the AxiSEM3D
code by opening a [pull request](https://guides.github.com/introduction/flow/).
This will give others a chance to review your code.

We follow the philosophy that no pull request (independent of the author) is
merged without a review from one other member of the community, and approval of
one of the maintainers. This applies to maintainers as well as to first-time
contributors. We know that a review can be a daunting process, but pledge to
keep all comments friendly and supportive, as you can see in the list of 
[open
pull requests](https://github.com/AxiSEMunity/AxiSEM3D/pulls)! We are as
interested in making AxiSEM3D better as you are!

While this seems very
formal, keeping all of the code review in one place makes it easier to
coordinate changes to the code (and there are usually several people making
changes to the code at once). This process is described at length in the
deal.II video [lecture
32.8](http://www.math.colostate.edu/~bangerth/videos.676.32.8.html).  Please do
not hesitate to ask questions about the workflow on the forum if you are
not sure what to do.

### Coding conventions


We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to keep the coding style consistent. A `.clang-format` files is included in the repository and all source code can be
automatically formatted by running

### Changelog entries
If your new pull request creates a change that is noticeable to AxiSEM3D users,
please add an entry to our
[changelog](https://github.com/AxiSEMunity/AxiSEM3D/blob/master/CHANGELOG.md). 
Start your newline with an action verb such as added, implemened, changed, fixed, removed, etc., explain the nature of your change,
include your name and date in the file as shown (this ensures you will get credit for your work), and lastly, include the pull request number.

## Acknowledgment of contributions

The AxiSEM3D community is grateful for every contribution! 
Your contribution is *formally* acknowledged in several ways:

- Every commit that is merged into the AxiSEM3D repository makes you part of
  the growing group of AxiSEM3D
  [contributors](https://github.com/AxiSEMunity/AxiSEM3D/graphs/contributors).
- For every release the most significant entries of our
  [changelog](https://github.com/AxiSEMunity/AxiSEM3D/blob/master/CHANGELOG.md)
  are selected to generate our release announcements. Additionally, all entries
  remain available for all previous releases of AxiSEM3D inside this file.
- If you contributed a significant part of the manual (such as a new cookbook,
  benchmark, or subsection), you will be acknowledged in [Authors.md](https://github.com/AxiSEMunity/AxiSEM3D/blob/master/AUTHORS.md).
- The Principal Developers of AxiSEM3D meet regularly to discuss the direction
  of the project and to consider inviting new members to join the group of
  Principal Developers.
  New members are typically invited based on the following criteria:
  - A deep understanding of AxiSEM3D's structure, design, and long-term vision;
  - A demonstrated commitment to advancing the project's goals and supporting
    the AxiSEM3D community;
  - Meaningful contributions to AxiSEM3D. These contributions may of course be
    to the source code, but can also include documentation, tutorials,
    benchmarks, or active engagement on the forum;
  - Consistent and active involvement in the project for over a year,
    beyond participation in user meetings;
  - A collaborative mindset and a constructive, solution-oriented approach
    to discussions.

  The group of current Principal Developers is listed in the
  [AUTHORS.md](https://github.com/AxiSEMunity/AxiSEM3D/blob/master/AUTHORS.md)
  file in the main repository.

## License

AxiSEM3D is published under the
[MIT license](https://github.com/AxiSEMunity/AxiSEM3D/blob/master/LICENSE).
While you will retain copyright on your contributions, all changes to the code
must be provided under this common license.
