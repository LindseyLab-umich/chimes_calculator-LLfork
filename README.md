SPECIAL FOR LL DEVELOPERS
==========================

Enabling clone to LC computers from GitLab
------------------------------------------

The GitLab repo lives [here](https://lc.llnl.gov/gitlab/chimes/chimes_calculator/).

To access the repo, Becky must add you. She will need your GitLab username, which can be found by clicking your user icon at the top right of the page. It will start with an `@`, e.g., `@rlindsey`.

Once you get repo access, you need to configure your SSH/SSH keys so you can clone the repo. To do so, check if you have a public key (`ls ~/.ssh/*pub`). If you don't have one, generate one. If you do have one, make sure your `~/.ssh/config` tells GitLab where to look for it, e.g. add something like:

Host czgitlab.llnl.gov\
&nbsp;&nbsp;AddKeysToAgent yes\
&nbsp;&nbsp;IdentityFile ~/.ssh/id_ed25519\
&nbsp;&nbsp;PubkeyAcceptedKeyTypes +ssh-ed25519

Next, add it to your SSH keys on Gitlab. To do so, enter "ssh" in the search bar. From there, it is just like Github.

Enabling push/pull between LC and Git**Hub**
---------------------------------------------

First, ensure you have an ssh key for the LC machine defined in GitHub. Next, configure the SSO. To do this, on the GitHub SSH key page, fine your LC key, click the "Configure SSO" button to the left, and click "Approve". Once you do this, you will be able to push/pull.


Notes **for Becky**: Push/pull: LC <---> Git**Lab** vs LC <---> Git**Hub**
---------------------------------------------

**This bit is already done:**


To begin, clone the Git**Lab** repo to the LC with:

`git clone ssh://git@czgitlab.llnl.gov:7999/chimes/chimes_calculator.git`

Verify the branches, and check out the ones that matter (main and develop for this repo):

`git branch -r `\
`git checkout develop`\
`git checkout main`

Next, create a new empty repo on Git**Hub**. Github will provide you with the URL for the repo. We'll use this URL to set the GitHub repo as another remote:

`git remote add LLgithub git@github.com:LindseyLab-umich/chimes_calculator-260324-LLNL_Merge.git`

Confirm it worked with `git remote -v`. You should see two sets of fetch/push URLS.

Now we push the develop and main branches into this Git**Hub** remote. Remember, we're in main now:

`git push LLgithub main`\
`git checkout develop`\
`git push LLgithub develop`

Now we can operate on the repo in GitHub as usual, forking, etc. 

Eventually, we'll need to communicate the edits in LLgithub back to the Git**Lab** repo. Assuming you're on the branch you want to push/pull from, the way to do this will be, e.g.:

1. On the LC, pull changes from LLgithub, e.g., `git pull LLgithub`
2. On the LC, push local changes (just pulled from LLgithub) to GitLab, e.g., `git push origin`



*********************************************************************************************************

<p style="text-align:center;">
    <img src="./doc/ChIMES_Github_logo-2.png" alt="" width="250"/>
</p>
<hr>

The Chebyshev Interaction Model for Efficient Simulation (ChIMES) is a machine-learned interatomic potential that can target chemical reactivity. ChIMES models are able to approach quantum-accuracy through a systematically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. ChIMES has successfully been applied to a number of condensed phase systems, including water under ambient and extreme conditions, molten carbon, and liquid carbon monoxide under planetary interior conditions. ChIMES can also be used as a many-body repulsive energy for the density functional based tight binding (DFTB) method.

The ChIMES calculator comprises a flexible tool set for evaluating ChIMES interactions (e.g. in simulations, single point calculations, etc). Users have the option of directly embedding the ChIMES calculator within their codes (e.g. see ‘’The ChIMES Calculator,’’ in the documentation for advanced users), or evaluating interactions through the beginner-friendly serial interface, each of which have Python, C++, C, and FORTRAN API’s.


<hr>

Documentation
----------------

[**Full documentation**](https://chimes-calculator.readthedocs.io/en/latest/) is available.

<hr>

Community
------------------------

Questions, discussion, and contributions (e.g. bug fixes, documentation, and extensions) are welcome. 

Additional Resources: [ChIMES Google group](https://groups.google.com/g/chimes_software).

<hr>

Contributing
------------------------

Contributions to the ChIMES calculator should be made through a pull request, with ``develop`` as the destination branch. A test suite log file should be attached to the PR. For additional contributing guidelines, see the [documentation](https://chimes-calculator.readthedocs.io/en/latest/contributing.html).

The ChIMES calculator `develop` branch has the latest contributions. Pull requests should target `develop`, and users who want the latest package versions,
features, etc. can use `develop`.

<hr>

Releases
--------

For most users, we recommend using the ChIMES calculator [stable releases](https://github.com/rk-lindsey/chimes_calculator/releases).

Each ChIMES calculator release series also has a corresponding branch, e.g. `releases/v0.14` has `0.14.x` versions, and `releases/v0.13` has `0.13.x` versions. We back-port important bug fixes to these branches but we do not advance the package versions or make changes that would otherwise change the way ChIMES calculator is used or behaves. So, you can base your ChIMES deployment on a release branch and `git pull` to get fixes, without the continuous changes that comes with `develop`.  The latest release is always available with the `releases/latest` tag.

<hr>

Authors
----------------

The ChIMES calculator was developed by Rebecca K. Lindsey, Nir Goldman, and Laurence E Fried.

Contributors can be found [here](https://github.com/rk-lindsey/chimes_calculator/graphs/contributors).

<hr>

Citing
----------------

See [the documentation](https://chimes-calculator.readthedocs.io/en/latest/citing.html) for guidance on referencing ChIMES and the ChIMES calculator in a publication.

<hr>

License
----------------

The ChIMES calculator is distributed under terms of [LGPL v3.0 License](https://github.com/rk-lindsey/chimes_calculator/blob/main/LICENSE).

LLNL-CODE- 817533
