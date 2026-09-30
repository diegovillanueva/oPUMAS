PUMAS
======
Parameterization for Unified Microphysics Across Scales

This repository contains the open source code for most versions of the Morrison-Gettleman (MG) microphysics as well as the most recent releases called the Parameterization for Unified Microphysics Across Scales, or PUMAS.

 Checking out and running PUMAS
================================
 Running with an ESCOMP/CAM branch
-----------------------------------
The code in this repository is not sufficient for completing a simulation on its own. To run PUMAS within the Community Atmosphere Model (CAM), you will need to check out a branch of CAM.

Checking out any CAM branch since cam6_3_046 will give you access to the PUMAS microphysics scheme using the cam_dev or cam7 physics package.

For more information on checking out and running CAM, please see:
https://github.com/ESCOMP/CAM

 Modifying PUMAS code
-----------------------------------
PUMAS and CAM are community supported scientific projects. We welcome feedback in the Issues section of our Github repository and Pull Request from contributing community members. For information on the Github workflow for NCAR models and parameterizations, read the wiki description here:
https://github.com/ESCOMP/CAM/wiki

Because the PUMAS microphysics parameterization is managed in a seperate repository from the NCAR CAM model, any development will need two Github forks. You will need (1) A fork of https://github.com/ESCOMP/CAM and (2) A fork of https://github.com/ESCOMP/PUMAS

These are instructions for checking out a branch of CAM and pointing to a differnt branch of PUMAS then what you recieve from a checked out CAM branch or tag.

Clone a repository that contains a CAM source tree. This will likely be your fork, but could be a group fork like the PUMASDevelopment fork. ::
```
git clone https://github.com/PUMASDevelopment/CAM.git Github_CAM_PUMAS_Clone
```
This will create a directory Github_CAM_PUMAS_Clone in your current working directory.

Go into the newly created directory and checkout a modified branch. ::
```
cd Github_CAM_PUMAS_Clone git checkout username/my_new_branch
```
From the root of the CAM clone, run the script manage_externals/checkout_externals. ::
```
./manage_externals/checkout_externals
```
The checkout_externals script will populate the cam directory with the relevant versions of each of the components along with the CIME infrastructure code.

At this point you have all of the code needed for CAM with the PUMAS microphysics available.

To make changes to PUMAS, first create a branch for your work in the ESCOMP/PUMAS Github repository (or your fork of this repository) by clicking on the "Branch:master" drop down box on the middle left part of the main page (just below the purble line), and type the name of your new branch into the "Find or create a branch..." text area.

For new branch names, it is generally a good idea to put your Github name first, and then the goal of the branch after a slash. So, a name for a branch to improve graupel tendencies might be "katetc/graupel_upgrade". Hit enter, and your new branch is now shown in Github.

The second step for making changes in PUMAS is to update the source code in the pumas subdirectory to work with this branch. The easiest way to manage this is by using manage_externals. Edit the `Externals_CAM.cfg` file to point to your fork and branch under the `[pumas]` section. Then run `manage_externals/checkout_externals` from the root of your CAM checkout. You can go into your PUMAS branch to make sure that it is checked out correctly. ::
```
cd Github_CAM_PUMAS_Clone/src/physics/pumas
git status
```
Once you have your own branch checked out, you can make local changes to the code. When it's time to commit them, you will need to ::
```
git add filename.F90
git commit -m "Commit message"
git push
```
This will push your changes to your remote branch. When you have finished with ALL of your changes (can be multiple commits), then it will be time to merge your branch changes back to the master branch in the ESCOMP/PUMAS repo. You can do this by clicking "Pull requests" in the PUMAS github repo (or fork) main page, and then the big green "New Pull Request" button. In the gray bar, for the "compare" pull down, click and select your branch. This will show all of the differences between your branch and the main_cam PUMAS branch. If you are happy with these changes, click the green "Create pull request" button again, and that will take you to a pull request form. Fill out a description of your changes and issue the pull request. This request will be reviewed by a software engineer and then merged into the main branch. Your development branch can then be deleted, and you can start a new one for the next issue.


References
===========
Sun, J., Dennis, J. M., Mickelson, S. A., Vanderwende, B., Gettelman, A., & Thayer-Calder, K. (2023). Acceleration of the Parameterization of Unified Microphysics Across Scales (PUMAS) on the graphics processing unit (GPU) with directive-based methods. Journal of Advances in Modeling Earth Systems, 15, e2022MS003515. https://doi.org/10.1029/2022MS003515

Gettelman, A., Morrison, H., Eidhammer, T., Thayer-Calder, K., Sun, J., Forbes, R., McGraw, Z., Zhu, J., Storelvmo, T., and Dennis, J.: Importance of ice nucleation and precipitation on climate with the Parameterization of Unified Microphysics Across Scales version 1 (PUMASv1), Geosci. Model Dev., 16, 1735–1754, https://doi.org/10.5194/gmd-16-1735-2023, 2023.

Gettelman, A., H. Morrison, K. Thayer‐Calder, and C. M. Zarzycki. 2019. The Impact of Rimed Ice Hydrometeors on Global and Regional Climate. Journal of Advances in Modeling Earth Systems. https://doi.org/10.1029/2018MS001488.

Gettelman, A. and H. Morrison, Advanced Two-Moment Microphysics for Global Models. Part I: Off line tests and comparisons with other schemes. J. Climate, 28, 1268-1287. doi: 10.1175/JCLI-D-14-00102.1, 2015.

Gettelman, A., H. Morrison, S. Santos, P. Bogenschutz and P. H. Caldwell. Advanced Two-Moment Microphysics for Global Models. Part II: Global model solutions and Aerosol-Cloud Interactions. J. Climate, 28, 1288-1307. doi:10.1175/JCLI-D-14-00103.1 , 2015.
