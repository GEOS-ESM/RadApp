# RadApp Fixture

## How to build RadApp

### Preliminary Steps

#### Load Build Modules

In your `.bashrc` or `.tcshrc` or other rc file add a line:

##### NCCS
```
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES12
```

##### NAS
```
module use -a /nobackup/gmao_SIteam/modulefiles
```

##### GMAO Desktops
On the GMAO desktops, the SI Team modulefiles should automatically be
part of running `module avail` but if not, they are in:

```
module use -a /ford1/share/gmao_SIteam/modulefiles
```

Also do this in any interactive window you have. This allows you to get module files needed to correctly checkout and build the model.

Now load the `GEOSenv` module:
```
module load GEOSenv
```
which obtains the latest `git`, `CMake`, and `mepo` modules.

#### Obtain the Model

On GitHub, there are three ways to clone the model: SSH, HTTPS, or GitHub CLI.
The first two are "git protocols" which determine how `git` communicates with
GitHub: either through https or ssh. (The latter is a CLI that uses either ssh or
https protocol underneath.)

For developers of GEOS, the SSH git protocol is recommended as it can avoid some issues if
[two-factor authentication
(2FA)](https://docs.github.com/en/github/authenticating-to-github/securing-your-account-with-two-factor-authentication-2fa)
is enabled on GitHub.

##### SSH

To clone the RadApp using the SSH url (starts with `git@github.com`), you run:
```
git clone git@github.com:GEOS-ESM/RadApp.git
```

###### Permission denied (publickey)

If this is your first time using GitHub with any SSH URL, you might get this
error:
```
Permission denied (publickey).
fatal: Could not read from remote repository.

Please make sure you have the correct access rights
and the repository exists.
```

If you do see this, you need to [upload an ssh
key](https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account)
to your GitHub account. This needs to be done on any machine that you want to
use the SSH URL through.


##### HTTPS

To clone the model through HTTPS you run:

```
git clone https://github.com/GEOS-ESM/RadApp.git
```

Note that if you use the HTTPS URL and have 2FA set up on GitHub, you will need
to use [personal access
tokens](https://docs.github.com/en/github/authenticating-to-github/accessing-github-using-two-factor-authentication#authenticating-on-the-command-line-using-https)
as a password.

##### GitHub CLI

You can also use the [GitHub CLI](https://cli.github.com/) with:
```
gh repo clone GEOS-ESM/RadApp
```

Note that when you first use `gh`, it will ask what your preferred git protocol
is (https or ssh) to use "underneath". The caveats above will apply to whichever
you choose.

---

### Single Step Building of the Model

If all you wish is to build the model, you can run `parallel_build.csh` from a head node. Doing so will checkout all the external repositories of the model and build it. When done, the resulting model build will be found in `build/` and the installation will be found in `install/` with setup scripts like `gcm_setup` and `fvsetup` in `install/bin`.

#### Develop Version of RadApp

`parallel_build.csh` provides a special flag for checking out the
development branches of GEOSgcm_GridComp and GMAO_Shared. If you run:

```
parallel_build.csh -develop
```
then `mepo` will run:

```
mepo develop GEOSgcm_GridComp GMAO_Shared
```

#### Debug Version of RadApp

To obtain a debug version, you can run `parallel_build.csh -debug` which will build with debugging flags. This will build in `build-Debug/` and install into `install-Debug/`.

#### Do not create and install source tarfile with parallel_build

Note that running with `parallel_build.csh` will create and install a tarfile of the source code at build time. If you wish to avoid
this, run `parallel_build.csh` with the `-no-tar` option.

---

### Multiple Steps for Building the Model

The steps detailed below are essentially those that `parallel_build.csh` performs for you. Either method should yield identical builds.

#### Mepo

The RadApp is comprised of a set of sub-repositories. These are
managed by a tool called [mepo](https://github.com/GEOS-ESM/mepo). To
clone all the sub-repos, you can run `mepo clone` inside the fixture:

```
cd RadApp
mepo clone
```

The first command initializes the multi-repository and the second one
clones and assembles all the sub-repositories according to
`components.yaml`

#### Checking out develop branches of GEOSgcm_GridComp  and GMAO_Shared

To get development branches of GEOSgcm_GridComp  and GMAO_Shared (a la
the `-develop` flag for `parallel_build.csh`, one needs to run the
equivalent `mepo` command. As mepo itself knows (via `components.yaml`) what the development branch of each
subrepository is, the equivalent of `-develop` for `mepo` is to
checkout the development branches of GEOSgcm_GridComp and GMAO_Shared:
```
mepo develop GEOSgcm_GridComp GMAO_Shared
```

This must be done *after* `mepo clone` as it is running a git command in
each sub-repository.

#### Build the Model

##### Load Compiler, MPI Stack, and Baselibs
On tcsh:
```
source @env/g5_modules
```
or on bash:
```
source @env/g5_modules.sh
```

##### Create Build Directory
We currently do not allow in-source builds of RadApp. So we must make a directory:
```
mkdir build
```
The advantages of this is that you can build both a Debug and Release version with the same clone if desired.

##### Run CMake
CMake generates the Makefiles needed to build the model.
```
cd build
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_INSTALL_PREFIX=../install
```
This will install to a directory parallel to your `build` directory. If you prefer to install elsewhere change the path in:
```
-DCMAKE_INSTALL_PREFIX=<path>
```
and CMake will install there.

###### Create and install source tarfile

Note that running with `parallel_build.csh` will create and install a tarfile of the source code at build time. But if CMake is run by hand, this is not the default action (as many who build with CMake by hand are developers and not often running experiments). In order to enable this at install time, add:
```
-DINSTALL_SOURCE_TARFILE=ON
```
to your CMake command.

##### Build and Install with Make
```
make -jN install
```
where `N` is the number of parallel processes. On discover head nodes, this should only be as high as 2 due to limits on the head nodes. On a compute node, you can set `N` has high as you like, though 8-12 is about the limit of parallelism in our model's make system.

### Run the GCM

Once the model has built successfully, you will have an `install/` directory in your checkout. To run `gcm_setup` go to the `install/bin/` directory and run it there:
```
cd install/bin
./gcm_setup
```

## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3](LICENSE-NOSA).
