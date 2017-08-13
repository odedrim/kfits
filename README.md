# kfits
Kinetics Fitting Software

This software package was written by Oded Rimon, at Dr. Dana Reichmann's Lab, HUJI, Israel.

The Lab is Located at:  
Dept. of Biological Chemistry, room 1-626  
The Alexander Silberman Institute of Life Science  
The Hebrew University of Jerusalem  
Givat Ram Campus, Jerusalem, 91904, Israel  
Tel: (+972) 2-658-5703

Email us at:  
Oded: oded.rimon AT mail.huji.ac.il  
Dana: danare AT mail.huji.ac.il

---

# Usage

Once you have *Kfits* installed on your system, simply type `kfits` in a command line / terminal window (on Windows: Command Prompt or the Windows Run command).  
To act as a server for other computers on the local network, use `kfits_server` instead. First, edit the file &lt;Kfits Directory&gt;/afgui/afgui/settings.py and add your IP address(es) to ALLOWED_HOSTS. Then run `kfits_server` and type in the IP interface you would like to serve on (0.0.0.0 for all interfaces) and the port to serve on (8000 is usually fine) as prompted.  
If you use `kfits` (or `kfits_server` with port 8000), *Kfits* will be accessible *via* the URL `http://127.0.0.1:8000/fitter/`. In either case, a browser window should open at the correct address, and the graphical user interface will guide you through the workflow of *Kfits*. The 'Instructions' button on the top of the screen will lead you to a detailed description of the workflow, as well as instructions on how to extend the capabilities of the software.

---

# Installation Instructions

***Executive Summary*** First of all try: `pip install kfits`  
If that failed, follow the instructions below.

## 1 Get Python 2.7

If you already have Python 2 (2.6 or later, **but not 3.x**) installed, skip to step 2. If you have both Python 2 and Python 3 installed, that's fine too. If you're working on a linux system, chances are you already have Python 2.7 installed.  
Otherwise, follow the instructions on [this official download page](https://www.python.org/downloads/) to download and install Python 2.7. Note that Python 3 will not work for this software, so don't be tempted ;)  
For Windows and MacOS, you may choose to work with [ActiveState ActivePython](https://www.activestate.com/activepython/downloads) or another distribution. As long as it's python 2.7, it's fine.

## 2 Install pip - the python package installer
If you are certain you have pip installed, you may skip to step 3.
Otherwise, follow the instructions on the [pip installation page](https://pip.pypa.io/en/stable/installing/).  
Shortly, depending on your system, it should be as easy as running (on Fedora <22):
```
sudo yum upgrade python-setuptools
sudo yum install python-pip python-wheel
```
Or (on Fedora >=22):
```
sudo dnf upgrade python-setuptools
sudo dnf install python-pip python-wheel
````
Or (on Debian / Ubuntu):
```
sudo apt-get install python-pip
```
Or (on Windows or any other system):
Download [get_pip.py](https://bootstrap.pypa.io/get-pip.py) and run `python get_pip.py` in the folder you downloaded it to.

## 3 Get the numerical (NumPy) and scientific (SciPy) python packages
Depending on your specific system, environment and sheer luck, this step may be smooth and easy or long and tedious. We sincerely apologise if the latter applies to you - unfortunately it is out of our hands (we even tried to rewrite the code without any use of scipy, but it is just so darn useful..!).  
First, check if you already have scipy. Open a python terminal (either IDLE or by simply typing `python` in a command line window) and try: `import scipy`. If no error message was presented - you're in luck! Skip to step 4!  
No luck..?

### 3.1 Install NumPy
Well, this is where things start to get a bit complicated.  
There are many ways to install NumPy, and frankly I haven't seen one that works on all systems. If you're working on a **linux or linux-based system**, the best way IMHO is using `pip install numpy`. However, you will first need to install `python-dev`, `liblapack-dev` and `gfortran`, and if you don't have a suitable C/C++ compiler installed, you would need one of those as well. For me, on an **Ubuntu** system, the following two commands were sufficient:
```
sudo apt-get install python-dev liblapack-dev gfortran
sudo pip install numpy
```
However, your results may vary (and if you're using another linux distribution, remember to replace `apt-get install` with the appropriate installation command, e.g. `yum install`).  
On a **Windows** system, you may prefer to download pre-compiled NumPy from [this University of California page](http://www.lfd.uci.edu/~gohlke/pythonlibs/) and follow the instructions there; or you may choose to install [Visual C++ for Python 2.7](http://www.microsoft.com/en-us/download/details.aspx?id=44266) and use `pip install numpy` as you would on linux-based systems.  
Once you chose your method and you think you have NumPy installed, open a python terminal and type `import numpy`. If no error message appears - you're finally done with this step, and almost ready to install _Kfits_ :)

### 3.2 Install SciPy
Once you successfuly installed numpy, installing scipy should be a piece of cake. Whatever was your last step for installing numpy (e.g. `pip install numpy` or downloading from [the unofficial binaries page](http://www.lfd.uci.edu/~gohlke/pythonlibs/)) - do exactly the same for scipy. On linux-based systems, `sudo apt-get install python-scipy` may be your best bet.  
Once you have scipy successfully installed, you should be able to type in `import scipy` in a python terminal and receive no error messages when you press ENTER.

## 4 Download and Install *Kfits*
Finally, everything is prepared and you can download and install *Kfits*.  
There is more than one way to do this, and you may choose whichever you prefer. We think the first is the simplest and most straight-forward method, but really, feel free to use others. We won't be offended. Really.

### *Via* pip
```
pip install kfits
```
Seriously. That's it.

### From source
You can download the source code of *Kfits* from [github](https://github.com/odedrim/kfits). Either clone it with git (go to a directory of your choice, type in `git clone https://github.com/odedrim/kfits`, and you'll have a subdirectory named `kfits` containing the source code) or download as ZIP and unzip it to a directory of your choice.  
Once you have a directory which contains the source code of kfits, change dir (`cd`) to that directory and type `python setup.py install` or `sudo python setup.py install` (to install with root privileges on a linux-based system)... And you're done.

### Precompiled Packages (binary distributions)
We will do our best to precompile binary packages for as many operating systems and environments as we can - but eventually we are just Biochemists in the academy, we don't really have access to that many different systems.  
The binary distributions will be made available at [github](https://github.com/odedrim/kfits), in the [install directory](https://github.com/odedrim/kfits/tree/master/install). To install them, simply download the one best suited to your setup, and double click / run it to install.

