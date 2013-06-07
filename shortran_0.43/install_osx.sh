### script to assist in installing shortran dependencies on MacOSX

### WARNING! This script installs python 2.6 and selects this at the default python instance on your computer.
### If you prefer to keep the original MacOS supplied python, please install the dependencies manually instead.

### NOTE! Make sure that Macports and MySQL are installed before running this install script


### Set MySQL package name variable
mysql_version='mysql-5.5.24-osx10.6-x86_64'


### set exit on error
set -e

### make sure the .bash_profile file exists
touch ~/.bash_profile


### make sure you start from the shortran directory
ls shortran.py || exit "Please change to the shortran directory and run the install script again"


### add bowtie, and MySQL to the path
echo export PATH=$PATH:$PWD/bowtie-0.12.7/:$PWD/echo_v1_12/:/usr/local/"$mysql_version"/bin/ >> ~/.bash_profile
export PATH=$PATH:$PWD/bowtie-0.12.7/:$PWD/echo_v1_12/:/usr/local/"$mysql_version"/bin/



### install other dependencies
nice -n 19 sudo port install python26
nice -n 19 sudo port select --set python python26
nice -n 19 sudo port install py26-wxpython py26-numpy py26-matplotlib
nice -n 19 sudo port install py26-scipy py26-ipython
nice -n 19 sudo port install ncbi_tools

### start the mysql server
sudo mysqld_safe > server.log &


### make a basic mysql setup
### this will create a database called "profiles" and give a user called "user" full access to this database. No password is set.
sudo mysql < "mysql_setup.txt"

